/* Ensembl REST API client: throttled (Ensembl allows ~15 req/s), cached, retries on 429/5xx/network errors. */
const API = (() => {
  const BASE = 'https://rest.ensembl.org';
  const MAX_CONCURRENT = 5;
  const cache = new Map();
  const queue = [];
  let active = 0;

  const sleep = ms => new Promise(r => setTimeout(r, ms));
  // Ensembl has transient 5xx outages lasting tens of seconds; back off 1.5, 3, 6, 12, 15, 15 s (~50 s total).
  const backoff = attempt => Math.min(1500 * 2 ** attempt, 15000);

  function schedule(fn) {
    return new Promise((resolve, reject) => { queue.push({ fn, resolve, reject }); pump(); });
  }
  function pump() {
    while (active < MAX_CONCURRENT && queue.length) {
      const { fn, resolve, reject } = queue.shift();
      active++;
      fn().then(resolve, reject).finally(() => { active--; pump(); });
    }
  }

  function request(path, { method = 'GET', body = null, retries = 6 } = {}) {
    const key = method + ' ' + path + (body ? JSON.stringify(body) : '');
    if (cache.has(key)) return cache.get(key);
    const p = schedule(async () => {
      for (let attempt = 0; ; attempt++) {
        let url = BASE + path;
        const init = { method, headers: { Accept: 'application/json' } };
        if (method === 'GET') {
          url += (url.includes('?') ? '&' : '?') + 'content-type=application/json';
        } else {
          // text/plain keeps this a CORS "simple" request (no OPTIONS preflight, which Ensembl often
          // answers with a 5xx lacking CORS headers); Ensembl still parses the JSON body.
          init.headers['Content-Type'] = 'text/plain';
          init.body = JSON.stringify(body);
        }
        let r;
        try { r = await fetch(url, init); }
        catch (e) {
          if (attempt < retries) { await sleep(backoff(attempt)); continue; }
          throw new Error('Network error contacting Ensembl (' + e.message + ')');
        }
        if ((r.status === 429 || r.status >= 500) && attempt < retries) {
          const ra = parseFloat(r.headers.get('Retry-After'));
          await sleep(isFinite(ra) ? ra * 1000 + 100 : backoff(attempt));
          continue;
        }
        if (!r.ok) {
          let msg = r.status + ' ' + r.statusText;
          try { const j = await r.json(); if (j.error) msg = j.error; } catch (_) {}
          const err = new Error(msg); err.status = r.status; throw err;
        }
        return r.json();
      }
    });
    cache.set(key, p);
    p.catch(() => cache.delete(key));
    return p;
  }

  const enc = encodeURIComponent;
  const isEnsGene = q => /^ENS[A-Z]*G\d{6,}(\.\d+)?$/i.test(q);
  const stripVersion = id => id.replace(/\.\d+$/, '');

  async function lookupGene(species, query) {
    query = query.trim();
    const opts = 'expand=1&utr=1&mane=1';
    if (isEnsGene(query)) return request(`/lookup/id/${enc(stripVersion(query))}?${opts}`);
    try {
      return await request(`/lookup/symbol/${species}/${enc(query)}?${opts}`);
    } catch (e) {
      // fall back to synonym / xref search
      const hits = await request(`/xrefs/symbol/${species}/${enc(query)}?object_type=gene`).catch(() => []);
      const g = (hits || []).find(h => h.type === 'gene' && /^ENS/.test(h.id));
      if (!g) throw new Error(`No ${species.replace('_', ' ')} gene found for "${query}"`);
      return request(`/lookup/id/${enc(g.id)}?${opts}`);
    }
  }

  return {
    lookupGene,
    overlapTranscripts: geneId => request(`/overlap/id/${enc(geneId)}?feature=transcript`),
    overlapGenes: (species, chr, s, e) => request(`/overlap/region/${species}/${chr}:${Math.max(1, Math.round(s))}-${Math.round(e)}?feature=gene`),
    xrefs: id => request(`/xrefs/id/${enc(id)}?all_levels=1`),
    geneXrefs: (id, db) => request(`/xrefs/id/${enc(id)}?external_db=${enc(db)}`),
    sequenceId: (id, type) => request(`/sequence/id/${enc(id)}?type=${type}`),
    /** regions: array of "chr:start..end:strand"; returns Map(region -> sequence) */
    async sequenceRegions(species, regions) {
      const out = new Map();
      const chunks = [];
      // small chunks: Ensembl takes ~20 s per 50 regions under load, and long requests are the ones that time out
      for (let i = 0; i < regions.length; i += 25) chunks.push(regions.slice(i, i + 25));
      await Promise.all(chunks.map(async chunk => {
        const res = await request(`/sequence/region/${species}`, { method: 'POST', body: { regions: chunk } });
        res.forEach((r, i) => out.set(chunk.includes(r.query) ? r.query : chunk[i], r.seq));
      }));
      return out;
    },
  };
})();
