/* Main controller: search, gene summary, transcript table, viewer, sequence panel. */
const SPECIES = {
  human: { api: 'homo_sapiens', ens: 'Homo_sapiens', ucsc: 'hg38', name: 'Human' },
  mouse: { api: 'mus_musculus', ens: 'Mus_musculus', ucsc: 'mm39', name: 'Mouse' },
};

const $ = s => document.querySelector(s);
const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = n => Number(n).toLocaleString('en-US');
const ext = (url, text, title = '') => `<a href="${esc(url)}" target="_blank" rel="noopener" title="${esc(title)}">${esc(text)}</a>`;

const state = { sp: 'human', gene: null, txs: [], selected: null, sort: { key: null, dir: 1 }, flank: 100, token: 0, seqToken: 0, segs: null };

const viewer = new GenomeViewer($('#viewer'), {
  onSelect: id => selectTranscript(id),
  onGeneClick: id => { if (id !== state.gene?.id) location.hash = `#/${state.sp}/${id}`; },
});

/* ---------------- routing ---------------- */
$('#searchForm').addEventListener('submit', e => {
  e.preventDefault();
  const q = $('#query').value.trim();
  if (q) location.hash = `#/${$('#species').value}/${encodeURIComponent(q)}`;
});
window.addEventListener('hashchange', route);
route();

function route() {
  const m = location.hash.match(/^#\/(human|mouse)\/([^/]+)(?:\/([^/]+))?/);
  if (!m) return;
  const [, sp, qRaw, tx] = m;
  const q = decodeURIComponent(qRaw);
  $('#species').value = sp;
  $('#query').value = q;
  const g = state.gene;
  if (g && state.sp === sp && (g.id === q || (g.display_name || '').toLowerCase() === q.toLowerCase())) {
    if (tx && tx !== state.selected) selectTranscript(tx, false);
    return;
  }
  loadGene(sp, q, tx);
}

function setStatus(msg, isErr = false) {
  const el = $('#status');
  if (!msg) return el.classList.add('hidden');
  el.className = 'status' + (isErr ? ' error' : '');
  el.innerHTML = msg;
}

/* ---------------- gene loading ---------------- */
async function loadGene(sp, query, txId) {
  const token = ++state.token;
  state.sp = sp; state.selected = null; state.segs = null; state.entrez = null; state.sort = { key: null, dir: 1 };
  setStatus(`<span class="spinner"></span>Looking up <b>${esc(query)}</b> in ${SPECIES[sp].name}…`);
  $('#welcome').classList.add('hidden');
  ['#genePanel', '#txPanel', '#viewerPanel', '#seqPanel'].forEach(s => $(s).classList.add('hidden'));

  let gene;
  try { gene = await API.lookupGene(SPECIES[sp].api, query); }
  catch (e) {
    if (token !== state.token) return;
    setStatus(`Could not load <b>${esc(query)}</b>: ${esc(e.message)}`, true);
    return;
  }
  if (token !== state.token) return;
  setStatus('');
  state.gene = gene;
  state.txs = (gene.Transcript || []).map(t => normTx(t, gene));
  document.title = `${gene.display_name || gene.id} — EnsemblScope`;

  renderGene();
  renderTable();
  ['#genePanel', '#txPanel', '#viewerPanel'].forEach(s => $(s).classList.remove('hidden'));
  viewer.setData({ chr: gene.seq_region_name, start: gene.start, end: gene.end, txs: visibleTxs() });

  // async enrichments
  enrichTags(token);
  enrichXrefs(token);
  API.geneXrefs(gene.id, 'EntrezGene').then(x => {
    if (token !== state.token || !x?.length) return;
    state.entrez = x[0].primary_id; renderGene();
  }).catch(() => {});
  const span = gene.end - gene.start;
  if (span < 4.8e6) {
    const pad = Math.min(span * 2.5 + 5000, (4.9e6 - span) / 2);
    API.overlapGenes(SPECIES[sp].api, gene.seq_region_name, gene.start - pad, gene.end + pad)
      .then(gs => { if (token === state.token) viewer.setGenes(gs); }).catch(() => {});
  }

  const want = txId && state.txs.find(t => t.id === txId) ? txId
    : (state.txs.find(t => t.mane) || state.txs.find(t => t.canonical) || {}).id;
  if (want) selectTranscript(want, !!txId);
}

function normTx(t, gene) {
  const strand = t.strand;
  const exons = [...(t.Exon || [])].sort((a, b) => strand > 0 ? a.start - b.start : b.end - a.end);
  const length = exons.reduce((s, e) => s + e.end - e.start + 1, 0);
  const tr = t.Translation;
  const maneArr = Array.isArray(t.MANE) ? t.MANE : (t.MANE ? [t.MANE] : []);
  const mSel = maneArr.find(m => /select/i.test(m.type || '')) || null;
  const mPc = maneArr.find(m => /clinical/i.test(m.type || '')) || null;
  return {
    id: t.id, version: t.version, name: t.display_name || t.id, biotype: t.biotype,
    geneId: gene.id, chr: t.seq_region_name || gene.seq_region_name,
    start: t.start, end: t.end, strand, exons, length,
    tl: tr ? { id: tr.id, start: tr.start, end: tr.end, aa: tr.length } : null,
    canonical: !!t.is_canonical,
    mane: !!mSel, manePC: !!mPc, maneRefseq: (mSel || mPc)?.refseq_match || null,
    tags: [], tsl: null, ccds: null, refseq: null, uniprot: null,
  };
}

function enrichTags(token) {
  API.overlapTranscripts(state.gene.id).then(list => {
    if (token !== state.token) return;
    const byId = new Map((list || []).map(f => [f.transcript_id || f.id, f]));
    state.txs.forEach(t => {
      const f = byId.get(t.id); if (!f) return;
      const tags = Array.isArray(f.tag) ? f.tag : (f.tag ? String(f.tag).split(',') : []);
      t.tags = tags;
      if (tags.includes('MANE_Select')) t.mane = true;
      if (tags.includes('MANE_Plus_Clinical')) t.manePC = true;
      if (tags.includes('Ensembl_canonical')) t.canonical = true;
      if (f.transcript_support_level) t.tsl = String(f.transcript_support_level).split(' ')[0];
      if (f.ccdsid && !t.ccds) t.ccds = f.ccdsid;
    });
    renderTable();
    viewer.txs = visibleTxs(); viewer.render();
  }).catch(() => {});
}

function enrichXrefs(token) {
  let pendingRender = null;
  const scheduleRender = () => { if (!pendingRender) pendingRender = setTimeout(() => { pendingRender = null; if (token === state.token) renderTable(); }, 150); };
  state.txs.forEach(t => {
    API.xrefs(t.id).then(xs => {
      if (token !== state.token) return;
      const refseq = new Set(), uni = new Set(), trembl = new Set();
      (xs || []).forEach(x => {
        const db = x.dbname || '';
        if (db === 'CCDS' && !t.ccds) t.ccds = x.primary_id;
        else if (/^RefSeq_(mRNA|ncRNA)/.test(db)) refseq.add(x.display_id || x.primary_id);
        else if (db === 'Uniprot/SWISSPROT') uni.add(x.primary_id);
        else if (db === 'Uniprot/SPTREMBL') trembl.add(x.primary_id);
      });
      if (t.maneRefseq) refseq.add(t.maneRefseq);
      t.refseq = [...refseq].sort((a, b) => (b === t.maneRefseq) - (a === t.maneRefseq) || a.localeCompare(b));
      t.uniprot = [...uni]; t.trembl = [...trembl];
      scheduleRender();
    }).catch(() => { t.refseq = t.refseq || []; t.uniprot = t.uniprot || []; scheduleRender(); });
  });
}

/* ---------------- gene panel ---------------- */
function renderGene() {
  const g = state.gene, sp = SPECIES[state.sp];
  const desc = (g.description || '').replace(/\s*\[Source:.*?\]/, '');
  const chr = g.seq_region_name;
  const ucscChr = chr === 'MT' ? 'chrM' : 'chr' + chr;
  $('#geneTitle').innerHTML = `${esc(g.display_name || g.id)} <span class="muted">${esc(desc)}</span>`;
  const links = [
    ext(`https://www.ensembl.org/${sp.ens}/Gene/Summary?g=${g.id}`, 'Ensembl'),
    state.entrez ? ext(`https://www.ncbi.nlm.nih.gov/gene/${state.entrez}`, 'NCBI Gene') :
      ext(`https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent((g.display_name || g.id) + '[sym] AND ' + sp.name + '[orgn]')}`, 'NCBI Gene'),
    ext(`https://genome.ucsc.edu/cgi-bin/hgTracks?db=${sp.ucsc}&position=${ucscChr}:${g.start}-${g.end}`, 'UCSC'),
    state.sp === 'human' ? ext(`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(g.display_name || '')}`, 'GeneCards')
      : ext(`https://www.informatics.jax.org/searchtool/Search.do?query=${encodeURIComponent(g.display_name || g.id)}`, 'MGI'),
  ];
  $('#geneLinks').innerHTML = links.join('');
  $('#geneInfo').innerHTML = `
    <dt>Gene ID</dt><dd>${esc(g.id)}${g.version ? '.' + g.version : ''}</dd>
    <dt>Location</dt><dd>Chromosome ${esc(chr)}: ${fmt(g.start)}-${fmt(g.end)} ${g.strand > 0 ? 'forward' : 'reverse'} strand · ${fmt(g.end - g.start + 1)} bp · ${esc(g.assembly_name || '')}</dd>
    <dt>Biotype</dt><dd>${esc(g.biotype)}</dd>
    <dt>Transcripts</dt><dd>${state.txs.length} (${state.txs.filter(t => t.tl).length} protein-coding)</dd>
    ${g.source ? `<dt>Annotation</dt><dd>${esc(g.source)}${g.logic_name ? ' · ' + esc(g.logic_name) : ''}</dd>` : ''}`;
}

/* ---------------- transcript table ---------------- */
function defaultRank(t) { return t.mane ? 0 : t.manePC ? 1 : t.canonical ? 2 : t.tl ? 3 : 4; }
function sortedTxs() {
  const { key, dir } = state.sort;
  const arr = [...state.txs];
  if (!key) return arr.sort((a, b) => defaultRank(a) - defaultRank(b) || (b.tl?.aa || 0) - (a.tl?.aa || 0) || b.length - a.length || a.name.localeCompare(b.name));
  const val = t => key === 'aa' ? (t.tl?.aa || 0) : key === 'ccds' ? (t.ccds || '~') : t[key];
  return arr.sort((a, b) => {
    const va = val(a), vb = val(b);
    return (typeof va === 'number' ? va - vb : String(va).localeCompare(String(vb), undefined, { numeric: true })) * dir;
  });
}

$('#txTable thead').addEventListener('click', e => {
  const th = e.target.closest('th[data-sort]'); if (!th) return;
  const k = th.dataset.sort;
  state.sort = state.sort.key === k ? { key: k, dir: -state.sort.dir } : { key: k, dir: (k === 'length' || k === 'aa') ? -1 : 1 };
  renderTable();
  viewer.txs = visibleTxs(); viewer.render();
});
$('#onlyCoding').addEventListener('change', () => { renderTable(); viewer.txs = visibleTxs(); viewer.render(); });
function visibleTxs() { return sortedTxs().filter(t => !$('#onlyCoding').checked || t.biotype === 'protein_coding'); }

function renderTable() {
  const sp = SPECIES[state.sp];
  document.querySelectorAll('#txTable th').forEach(th => {
    th.classList.toggle('sorted', th.dataset.sort === state.sort.key);
    th.classList.toggle('asc', th.dataset.sort === state.sort.key && state.sort.dir > 0);
  });
  const rows = visibleTxs();
  $('#txCount').textContent = `(${rows.length}${rows.length !== state.txs.length ? ' of ' + state.txs.length : ''})`;
  const pend = '<span class="pending">…</span>';
  $('#txTable tbody').innerHTML = rows.map(t => {
    const refseq = t.refseq == null ? pend : t.refseq.length ? t.refseq.slice(0, 3).map(r =>
      ext(`https://www.ncbi.nlm.nih.gov/nuccore/${r}`, r, 'Open RefSeq record at NCBI (CDS annotated)')).join(', ')
      + (t.refseq.length > 3 ? ` <span class="muted">+${t.refseq.length - 3}</span>` : '') : '-';
    const ccds = t.ccds ? ext(`https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=${t.ccds}`, t.ccds, 'Consensus CDS at NCBI')
      : (t.refseq == null && t.tl ? pend : '-');
    const uni = t.uniprot == null ? pend : (t.uniprot.length ? t.uniprot : (t.trembl || []).slice(0, 1)).map(u =>
      ext(`https://www.uniprot.org/uniprotkb/${u}/entry`, u)).join(', ') || '-';
    const flags = [
      t.mane ? '<span class="badge mane" title="MANE Select">MANE Select</span>' : '',
      t.manePC ? '<span class="badge manepc">MANE Plus Clinical</span>' : '',
      t.canonical ? '<span class="badge canon" title="Ensembl canonical">Canonical</span>' : '',
      t.tags.includes('basic') ? '<span class="badge basic">basic</span>' : '',
      t.tsl ? `<span class="badge tsl" title="Transcript support level">TSL:${esc(t.tsl)}</span>` : '',
      t.tags.includes('appris_principal_1') ? '<span class="badge">APPRIS P1</span>' : '',
    ].join('');
    return `<tr data-tx="${esc(t.id)}" class="${t.id === state.selected ? 'selected' : ''}">
      <td class="tx-name">${esc(t.name)}</td>
      <td>${ext(`https://www.ensembl.org/${sp.ens}/Transcript/Summary?t=${t.id}`, t.id + (t.version ? '.' + t.version : ''), 'Open in Ensembl')}</td>
      <td class="num">${fmt(t.length)}</td>
      <td class="num">${t.tl ? fmt(t.tl.aa) + ' aa' : '<span class="muted">no protein</span>'}</td>
      <td>${esc(t.biotype.replace(/_/g, ' '))}</td>
      <td>${ccds}</td><td>${refseq}</td><td>${uni}</td><td>${flags}</td></tr>`;
  }).join('');
}

$('#txTable tbody').addEventListener('click', e => {
  if (e.target.closest('a')) return;
  const tr = e.target.closest('tr[data-tx]');
  if (tr) selectTranscript(tr.dataset.tx);
});

/* ---------------- transcript detail ---------------- */
async function selectTranscript(id, scroll = true) {
  const tx = state.txs.find(t => t.id === id);
  if (!tx) return;
  state.selected = id;
  document.querySelectorAll('#txTable tbody tr').forEach(r => r.classList.toggle('selected', r.dataset.tx === id));
  viewer.setSelected(id);
  history.replaceState(null, '', `#/${state.sp}/${state.gene.display_name || state.gene.id}/${id}`);

  const sp = SPECIES[state.sp];
  $('#seqPanel').classList.remove('hidden');
  $('#seqTitle').innerHTML = `${esc(tx.name)} <span class="muted">${esc(tx.id)} · ${esc(tx.biotype.replace(/_/g, ' '))}${tx.mane ? ' · MANE Select' : ''}</span>`;
  const dl = [
    `<button class="btn" data-dl="cdna">cDNA FASTA</button>`,
    tx.tl ? `<button class="btn" data-dl="cds">CDS FASTA</button>` : '',
    tx.tl ? `<button class="btn" data-dl="protein">Protein FASTA</button>` : '',
    `<button class="btn" data-dl="annot">Annotated .txt</button>`,
    ext(`https://www.ensembl.org/${sp.ens}/Transcript/Exons?t=${tx.id}`, 'Ensembl exons'),
  ];
  $('#seqTools').innerHTML = dl.join('');
  $('#seqTools').onclick = e => { const b = e.target.closest('[data-dl]'); if (b) download(tx, b.dataset.dl, b); };

  if (scroll) $('#seqPanel').scrollIntoView({ behavior: 'smooth', block: 'start' });
  await loadSequence(tx);
}

async function loadSequence(tx) {
  const token = ++state.seqToken;
  const body = $('#seqBody');
  body.innerHTML = `<p><span class="spinner"></span>Fetching exon and intron sequences from Ensembl…</p>`;
  try {
    const segs = await SeqView.build(SPECIES[state.sp].api, tx, state.flank);
    if (token !== state.seqToken) return;
    state.segs = segs;
    SeqView.render(body, tx, state.gene, segs, {
      flank: state.flank,
      onFlank: v => { state.flank = v; loadSequence(tx); },
    });
  } catch (e) {
    if (token !== state.seqToken) return;
    body.innerHTML = `<div class="status error">Could not load sequence: ${esc(e.message)}</div>`;
  }
}

async function download(tx, type, btn) {
  const save = (text, name) => {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
    a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  };
  if (type === 'annot') {
    if (state.segs) save(SeqView.toText(tx, state.segs), `${tx.name}_annotated.txt`);
    return;
  }
  const old = btn.textContent; btn.textContent = '…';
  try {
    const r = await API.sequenceId(tx.id, type);
    const seq = (r.seq || '').replace(/(.{60})/g, '$1\n');
    save(`>${tx.id} ${tx.name} ${type}\n${seq}\n`, `${tx.name}_${type}.fa`);
  } catch (e) { alert('Download failed: ' + e.message); }
  btn.textContent = old;
}
