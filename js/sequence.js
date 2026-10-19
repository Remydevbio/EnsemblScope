/* Transcript sequence view: 5'UTR / CDS / 3'UTR parts per exon, intron heads & tails. */
const SeqView = (() => {
  const fmt = n => Number(n).toLocaleString('en-US');
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const PART_INFO = {
    u5: { label: "5′ UTR", txt: "5'UTR" },
    cd: { label: 'CDS', txt: 'CDS' },
    u3: { label: "3′ UTR", txt: "3'UTR" },
    nc: { label: 'non-coding exon', txt: 'exon' },
  };

  /** Fetch exon sequences and intron flanks; returns ordered segments. */
  async function build(species, tx, flank) {
    const st = tx.strand, chr = tx.chr;
    const reg = (s, e) => `${chr}:${s}..${e}:${st}`;
    const regions = [];
    const segs = [];
    let cdna = 1;

    tx.exons.forEach((ex, i) => {
      const len = ex.end - ex.start + 1;
      segs.push({ type: 'exon', n: i + 1, id: ex.id, start: ex.start, end: ex.end, len, cdnaStart: cdna, region: reg(ex.start, ex.end) });
      regions.push(reg(ex.start, ex.end));
      cdna += len;
      if (i < tx.exons.length - 1) {
        const nx = tx.exons[i + 1];
        const is = st > 0 ? ex.end + 1 : nx.end + 1;
        const ie = st > 0 ? nx.start - 1 : ex.start - 1;
        const ilen = ie - is + 1;
        const seg = { type: 'intron', n: i + 1, start: is, end: ie, len: ilen };
        if (ilen > 0) {
          if (ilen <= 2 * flank) {
            seg.full = reg(is, ie); regions.push(seg.full);
          } else {
            const head = st > 0 ? [is, is + flank - 1] : [ie - flank + 1, ie];
            const tail = st > 0 ? [ie - flank + 1, ie] : [is, is + flank - 1];
            seg.headR = reg(...head); seg.tailR = reg(...tail);
            regions.push(seg.headR, seg.tailR);
          }
        }
        segs.push(seg);
      }
    });

    const seqs = await API.sequenceRegions(species, regions);
    let cdsCount = 0;
    for (const s of segs) {
      if (s.type === 'exon') {
        s.seq = (seqs.get(s.region) || '').toUpperCase();
        s.parts = exonParts(s, tx);
        for (const p of s.parts) {
          if (p.cls === 'cd') { p.cdsStart = cdsCount + 1; cdsCount += p.seq.length; }
        }
      } else if (s.full) {
        s.seq = (seqs.get(s.full) || '').toLowerCase();
      } else if (s.headR) {
        s.head = (seqs.get(s.headR) || '').toLowerCase();
        s.tail = (seqs.get(s.tailR) || '').toLowerCase();
      }
    }
    return segs;
  }

  /** Split an exon into runs of 5'UTR / CDS / 3'UTR (transcript orientation). */
  function exonParts(seg, tx) {
    const tl = tx.tl, st = tx.strand, seq = seg.seq;
    const clsAt = k => {
      if (!tl) return 'nc';
      const g = st > 0 ? seg.start + k : seg.end - k;
      if (g >= tl.start && g <= tl.end) return 'cd';
      const before = st > 0 ? g < tl.start : g > tl.end;
      return before ? 'u5' : 'u3';
    };
    const parts = [];
    for (let k = 0; k < seq.length; k++) {
      const c = clsAt(k);
      const last = parts[parts.length - 1];
      if (last && last.cls === c) last.to = k; else parts.push({ cls: c, from: k, to: k });
    }
    parts.forEach(p => { p.seq = seq.slice(p.from, p.to + 1); p.cdnaStart = seg.cdnaStart + p.from; });
    return parts;
  }

  /** Render a sequence as numbered lines of 60, blocks of 10. clsFn(k) -> css class. */
  function seqLines(seq, startNum, clsFn, lineLen = 60) {
    let html = '';
    for (let i = 0; i < seq.length; i += lineLen) {
      let line = '', cur = null, buf = '';
      const end = Math.min(seq.length, i + lineLen);
      for (let k = i; k < end; k++) {
        const c = clsFn(k);
        if (c !== cur) { if (buf) line += `<span class="${cur}">${buf}</span>`; buf = ''; cur = c; }
        buf += seq[k];
        if ((k - i + 1) % 10 === 0 && k < end - 1) { line += `<span class="${cur}">${buf}</span> `; buf = ''; }
      }
      if (buf) line += `<span class="${cur}">${buf}</span>`;
      html += `<div class="sl"><span class="sn">${startNum + i}</span><span>${line}</span></div>`;
    }
    return html;
  }

  function render(container, tx, gene, segs, ctx) {
    const exons = segs.filter(s => s.type === 'exon');
    const sum = { u5: 0, cd: 0, u3: 0, nc: 0 };
    exons.forEach(e => e.parts.forEach(p => sum[p.cls] += p.seq.length));
    const cdna = exons.reduce((a, e) => a + e.len, 0);
    const copyStore = [];
    const copyBtn = (text) => { copyStore.push(text); return `<button class="btn" data-copy="${copyStore.length - 1}">Copy</button>`; };

    let html = `
      <div class="controls">
        <label>Intron bases shown at each end:
          <select id="flankSel">${[50, 100, 200, 500].map(v => `<option ${v === ctx.flank ? 'selected' : ''}>${v}</option>`).join('')}</select>
        </label>
      </div>
      <div class="seq-summary">
        <span>Location <b>${esc(tx.chr)}:${fmt(tx.start)}-${fmt(tx.end)}</b> (${tx.strand > 0 ? 'forward' : 'reverse'} strand)</span>
        <span>Exons <b>${exons.length}</b></span>
        <span>Transcript <b>${fmt(cdna)}</b> bp</span>
        ${tx.tl ? `<span>5′UTR <b>${fmt(sum.u5)}</b> bp</span><span>CDS <b>${fmt(sum.cd)}</b> bp</span><span>3′UTR <b>${fmt(sum.u3)}</b> bp</span><span>Protein <b>${fmt(tx.tl.aa)}</b> aa</span>` : `<span>Biotype <b>${esc(tx.biotype)}</b></span>`}
      </div>
      <div class="legend">
        ${tx.tl ? `<span class="u5">5′UTR</span><span class="cd">CDS</span><span class="u3">3′UTR</span>` : `<span class="nc">EXON</span>`}
        <span class="in" style="background:#f3f4f6">intron</span>
        <span class="ss" style="background:#f3f4f6">gt…ag splice sites</span>
        <span style="color:#6b7785">Numbers = cDNA position (exons) / intron position (introns)</span>
      </div>
      <div class="chips">${segs.map((s, i) => {
        if (s.type === 'exon') {
          const dom = s.parts.some(p => p.cls === 'cd') ? 'cd' : s.parts[0]?.cls || 'nc';
          return `<span class="chip ex ${dom}" data-jump="seg${i}">E${s.n}</span>`;
        }
        return `<span class="chip in" data-jump="seg${i}">I${s.n}</span>`;
      }).join('')}</div>`;

    segs.forEach((s, i) => {
      if (s.type === 'exon') {
        const desc = s.parts.map(p => `${PART_INFO[p.cls].txt} ${fmt(p.seq.length)}`).join(' + ');
        html += `<div class="seg exon" id="seg${i}">
          <div class="seg-head"><span class="t">Exon ${s.n}</span><span>${fmt(s.len)} bp</span>
            <span class="c">${esc(s.id)}</span><span class="c">${esc(tx.chr)}:${fmt(s.start)}-${fmt(s.end)}</span>
            <span class="muted">${desc}</span>${copyBtn(s.seq)}</div>`;
        s.parts.forEach(p => {
          const info = PART_INFO[p.cls];
          let lab = info.label;
          if (p.cls === 'cd') lab += ` · c.${fmt(p.cdsStart)}–${fmt(p.cdsStart + p.seq.length - 1)}`;
          lab += ` · ${fmt(p.seq.length)} bp`;
          html += `<div class="part"><div class="part-label lab-${p.cls}">${lab}</div>
            <div class="seqblock">${seqLines(p.seq, p.cdnaStart, () => p.cls)}</div></div>`;
        });
        html += `</div>`;
      } else {
        const loc = s.len > 0 ? `${esc(tx.chr)}:${fmt(s.start)}-${fmt(s.end)}` : '';
        let body = '', copyText = '';
        const L = s.len;
        if (s.len <= 0) {
          body = `<div class="gap">no intronic bases (adjacent exons)</div>`;
        } else if (s.seq != null) {
          body = `<div class="seqblock">${seqLines(s.seq, 1, k => (k < 2 || k >= L - 2) ? 'ss' : 'in')}</div>`;
          copyText = s.seq;
        } else {
          const h = s.head, t = s.tail;
          body = `<div class="seqblock">${seqLines(h, 1, k => k < 2 ? 'ss' : 'in')}</div>
            <div class="gap">··· ${fmt(L - h.length - t.length)} bp not shown ···</div>
            <div class="seqblock">${seqLines(t, L - t.length + 1, k => k >= t.length - 2 ? 'ss' : 'in')}</div>`;
          copyText = h + '\n' + t;
        }
        html += `<div class="seg intron" id="seg${i}">
          <div class="seg-head"><span class="t">Intron ${s.n}-${s.n + 1}</span><span>${fmt(Math.max(0, L))} bp</span>
            <span class="c">${loc}</span>${copyText ? copyBtn(copyText) : ''}</div>
          <div class="part">${body}</div></div>`;
      }
    });

    container.innerHTML = html;

    container.querySelector('#flankSel').addEventListener('change', e => ctx.onFlank(+e.target.value));
    container.querySelectorAll('[data-jump]').forEach(el => el.addEventListener('click', () =>
      document.getElementById(el.dataset.jump)?.scrollIntoView({ behavior: 'smooth', block: 'start' })));
    container.querySelectorAll('[data-copy]').forEach(el => el.addEventListener('click', async () => {
      try { await navigator.clipboard.writeText(copyStore[+el.dataset.copy]); el.textContent = 'Copied ✓'; }
      catch (_) { el.textContent = 'Copy failed'; }
      setTimeout(() => el.textContent = 'Copy', 1400);
    }));
  }

  /** Plain-text export of the annotated view. */
  function toText(tx, segs) {
    const wrap = s => s.replace(/(.{60})/g, '$1\n').trim();
    let out = `# ${tx.name} ${tx.id}  ${tx.chr}:${tx.start}-${tx.end} strand ${tx.strand > 0 ? '+' : '-'}\n`;
    out += `# Exons UPPERCASE, introns lowercase (intron ends only when long)\n\n`;
    segs.forEach(s => {
      if (s.type === 'exon') {
        s.parts.forEach(p => out += `>Exon_${s.n}_${PART_INFO[p.cls].txt} ${p.seq.length}bp\n${wrap(p.seq)}\n`);
      } else if (s.len > 0) {
        if (s.seq != null) out += `>Intron_${s.n} ${s.len}bp (full)\n${wrap(s.seq)}\n`;
        else out += `>Intron_${s.n} ${s.len}bp (first ${s.head.length} / last ${s.tail.length})\n${wrap(s.head)}\n...\n${wrap(s.tail)}\n`;
      }
    });
    return out;
  }

  return { build, render, toText };
})();
