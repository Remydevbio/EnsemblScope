/* SVG genome viewer: ruler, neighbouring genes track, one row per transcript. */
class GenomeViewer {
  constructor(el, { onSelect, onGeneClick } = {}) {
    this.el = el;
    this.onSelect = onSelect;
    this.onGeneClick = onGeneClick;
    this.txs = []; this.genes = []; this.sel = null; this.view = [0, 1]; this.full = [0, 1];
    el.innerHTML = `
      <div class="gv-toolbar">
        <button class="btn" data-act="left" title="Pan left">◀</button>
        <button class="btn" data-act="right" title="Pan right">▶</button>
        <button class="btn" data-act="in" title="Zoom in">＋ Zoom in</button>
        <button class="btn" data-act="out" title="Zoom out">－ Zoom out</button>
        <button class="btn" data-act="reset" title="Fit gene">Fit gene</button>
        <button class="btn" data-act="seltx" title="Fit selected transcript">Fit selected</button>
        <span class="gv-loc"></span>
      </div>
      <div class="gv-canvas"><svg xmlns="http://www.w3.org/2000/svg"></svg><div class="gv-tip"></div></div>
      <div class="gv-legend">
        <span><i style="background:#b3262b"></i>coding exon (CDS)</span>
        <span><i style="background:#e9a3a3;height:6px"></i>UTR</span>
        <span><i style="background:#3d7fc1"></i>non-coding exon</span>
        <span><i style="background:#f2b631"></i>MANE Select</span>
      </div>`;
    this.svg = el.querySelector('svg');
    this.tip = el.querySelector('.gv-tip');
    this.locEl = el.querySelector('.gv-loc');
    this.canvas = el.querySelector('.gv-canvas');

    el.querySelector('.gv-toolbar').addEventListener('click', e => {
      const act = e.target.closest('[data-act]')?.dataset.act; if (!act) return;
      const [a, b] = this.view, w = b - a;
      if (act === 'left') this.setView(a - w * 0.3, b - w * 0.3);
      if (act === 'right') this.setView(a + w * 0.3, b + w * 0.3);
      if (act === 'in') this.zoom(0.5);
      if (act === 'out') this.zoom(2);
      if (act === 'reset') this.setView(...this.full);
      if (act === 'seltx') {
        const t = this.txs.find(t => t.id === this.sel);
        if (t) { const pad = (t.end - t.start) * 0.04 + 20; this.setView(t.start - pad, t.end + pad); }
      }
    });

    // drag to pan
    let drag = null;
    this.svg.addEventListener('pointerdown', e => {
      drag = { x: e.clientX, view: [...this.view], moved: false };
      this.svg.setPointerCapture(e.pointerId);
    });
    this.svg.addEventListener('pointermove', e => {
      if (drag) {
        const dx = e.clientX - drag.x;
        if (Math.abs(dx) > 3) { drag.moved = true; this.svg.classList.add('dragging'); this.hideTip(); }
        if (drag.moved) {
          const bp = dx / this.scale;
          this.view = [drag.view[0] - bp, drag.view[1] - bp];
          this.render();
        }
        return;
      }
      this.showTip(e);
    });
    const endDrag = e => {
      if (!drag) return;
      this.dragged = drag.moved; drag = null;
      this.svg.classList.remove('dragging');
      if (this.svg.hasPointerCapture?.(e.pointerId)) this.svg.releasePointerCapture(e.pointerId);
    };
    this.svg.addEventListener('pointerup', endDrag);
    this.svg.addEventListener('pointercancel', endDrag);
    this.svg.addEventListener('pointerleave', () => this.hideTip());

    this.svg.addEventListener('click', e => {
      if (this.dragged) { this.dragged = false; return; }
      const g = e.target.closest('[data-gene]');
      if (g && this.onGeneClick) return this.onGeneClick(g.dataset.gene);
      const t = e.target.closest('[data-tx]');
      if (t && this.onSelect) this.onSelect(t.dataset.tx);
    });

    this.svg.addEventListener('wheel', e => {
      e.preventDefault();
      const rect = this.svg.getBoundingClientRect();
      const px = e.clientX - rect.left;
      const center = this.view[0] + (px - this.LW) / this.scale;
      this.zoom(e.deltaY > 0 ? 1.25 : 0.8, center);
    }, { passive: false });

    new ResizeObserver(() => this.render()).observe(this.canvas);
  }

  setData({ chr, start, end, txs }) {
    this.chr = chr; this.txs = txs; this.genes = [];
    const pad = Math.max(500, (end - start) * 0.05);
    this.full = [start - pad, end + pad];
    this.view = [...this.full];
    this.render();
  }
  setGenes(genes) { this.genes = genes || []; this.render(); }
  setSelected(id) { this.sel = id; this.render(); }
  setView(a, b) { this.view = [a, b]; this.render(); }

  zoom(f, center) {
    const [a, b] = this.view;
    const c = center ?? (a + b) / 2;
    const span = (b - a) * f;
    if (span < 40 && f < 1) return;
    const maxSpan = (this.full[1] - this.full[0]) * 6;
    if (span > maxSpan && f > 1) return;
    this.setView(c - (c - a) * f, c + (b - c) * f);
  }

  render() {
    if (this._raf) return;
    this._raf = requestAnimationFrame(() => { this._raf = null; this._render(); });
  }

  _render() {
    const W = this.canvas.clientWidth || 900;
    const LW = this.LW = Math.round(Math.min(180, Math.max(110, W * 0.2)));
    const TW = W - LW - 8;
    const [a, b] = this.view;
    const sc = this.scale = TW / (b - a);
    const X = p => LW + (p - a) * sc;
    const clip = (s, e) => { // genomic [s,e] inclusive -> px [x1,x2] or null
      let x1 = X(s), x2 = X(e + 1);
      if (x2 < LW || x1 > W) return null;
      x1 = Math.max(x1, LW); x2 = Math.min(x2, W);
      if (x2 - x1 < 1) { const m = (x1 + x2) / 2; x1 = m - 0.5; x2 = m + 0.5; }
      return [x1, x2];
    };
    const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
    const fmt = n => Math.round(n).toLocaleString('en-US');
    const out = [];

    this.locEl.textContent = `${this.chr}:${fmt(Math.max(1, a))}-${fmt(b)} (${fmt(b - a)} bp)`;

    // ---- ruler
    const step = niceStep((b - a) / Math.max(3, Math.floor(TW / 110)));
    out.push(`<rect x="${LW}" y="0" width="${TW}" height="22" fill="#f4f6f9"/>`);
    out.push(`<text x="6" y="15" font-size="11" fill="#6b7785">chr ${esc(this.chr)}</text>`);
    for (let t = Math.ceil(a / step) * step; t <= b; t += step) {
      const x = X(t);
      out.push(`<line x1="${x}" x2="${x}" y1="14" y2="22" stroke="#8a95a3"/>`);
      out.push(`<text x="${x + 3}" y="12" font-size="10" fill="#4b5663">${fmtPos(t, step)}</text>`);
    }
    let y = 30;

    // ---- neighbouring genes track (packed)
    const vis = this.genes.filter(g => g.end >= a && g.start <= b);
    if (vis.length) {
      out.push(`<text x="6" y="${y + 10}" font-size="11" font-weight="600" fill="#22364d">Genes</text>`);
      const rows = [];
      const placed = vis.sort((p, q) => p.start - q.start).map(g => {
        const c = clip(g.start, g.end);
        if (!c) return null;
        const label = (g.external_name || g.id) + (g.strand > 0 ? ' >' : ' <');
        const right = Math.max(c[1], c[0] + label.length * 6.2) + 8;
        let r = rows.findIndex(end => end < c[0]);
        if (r < 0) { r = rows.length; rows.push(0); }
        rows[r] = right;
        return { g, c, r, label };
      }).filter(Boolean);
      placed.forEach(({ g, c, r, label }) => {
        const gy = y + r * 24;
        const col = g.biotype === 'protein_coding' ? '#b3262b' : '#3d7fc1';
        const cur = this.txs[0] && g.id === this.txs[0].geneId;
        out.push(`<g data-gene="${esc(g.id)}" data-tip="${esc(`${g.external_name || ''}  ${g.id}\n${g.biotype}\n${this.chr}:${fmt(g.start)}-${fmt(g.end)} (${g.strand > 0 ? '+' : '−'})${cur ? '' : '\nclick to open'}`)}" style="cursor:pointer">
          <rect x="${c[0]}" y="${gy}" width="${c[1] - c[0]}" height="7" fill="${col}" ${cur ? 'stroke="#f2b631" stroke-width="2"' : ''}/>
          <text x="${c[0]}" y="${gy + 18}" font-size="10.5" fill="${cur ? '#22364d' : '#555'}" font-weight="${cur ? 700 : 400}">${esc(label)}</text></g>`);
      });
      y += rows.length * 24 + 8;
      out.push(`<line x1="0" x2="${W}" y1="${y}" y2="${y}" stroke="#e3e7ec"/>`);
      y += 6;
    }

    // ---- transcripts
    const RH = 24;
    this.txs.forEach(tx => {
      const cy = y + RH / 2;
      const sel = tx.id === this.sel;
      const coding = !!tx.tl;
      const cdsCol = '#b3262b', utrCol = coding ? '#e9a3a3' : '#3d7fc1';
      out.push(`<g data-tx="${esc(tx.id)}" style="cursor:pointer">`);
      out.push(`<rect x="0" y="${y}" width="${W}" height="${RH}" fill="${sel ? '#fff5c2' : 'transparent'}"/>`);
      if (tx.mane) out.push(`<rect x="2" y="${y + 4}" width="4" height="${RH - 8}" fill="#f2b631"/>`);
      const name = tx.name.length > 22 ? tx.name.slice(0, 21) + '…' : tx.name;
      out.push(`<text x="10" y="${cy + 4}" font-size="11.5" font-weight="${sel ? 700 : 500}" fill="${coding ? '#7a1418' : '#1f4f86'}">${esc(name)}</text>`);
      const lc = clip(tx.start, tx.end);
      if (lc) {
        out.push(`<line x1="${lc[0]}" x2="${lc[1]}" y1="${cy}" y2="${cy}" stroke="#7d8794" stroke-width="1"/>`);
        // strand chevrons
        const d = tx.strand > 0 ? 3 : -3;
        let path = '';
        for (let x = lc[0] + 10; x < lc[1] - 6; x += 28) path += `M${x - d},${cy - 3}L${x},${cy}L${x - d},${cy + 3}`;
        if (path) out.push(`<path d="${path}" stroke="#7d8794" fill="none" stroke-width="1"/>`);
      }
      tx.exons.forEach((ex, i) => {
        const c = clip(ex.start, ex.end);
        if (!c) return;
        let tip = `${tx.name} — exon ${i + 1}/${tx.exons.length}\n${ex.id}\n${this.chr}:${fmt(ex.start)}-${fmt(ex.end)}\n${fmt(ex.end - ex.start + 1)} bp`;
        if (tx.tl) {
          const cs = Math.max(ex.start, tx.tl.start), ce = Math.min(ex.end, tx.tl.end);
          tip += cs <= ce ? `  (CDS ${fmt(ce - cs + 1)} bp)` : '  (UTR only)';
          out.push(`<rect x="${c[0]}" y="${cy - 4}" width="${c[1] - c[0]}" height="8" fill="${utrCol}" data-tip="${esc(tip)}"/>`);
          if (cs <= ce) {
            const cc = clip(cs, ce);
            if (cc) out.push(`<rect x="${cc[0]}" y="${cy - 7}" width="${cc[1] - cc[0]}" height="14" fill="${cdsCol}" data-tip="${esc(tip)}"/>`);
          }
        } else {
          out.push(`<rect x="${c[0]}" y="${cy - 6}" width="${c[1] - c[0]}" height="12" fill="${utrCol}" data-tip="${esc(tip)}"/>`);
        }
      });
      out.push(`</g>`);
      y += RH;
    });
    y += 4;
    // left label separator
    out.push(`<line x1="${LW - 2}" x2="${LW - 2}" y1="22" y2="${y}" stroke="#e3e7ec"/>`);

    this.svg.setAttribute('viewBox', `0 0 ${W} ${y}`);
    this.svg.setAttribute('height', y);
    this.svg.innerHTML = out.join('');
  }

  showTip(e) {
    const t = e.target.closest?.('[data-tip]');
    if (!t) return this.hideTip();
    const r = this.canvas.getBoundingClientRect();
    this.tip.textContent = t.dataset.tip;
    this.tip.style.display = 'block';
    let x = e.clientX - r.left + 14, y = e.clientY - r.top + 14;
    if (x + this.tip.offsetWidth > r.width) x = e.clientX - r.left - this.tip.offsetWidth - 10;
    this.tip.style.left = x + 'px';
    this.tip.style.top = y + 'px';
  }
  hideTip() { this.tip.style.display = 'none'; }
}

function niceStep(raw) {
  const p = Math.pow(10, Math.floor(Math.log10(Math.max(raw, 1))));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
function fmtPos(t, step) {
  if (step >= 1e6) return (t / 1e6).toFixed(0) + ' Mb';
  if (step >= 1e5) return (t / 1e6).toFixed(1) + ' Mb';
  if (step >= 1e4) return (t / 1e6).toFixed(2) + ' Mb';
  if (step >= 1e3) return (t / 1e3).toLocaleString('en-US', { maximumFractionDigits: 0 }) + ' kb';
  return Math.round(t).toLocaleString('en-US');
}
