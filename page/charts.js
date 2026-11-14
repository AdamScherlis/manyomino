// Minimal SVG charts: log/linear axes, points with error bars, lines, bands, hover, table view.
function fmtNum(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e15 || a < 1e-4) { const e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e); return (Math.abs(m - 1) < 1e-9 ? '' : m.toFixed(1) + '×') + '10' + String(e).replace(/-/, '⁻').replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]); }
  if (a >= 1e6) { const e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e); if (Math.abs(m - Math.round(m)) < 1e-9 && Math.round(m) === 1) return '10' + String(e).replace(/\d/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]); return (v / 1e6).toLocaleString('en-US', { maximumFractionDigits: 1 }) + 'M'; }
  if (a >= 100) return Math.round(v).toLocaleString('en-US');
  return String(+v.toPrecision(4));
}
function sig(v, n) { return (+v.toPrecision(n)).toLocaleString('en-US', { maximumSignificantDigits: n }); }

function chart(el, cfg) {
  const H = cfg.height || 340;
  const tip = document.createElement('div'); tip.className = 'tip'; tip.hidden = true;
  const holder = document.createElement('div'); holder.className = 'plot';
  el.innerHTML = ''; el.appendChild(holder); holder.appendChild(tip);
  if (cfg.series.filter(s => s.name).length >= 2) {
    const lg = document.createElement('div'); lg.className = 'legend';
    cfg.series.filter(s => s.name).forEach(s => {
      lg.insertAdjacentHTML('beforeend', `<span><svg width="22" height="10" aria-hidden="true">${s.pts ? `<circle cx="11" cy="5" r="4" fill="${s.color}"/>` : `<line x1="1" y1="5" x2="21" y2="5" stroke="${s.color}" stroke-width="2" ${s.dash ? 'stroke-dasharray="4 3"' : ''}/>`}</svg>${s.name}</span>`);
    });
    el.insertBefore(lg, holder);
  }
  // table view
  const det = document.createElement('details'); det.className = 'tableview';
  let th = `<summary>Show the numbers</summary><div class="scroll"><table><thead><tr><th>series</th><th>${cfg.xname || cfg.xlabel}</th><th>${cfg.yname || cfg.ylabel}</th><th>±</th></tr></thead><tbody>`;
  cfg.series.forEach(s => (s.pts || []).forEach(p => { th += `<tr><td>${s.name || ''}</td><td>${fmtNum(p[0])}</td><td>${sig(p[1], cfg.ysig || 5)}</td><td>${p[2] ? sig(p[2], 2) : ''}</td></tr>`; }));
  det.innerHTML = th + '</tbody></table></div>';
  if (!cfg.noTable) el.appendChild(det);

  const T = { x: v => cfg.xlog ? Math.log10(v) : v, y: v => cfg.ylog ? Math.log10(v) : v };
  function ticks(lo, hi, log, maxN) {
    if (log) {
      const out = []; const a = Math.ceil(lo - 1e-9), b = Math.floor(hi + 1e-9);
      const step = Math.max(1, Math.ceil((b - a + 1) / maxN));
      for (let e = a; e <= b; e += step) out.push(e);
      if (out.length < 3) { out.length = 0; for (let e = Math.floor(lo); e <= b; e++) for (const m of [1, 2, 5]) { const v = Math.log10(m) + e; if (v >= lo - 1e-9 && v <= hi + 1e-9) out.push(v); } }
      return out;
    }
    const span = hi - lo, raw = span / maxN, p = Math.pow(10, Math.floor(Math.log10(raw)));
    const st = [1, 2, 5, 10].map(m => m * p).find(s => s >= raw);
    const out = []; for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + 1e-9; v += st) out.push(+v.toFixed(12));
    return out;
  }
  let hot = [];
  function render() {
    const W = Math.max(280, holder.clientWidth), narrow = W < 380;
    const m = { l: 58, r: cfg.rightPad != null ? (narrow ? Math.min(cfg.rightPad, 20) : cfg.rightPad) : 18, t: 14, b: 44 };
    let xs = [], ys = [];
    cfg.series.forEach(s => { (s.pts || []).forEach(p => { xs.push(p[0]); ys.push(p[1] - (p[2] || 0), p[1] + (p[2] || 0)); }); (s.line || []).forEach(p => { xs.push(p[0]); ys.push(p[1]); }); });
    (cfg.bands || []).forEach(b => ys.push(b.y0, b.y1));
    let [x0, x1] = cfg.xdomain || [Math.min(...xs), Math.max(...xs)], [y0, y1] = cfg.ydomain || [Math.min(...ys), Math.max(...ys)];
    x0 = T.x(x0); x1 = T.x(x1); y0 = T.y(y0); y1 = T.y(y1);
    if (!cfg.xdomain) { const p = (x1 - x0) * 0.04; x0 -= p; x1 += p; }
    if (!cfg.ydomain) { const p = (y1 - y0) * 0.07; y0 -= p; y1 += p; }
    const X = v => m.l + (T.x(v) - x0) / (x1 - x0) * (W - m.l - m.r), Y = v => H - m.b - (T.y(v) - y0) / (y1 - y0) * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="${cfg.aria || ''}">`;
    const xt = ticks(x0, x1, cfg.xlog, narrow ? 4 : 7), yt = ticks(y0, y1, cfg.ylog, 6);
    const lab = (t, log, f) => {
      if (f) return f(log ? Math.pow(10, t) : t);
      if (log && Math.abs(t - Math.round(t)) < 1e-9 && (t >= 6 || t <= -4)) return `10<tspan dy="-5" font-size="8.5">${Math.round(t)}</tspan>`;
      return fmtNum(log ? Math.pow(10, t) : t);
    };
    yt.forEach(t => { const y = H - m.b - (t - y0) / (y1 - y0) * (H - m.t - m.b); s += `<line class="grid" x1="${m.l}" x2="${W - m.r}" y1="${y}" y2="${y}"/><text class="tick" x="${m.l - 8}" y="${y + 4}" text-anchor="end">${lab(t, cfg.ylog, cfg.yfmt)}</text>`; });
    xt.forEach(t => { const x = m.l + (t - x0) / (x1 - x0) * (W - m.l - m.r); s += `<line class="grid" x1="${x}" x2="${x}" y1="${m.t}" y2="${H - m.b}"/><text class="tick" x="${x}" y="${H - m.b + 17}" text-anchor="middle">${lab(t, cfg.xlog, cfg.xfmt)}</text>`; });
    s += `<line class="axis" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/>`;
    s += `<text class="alabel" x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle">${cfg.xlabel}</text>`;
    (cfg.bands || []).forEach(b => {
      s += `<rect x="${m.l}" width="${W - m.l - m.r}" y="${Y(b.y1)}" height="${Y(b.y0) - Y(b.y1)}" fill="${b.color}" opacity=".13"/>`;
      if (b.mid != null) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(b.mid)}" y2="${Y(b.mid)}" stroke="${b.color}" stroke-width="1.5"/>`;
      if (b.label) s += b.left ? `<text class="dlabel" x="${m.l + 8}" y="${Y(b.labelAt != null ? b.labelAt : b.y1) - 6}">${b.label}</text>` : `<text class="dlabel" x="${W - m.r - 6}" y="${Y(b.labelAt != null ? b.labelAt : b.y1) - 6}" text-anchor="end">${b.label}</text>`;
    });
    hot = [];
    cfg.series.forEach(se => {
      if (se.line) s += `<polyline fill="none" stroke="${se.color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round" ${se.dash ? 'stroke-dasharray="5 4"' : ''} points="${se.line.map(p => X(p[0]).toFixed(1) + ',' + Y(p[1]).toFixed(1)).join(' ')}"/>`;
      (se.pts || []).forEach(p => {
        const x = X(p[0]), y = Y(p[1]);
        if (p[2] && Math.abs(Y(p[1] + p[2]) - Y(Math.max(p[1] - p[2], 1e-300))) > 9) s += `<line x1="${x}" x2="${x}" y1="${Y(p[1] + p[2])}" y2="${Y(p[1] - p[2])}" stroke="${se.color}" stroke-width="1.5"/>`;
        s += `<circle cx="${x}" cy="${y}" r="${se.r || 4.5}" fill="${se.color}" stroke="var(--panel)" stroke-width="2"/>`;
        hot.push({ x, y, p, se });
      });
      if (se.label && !narrow) { const a = se.labelAt || (se.line || se.pts)[(se.line || se.pts).length - 1]; s += `<text class="dlabel" x="${X(a[0]) + (se.labelDx != null ? se.labelDx : 8)}" y="${Y(a[1]) + (se.labelDy != null ? se.labelDy : 4)}" text-anchor="${se.labelAnchor || 'start'}">${se.label}</text>`; }
    });
    s += `<circle class="ring" r="8" fill="none" stroke="var(--ink)" stroke-width="1.5" visibility="hidden"/></svg>`;
    holder.querySelectorAll('svg').forEach(e => e.remove());
    holder.insertAdjacentHTML('afterbegin', s);
  }
  holder.addEventListener('pointermove', e => {
    const r = holder.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
    let best = null, bd = 34 * 34;
    hot.forEach(h => { const d = (h.x - px) ** 2 + (h.y - py) ** 2; if (d < bd) { bd = d; best = h; } });
    const ring = holder.querySelector('.ring');
    if (!best) { tip.hidden = true; ring.setAttribute('visibility', 'hidden'); return; }
    ring.setAttribute('cx', best.x); ring.setAttribute('cy', best.y); ring.setAttribute('visibility', 'visible');
    tip.innerHTML = (cfg.tip ? cfg.tip(best.p, best.se) : `<b>${best.se.name || ''}</b><br>${cfg.xname || cfg.xlabel}: ${fmtNum(best.p[0])}<br>${cfg.yname || cfg.ylabel}: ${sig(best.p[1], cfg.ysig || 4)}${best.p[2] ? ' ± ' + sig(best.p[2], 2) : ''}`);
    tip.hidden = false;
    const tw = tip.offsetWidth, left = Math.min(Math.max(4, best.x + 12), r.width - tw - 4);
    tip.style.left = left + 'px'; tip.style.top = Math.max(4, best.y - tip.offsetHeight - 10) + 'px';
  });
  holder.addEventListener('pointerleave', () => { tip.hidden = true; const r = holder.querySelector('.ring'); if (r) r.setAttribute('visibility', 'hidden'); });
  new ResizeObserver(render).observe(holder);
  render();
}
