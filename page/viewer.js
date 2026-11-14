// Pan/zoom viewer for a polyomino shipped as one nibble per cell (see analysis/poly.py).
const RAMP = [[255, 236, 170], [255, 170, 96], [240, 96, 140], [160, 92, 236], [84, 116, 255]];
const PLATE = [13, 15, 22];
function rampColor(t) {
  t = Math.max(0, Math.min(1, t)) * (RAMP.length - 1);
  const i = Math.min(RAMP.length - 2, t | 0), f = t - i;
  return [0, 1, 2].map(c => Math.round(RAMP[i][c] * (1 - f) + RAMP[i + 1][c] * f));
}
const LUT = new Uint32Array(256);
for (let i = 0; i < 256; i++) { const c = rampColor(i / 255); LUT[i] = (255 << 24) | (c[2] << 16) | (c[1] << 8) | c[0]; }
const PLATE32 = (255 << 24) | (PLATE[2] << 16) | (PLATE[1] << 8) | PLATE[0];

function decodePoly(buf) {
  const dv = new DataView(buf), n = dv.getInt32(0, true);
  const b = new Uint8Array(buf, 12);
  const X = new Int32Array(n), Y = new Int32Array(n), D = new Int32Array(n);
  const sx = new Int32Array(n + 4), sy = new Int32Array(n + 4), sd = new Int32Array(n + 4);
  const DX = [1, 0, -1, 0], DY = [0, 1, 0, -1];
  let sp = 0, i = 0, maxD = 0;
  sx[0] = dv.getInt32(4, true); sy[0] = dv.getInt32(8, true); sd[0] = 0; sp = 1;
  while (sp > 0) {
    sp--; const x = sx[sp], y = sy[sp], d = sd[sp];
    X[i] = x; Y[i] = y; D[i] = d; if (d > maxD) maxD = d;
    const m = (b[i >> 1] >> ((i & 1) << 2)) & 15; i++;
    for (let k = 3; k >= 0; k--) if (m >> k & 1) { sx[sp] = x + DX[k]; sy[sp] = y + DY[k]; sd[sp] = d + 1; sp++; }
  }
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
  for (let j = 0; j < n; j++) { if (X[j] < x0) x0 = X[j]; if (X[j] > x1) x1 = X[j]; if (Y[j] < y0) y0 = Y[j]; if (Y[j] > y1) y1 = Y[j]; }
  const C = new Uint8Array(n);
  for (let j = 0; j < n; j++) C[j] = (D[j] * 255 / maxD) | 0;
  return { n, X, Y, D, C, maxD, x0, x1, y0, y1 };
}

function PolyViewer(root, poly, opts) {
  opts = opts || {};
  const canvas = root.querySelector('canvas'), ctx = canvas.getContext('2d');
  const readout = root.querySelector('[data-readout]');
  let w = 0, h = 0, dpr = 1, img = null, buf32 = null;
  let cx = (poly.x0 + poly.x1 + 1) / 2, cy = (poly.y0 + poly.y1 + 1) / 2, scale = 1, fitScale = 1;
  let pending = false, active = false;

  function resize() {
    const r = canvas.getBoundingClientRect();
    dpr = Math.min(2, window.devicePixelRatio || 1);
    const nw = Math.max(1, Math.round(r.width * dpr)), nh = Math.max(1, Math.round(r.height * dpr));
    const first = !img;
    if (nw === w && nh === h) return;
    const rel = first ? 1 : scale / fitScale;
    w = nw; h = nh; canvas.width = w; canvas.height = h;
    img = ctx.createImageData(w, h); buf32 = new Uint32Array(img.data.buffer);
    fitScale = Math.min(w / (poly.x1 - poly.x0 + 1), h / (poly.y1 - poly.y0 + 1)) * 0.94;
    scale = fitScale * rel;
    draw();
  }
  function fit() { cx = (poly.x0 + poly.x1 + 1) / 2; cy = (poly.y0 + poly.y1 + 1) / 2; scale = fitScale; request(); }
  function request() { if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; draw(); }); } }

  function draw() {
    if (!img) return;
    buf32.fill(PLATE32);
    const { n, X, Y, C } = poly, s = scale, ox = w / 2 - cx * s, oy = h / 2 + cy * s;
    if (s <= 1.25) {
      for (let i = 0; i < n; i++) {
        const px = (X[i] * s + ox) | 0, py = (oy - (Y[i] + 1) * s) | 0;
        if (px >= 0 && px < w && py >= 0 && py < h) buf32[py * w + px] = LUT[C[i]];
      }
    } else {
      const gap = s >= 7 ? 1 : 0;
      for (let i = 0; i < n; i++) {
        const fx = X[i] * s + ox, fy = oy - (Y[i] + 1) * s;
        if (fx <= -s || fx >= w || fy <= -s || fy >= h) continue;
        const xa = Math.max(0, Math.round(fx)), xb = Math.min(w, Math.round(fx + s) - gap);
        const ya = Math.max(0, Math.round(fy)), yb = Math.min(h, Math.round(fy + s) - gap), col = LUT[C[i]];
        for (let yy = ya; yy < yb; yy++) { const row = yy * w; for (let xx = xa; xx < xb; xx++) buf32[row + xx] = col; }
      }
    }
    ctx.putImageData(img, 0, 0);
    // scale bar: a round number of cells, 60-160 css px long
    const target = 110 * dpr / s; const p10 = Math.pow(10, Math.floor(Math.log10(target)));
    const len = [1, 2, 5, 10].map(m => m * p10).reduce((a, b) => Math.abs(b - target) < Math.abs(a - target) ? b : a);
    const L = len * s, bx = 14 * dpr, by = h - 14 * dpr;
    ctx.strokeStyle = 'rgba(233,235,242,.85)'; ctx.lineWidth = 1.5 * dpr; ctx.beginPath();
    ctx.moveTo(bx, by - 4 * dpr); ctx.lineTo(bx, by); ctx.lineTo(bx + L, by); ctx.lineTo(bx + L, by - 4 * dpr); ctx.stroke();
    ctx.fillStyle = 'rgba(233,235,242,.85)'; ctx.font = (11 * dpr) + 'px ui-monospace,Menlo,Consolas,monospace';
    ctx.fillText(len.toLocaleString('en-US') + (len === 1 ? ' cell' : ' cells'), bx, by - 8 * dpr);
    if (readout) {
      const z = s / dpr;
      readout.textContent = z >= 1 ? 'one cell = ' + (z >= 10 ? Math.round(z) : z.toFixed(1)) + ' px' : 'one px = ' + (1 / z >= 10 ? Math.round(1 / z) : (1 / z).toFixed(1)) + ' cells';
    }
  }
  function zoomAt(px, py, f) {     // px,py in canvas device pixels
    const ns = Math.max(fitScale * 0.6, Math.min(40 * dpr, scale * f));
    const wx = cx + (px - w / 2) / scale, wy = cy - (py - h / 2) / scale;
    scale = ns; cx = wx - (px - w / 2) / scale; cy = wy + (py - h / 2) / scale;
    request();
  }
  const pos = e => { const r = canvas.getBoundingClientRect(); return [(e.clientX - r.left) * dpr, (e.clientY - r.top) * dpr]; };
  function setActive(a) { active = a; root.classList.toggle('active', a); }

  // pointers: drag to pan, two fingers to pinch.  One finger scrolls the page unless active.
  const pts = new Map(); let pinch = 0;
  canvas.addEventListener('pointerdown', e => {
    if (e.pointerType === 'touch' && !active && pts.size === 0) { pts.set(e.pointerId, pos(e)); return; }
    setActive(true); pts.set(e.pointerId, pos(e)); canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener('pointermove', e => {
    if (!pts.has(e.pointerId)) return;
    const p = pos(e), q = pts.get(e.pointerId);
    if (pts.size === 1 && active) { cx -= (p[0] - q[0]) / scale; cy += (p[1] - q[1]) / scale; request(); }
    pts.set(e.pointerId, p);
    if (pts.size === 2) {
      const [a, b] = [...pts.values()], d = Math.hypot(a[0] - b[0], a[1] - b[1]);
      if (pinch) zoomAt((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, d / pinch);
      pinch = d; setActive(true);
    }
  });
  const up = e => { pts.delete(e.pointerId); pinch = 0; };
  canvas.addEventListener('pointerup', up); canvas.addEventListener('pointercancel', up);
  canvas.addEventListener('wheel', e => {
    if (!active && !e.ctrlKey) return;            // don't hijack page scroll until clicked
    e.preventDefault(); const p = pos(e); zoomAt(p[0], p[1], Math.exp(-e.deltaY * (e.ctrlKey ? 0.01 : 0.0015)));
  }, { passive: false });
  canvas.addEventListener('dblclick', e => { const p = pos(e); zoomAt(p[0], p[1], 2.5); });
  document.addEventListener('pointerdown', e => { if (!root.contains(e.target)) setActive(false); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape') setActive(false); });
  root.querySelectorAll('[data-zoom]').forEach(b => b.addEventListener('click', () => {
    const k = b.dataset.zoom; setActive(true);
    if (k === 'fit') fit(); else zoomAt(w / 2, h / 2, k === 'in' ? 2 : 0.5);
  }));
  root.querySelectorAll('[data-goto]').forEach(b => b.addEventListener('click', () => {
    // fly to the cell farthest from the centre (a tip) at single-cell resolution
    setActive(true);
    let best = 0; for (let i = 1; i < poly.n; i++) if (poly.D[i] > poly.D[best]) best = i;
    const tx = poly.X[best] + .5, ty = poly.Y[best] + .5, ts = 9 * dpr, s0 = scale, x0 = cx, y0 = cy, t0 = performance.now();
    const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
    (function step(t) {
      const u = reduce ? 1 : Math.min(1, (t - t0) / 2600), e = u < .5 ? 4 * u * u * u : 1 - Math.pow(-2 * u + 2, 3) / 2;
      scale = s0 * Math.pow(ts / s0, e);
      // pan so the target stays on a straight screen-space course
      const k = (1 / s0 - 1 / scale) / (1 / s0 - 1 / ts || 1);
      cx = x0 + (tx - x0) * k; cy = y0 + (ty - y0) * k; draw();
      if (u < 1) requestAnimationFrame(step);
    })(t0);
  }));
  new ResizeObserver(resize).observe(canvas);
  resize();
  return { fit, draw };
}
