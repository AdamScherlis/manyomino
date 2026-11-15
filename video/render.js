// Draws one frame of the music video per call to frame(f).  Driven by make.js.
const W = 1920, H = 1080, cv = document.getElementById('c'), ctx = cv.getContext('2d');
const small = new OffscreenCanvas(W / 4, H / 4), sctx = small.getContext('2d');
const RAMP = [[255, 236, 170], [255, 170, 96], [240, 96, 140], [160, 92, 236], [84, 116, 255]], PLATE = [13, 15, 22];
const NL = 64, LUT = [];
for (let i = 0; i < NL; i++) { const t = i / (NL - 1) * (RAMP.length - 1), k = Math.min(RAMP.length - 2, t | 0), g = t - k; LUT.push('rgb(' + [0, 1, 2].map(c => Math.round(RAMP[k][c] * (1 - g) + RAMP[k + 1][c] * g)).join(',') + ')'); }
const { NF, FPS, S16F } = SIM, sim = SIM.makeSim(Sampler), n = sim.n, S = sim.s;
const col = new Float32Array(n).fill(0.5);
// Camera path, computed ahead of time from a first pass of the same (deterministic) chain: for
// each frame the box that holds the shape over the next second or so, then smoothed in time.
const VH = H - 150;                    // the shape is framed above the caption strip
const CAM = (() => {
  const p = SIM.makeSim(Sampler), bx = new Float32Array(NF * 4);
  for (let f = 0; f < NF; f++) {
    p.step(f); const X = p.s.X, Y = p.s.Y; let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
    for (let i = 0; i < p.n; i++) { const x = X[i], y = Y[i]; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
    bx.set([x0, x1 + 1, y0, y1 + 1], f * 4);
  }
  const raw = new Float32Array(NF * 3);
  for (let f = 0; f < NF; f++) {
    let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
    for (let g = Math.max(0, f - 6); g <= Math.min(NF - 1, f + 20); g++) { x0 = Math.min(x0, bx[g * 4]); x1 = Math.max(x1, bx[g * 4 + 1]); y0 = Math.min(y0, bx[g * 4 + 2]); y1 = Math.max(y1, bx[g * 4 + 3]); }
    raw.set([(x0 + x1) / 2, (y0 + y1) / 2, Math.max(58, 1.1 * Math.max((x1 - x0) / 2 * VH / W, (y1 - y0) / 2))], f * 3);
  }
  const out = new Float32Array(NF * 3), SG = 11;
  for (let f = 0; f < NF; f++) {
    let w = 0, a = 0, b = 0, c = 0;
    for (let d = -3 * SG; d <= 3 * SG; d++) { const g = Math.max(0, Math.min(NF - 1, f + d)), k = Math.exp(-d * d / (2 * SG * SG)); w += k; a += k * raw[g * 3]; b += k * raw[g * 3 + 1]; c += k * raw[g * 3 + 2]; }
    out.set([a / w, b / w, c / w], f * 3);
  }
  return out;
})();
const cam = { x: 0, y: 0, h: 58 };
let ghosts = [], streaks = [], shake = 0, shakeF = -100, impact = 0;
const vig = ctx.createRadialGradient(W / 2, H / 2, H * 0.35, W / 2, H / 2, H * 1.05);
vig.addColorStop(0, 'rgba(0,0,0,0)'); vig.addColorStop(1, 'rgba(0,0,0,0.62)');

// captions: [first bar, last bar, text, kind]
const CAPS = [
  [0.3, 3.7, 'MANYOMINO', 'title'],
  [4.3, 7.5, 'move one cell at a time', 'cap'],
  [8.6, 13.6, 'or cut off a whole branch and paste it somewhere else', 'cap'],
  [32.5, 38.4, 'every move is exactly as likely as its reverse', 'cap'],
  [57.6, 63.4, 'so every shape of 2,500 cells is equally likely', 'cap'],
];
const smooth = (a, b, x) => { const t = Math.max(0, Math.min(1, (x - a) / (b - a))); return t * t * (3 - 2 * t); };

function frame(f, draw) {
  const r = sim.step(f), X = S.X, Y = S.Y, name = r.name, b = SIM.barOf(f), drop = name === 'drop1' || name === 'drop2';
  if (r.burst) {
    let k = 0;
    for (const e of r.burst) {
      if (!e.cells || k >= 8) break;
      const g = new Int32Array(e.cells.length * 2);
      e.cells.forEach((j, i) => { g[2 * i] = sim.pX[j]; g[2 * i + 1] = sim.pY[j]; });
      ghosts.push({ f, g, m: e.m });
      if (k < (drop ? 4 : 2) && e.m >= 12) streaks.push({ f, a: e.from, b: e.to, m: e.m });
      k++;
    }
    if (r.burst.length && r.burst[0].m >= 350) { shake = Math.min(1, r.burst[0].m / 900); shakeF = f; }
  }
  ghosts = ghosts.filter(g => f - g.f < 26); streaks = streaks.filter(s => f - s.f < 12);
  cam.x = CAM[f * 3]; cam.y = CAM[f * 3 + 1]; cam.h = CAM[f * 3 + 2];
  // ---- colour: graph distance from the middle of the shape, eased per cell
  const dist = S.depths(); let md = 1; for (let i = 0; i < n; i++) if (dist[i] > md) md = dist[i];
  const ease = f === 0 ? 1 : 0.1;
  for (let i = 0; i < n; i++) col[i] += (dist[i] / md - col[i]) * ease;
  if (!draw) return null;

  const beat = drop ? Math.exp(-((b * 4) % 1) * 4 * S16F / 4) : 0;       // decays over each quarter note
  for (const sct of SIM.SECTIONS) if (sct[0].startsWith('drop') && b >= sct[1] && b < sct[1] + 0.5) impact = Math.exp(-(b - sct[1]) * 16 * S16F / 5);
  const sAge = f - shakeF, sh = sAge < 14 ? shake * Math.exp(-sAge / 4) : 0;
  const s = (VH / 2) / cam.h * (1 + 0.014 * beat + 0.02 * sh);
  const ox = W / 2 - cam.x * s + sh * 9 * Math.sin(f * 2.4), oy = VH / 2 + 14 - cam.y * s + sh * 9 * Math.cos(f * 3.1);
  const gap = Math.max(0, Math.min(0.1 * s, (s - 4) * 0.12)), sz = s - gap + (gap > 0 ? 0 : 0.6);

  ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.globalCompositeOperation = 'source-over'; ctx.globalAlpha = 1; ctx.filter = 'none';
  ctx.fillStyle = `rgb(${PLATE})`; ctx.fillRect(0, 0, W, H);
  // ghosts: where a branch used to be
  for (const g of ghosts) {
    const a = 1 - (f - g.f) / 26; ctx.fillStyle = `rgba(120,150,255,${(0.34 * a * a).toFixed(3)})`;
    for (let i = 0; i < g.g.length; i += 2) ctx.fillRect(g.g[i] * s + ox, g.g[i + 1] * s + oy, sz, sz);
  }
  // cells, bucketed by colour
  const buckets = Array.from({ length: NL }, () => []);
  for (let i = 0; i < n; i++) buckets[Math.max(0, Math.min(NL - 1, (col[i] * (NL - 1)) | 0))].push(i);
  for (let k = 0; k < NL; k++) { const bk = buckets[k]; if (!bk.length) continue; ctx.fillStyle = LUT[k]; for (const i of bk) ctx.fillRect(X[i] * s + ox, Y[i] * s + oy, sz, sz); }
  // flashes: cells that just rode a transplanted branch, and cells that just relocated
  const sparkDur = r.R < 150 ? 16 : (r.R < 800 ? 9 : 5);
  const fb = Array.from({ length: 8 }, () => []);
  for (let i = 0; i < n; i++) {
    let a = 0; const ta = f - sim.tFrame[i], ra = f - sim.rFrame[i];
    if (ta < 40) a = Math.exp(-ta / 9) * (0.4 + 0.6 * Math.min(1, sim.tSize[i] / 7));
    if (ra < sparkDur) a = Math.max(a, 0.85 * (1 - ra / sparkDur));
    if (a > 0.04) fb[Math.min(7, (a * 8) | 0)].push(i);
  }
  for (let k = 0; k < 8; k++) { if (!fb[k].length) continue; ctx.fillStyle = `rgba(255,255,255,${((k + 0.5) / 8).toFixed(3)})`; for (const i of fb[k]) ctx.fillRect(X[i] * s + ox, Y[i] * s + oy, sz, sz); }
  // streaks: old place -> new place
  ctx.globalCompositeOperation = 'lighter'; ctx.lineCap = 'round';
  for (const st of streaks) {
    const a = 1 - (f - st.f) / 12; ctx.strokeStyle = `rgba(190,215,255,${(0.5 * a * a).toFixed(3)})`; ctx.lineWidth = Math.max(1.2, Math.min(3.5, 0.4 + Math.log2(st.m) * 0.3));
    const ax = (st.a[0] + 0.5) * s + ox, ay = (st.a[1] + 0.5) * s + oy, bx = (st.b[0] + 0.5) * s + ox, by = (st.b[1] + 0.5) * s + oy;
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.quadraticCurveTo((ax + bx) / 2 - (by - ay) * 0.22, (ay + by) / 2 + (bx - ax) * 0.22, bx, by); ctx.stroke();
  }
  // bloom
  sctx.globalCompositeOperation = 'copy'; sctx.drawImage(cv, 0, 0, W / 4, H / 4);
  ctx.filter = 'blur(10px)'; ctx.globalAlpha = 0.42 + 0.14 * beat + 0.2 * sh; ctx.drawImage(small, 0, 0, W, H);
  ctx.filter = 'blur(44px)'; ctx.globalAlpha = 0.3 + 0.1 * beat; ctx.drawImage(small, 0, 0, W, H);
  ctx.filter = 'none'; ctx.globalAlpha = 1; ctx.globalCompositeOperation = 'source-over';
  if (impact > 0.02) { ctx.globalCompositeOperation = 'lighter'; ctx.fillStyle = `rgba(200,215,255,${(0.5 * impact).toFixed(3)})`; ctx.fillRect(0, 0, W, H); ctx.globalCompositeOperation = 'source-over'; }
  impact = 0;
  ctx.fillStyle = vig; ctx.fillRect(0, 0, W, H);
  // text
  ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
  for (const [b0, b1, text, kind] of CAPS) {
    if (b < b0 || b > b1) continue;
    const a = smooth(b0, b0 + 0.35, b) * (1 - smooth(b1 - 0.35, b1, b));
    if (kind === 'title') {
      ctx.font = '300 62px "Inter Display", Inter, sans-serif'; ctx.letterSpacing = '36px'; ctx.fillStyle = `rgba(236,238,245,${(0.92 * a).toFixed(3)})`;
      ctx.fillText(text, W / 2 + 18, H - 120);
      ctx.font = '400 22px "JetBrains Mono", monospace'; ctx.letterSpacing = '5px'; ctx.fillStyle = `rgba(160,168,190,${(0.8 * a).toFixed(3)})`;
      ctx.fillText('a random walk through the polyominoes', W / 2, H - 72);
    } else {
      ctx.font = '300 38px "Inter Display", Inter, sans-serif'; ctx.letterSpacing = '1px';
      ctx.fillStyle = `rgba(0,0,0,${(0.5 * a).toFixed(3)})`; ctx.fillText(text, W / 2 + 1, H - 86 + 2);
      ctx.fillStyle = `rgba(236,238,245,${(0.94 * a).toFixed(3)})`; ctx.fillText(text, W / 2, H - 86);
    }
  }
  ctx.letterSpacing = '1px'; ctx.font = '400 19px "JetBrains Mono", monospace'; ctx.fillStyle = 'rgba(160,168,190,0.5)';
  ctx.textAlign = 'left'; ctx.fillText('moves ' + r.moves.toLocaleString('en-US'), 44, H - 38);
  ctx.textAlign = 'right'; ctx.fillText('n = 2,500', W - 44, H - 38);
  // fades
  const fade = Math.max(1 - smooth(0, 24, f), smooth(NF - 50, NF - 6, f));
  if (fade > 0) { ctx.fillStyle = `rgba(0,0,0,${fade.toFixed(3)})`; ctx.fillRect(0, 0, W, H); }
  return cv.toDataURL('image/jpeg', 0.95);
}
window.frame = frame;
