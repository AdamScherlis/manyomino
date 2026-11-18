// In-browser port of src/zoom.cpp, working on a window.
//
// A level is a G x G picture (one byte per site) of the part of the shape around the point being
// zoomed into.  Refining it: the middle of the picture is doubled (every cell becomes a 2x2
// block) and relaxed by the local walk that adds and removes single cells at the critical weight
// per cell, for T sweeps.  The doubled cells in an outer collar are held fixed; they stand in
// for the rest of the shape, which is never needed.  The middle G x G of the result is the next
// level.  Nothing here depends on how deep the zoom is, so it can go on indefinitely.
function Refiner(G, seed, opt) {
  opt = opt || {};
  const PAD = 32, COL = 32, W = G + 2 * (PAD + COL), HALF = W >> 1;
  const Z = 1 / 4.0625696, KMAX = opt.kmax || 48, T = opt.sweeps || 32;
  const g = new Uint8Array(W * W);                 // bit 0: occupied; bits 1-3: search label + 1
  let cells = new Int32Array(1 << 17), n = 0;
  const QCAP = (2 * KMAX + 1) * (2 * KMAX + 1) + 8;
  const q = [new Int32Array(QCAP), new Int32Array(QCAP), new Int32Array(QCAP), new Int32Array(QCAP)];
  const qn = [0, 0, 0, 0], qh = [0, 0, 0, 0], par = [0, 0, 0, 0];
  const DP = [1, W, -1, -W];
  const RP = [1, 1 + W, W, W - 1, -1, -1 - W, -W, 1 - W];
  let s = seed >>> 0;
  function rnd() { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }
  const find = a => { while (par[a] !== a) a = par[a]; return a; };
  const movable = p => { const x = p % W, y = (p / W) | 0; return x >= COL && x < W - COL && y >= COL && y < W - COL; };
  let sweepsDone = 0, movesDone = 0, added = 0, removed = 0;

  // are the occupied neighbours of site p joined to each other, without p, inside the box of
  // radius K around p?
  function cert(p, K) {
    let m = 0;
    for (let i = 0; i < 8; i++) if (g[p + RP[i]] & 1) m |= 1 << i;
    const deg = (m & 1) + ((m >> 2) & 1) + ((m >> 4) & 1) + ((m >> 6) & 1);
    if (deg <= 1) return true;
    let links = 0;
    for (let i = 0; i < 8; i += 2) if ((m >> i & 1) && (m >> ((i + 1) & 7) & 1) && (m >> ((i + 2) & 7) & 1)) links++;
    if (links >= deg - 1) return true;
    if (K <= 1) return false;
    const px = p % W, py = (p / W) | 0;
    const x0 = Math.max(1, px - K), x1 = Math.min(W - 2, px + K), y0 = Math.max(1, py - K), y1 = Math.min(W - 2, py + K);
    const was = g[p] & 1; g[p] &= 254;
    let ng = 0;
    for (let i = 0; i < 8; i += 2) if (m >> i & 1) { const j = p + RP[i]; q[ng][0] = j; qn[ng] = 1; qh[ng] = 0; g[j] |= (ng + 1) << 1; par[ng] = ng; ng++; }
    let ncls = ng, ok = false, dead = false;
    while (!ok && !dead) {
      let any = false;
      for (let a = 0; a < ng && !ok && !dead; a++) {
        if (qh[a] >= qn[a]) continue;
        const j = q[a][qh[a]++]; any = true;
        const x = j % W, y = (j / W) | 0;
        for (let d = 0; d < 4; d++) {
          const ax = x + (d === 0 ? 1 : d === 2 ? -1 : 0), ay = y + (d === 1 ? 1 : d === 3 ? -1 : 0);
          if (ax < x0 || ax > x1 || ay < y0 || ay > y1) continue;
          const b = j + DP[d], v = g[b];
          if (!(v & 1)) continue;
          const lb = (v >> 1) & 7;
          if (!lb) { g[b] |= (a + 1) << 1; q[a][qn[a]++] = b; }
          else { const ra = find(a), rb = find(lb - 1); if (ra !== rb) { par[ra] = rb; if (--ncls === 1) ok = true; } }
        }
        if (!ok) {                                  // a class whose searches have run dry can never join the others
          const ra = find(a); let alive = false;
          for (let c = 0; c < ng; c++) if (find(c) === ra && qh[c] < qn[c]) alive = true;
          if (!alive) dead = true;
        }
      }
      if (!any) break;
    }
    for (let a = 0; a < ng; a++) { const qa = q[a]; for (let i = 0, e = qn[a]; i < e; i++) g[qa[i]] &= 1; }
    if (was) g[p] |= 1;
    return ok;
  }

  function move() {
    const i = (rnd() * n) | 0, c = cells[i], r = (rnd() * 8) | 0;
    if (r < 4) {
      const t = c + DP[r];
      if ((g[t] & 1) || !movable(t)) return;
      let k = 0; for (let d = 0; d < 4; d++) { const u = t + DP[d]; if ((g[u] & 1) && movable(u)) k++; }
      if (k > 1 && ((rnd() * k) | 0) !== 0) return;
      const a = 4 * Z * n / (n + 1);
      if (a < 1 && rnd() >= a) return;
      const u = rnd(), K = u * KMAX < 1 ? KMAX : (1 / u) | 0;
      if (!cert(t, K)) return;
      g[t] |= 1; if (n === cells.length) { const c2 = new Int32Array(n * 2); c2.set(cells); cells = c2; }
      cells[n++] = t; added++;
    } else {
      if (n <= 1) return;
      let k = 0; for (let d = 0; d < 4; d++) { const u = c + DP[d]; if ((g[u] & 1) && movable(u)) k++; }
      if (!k) return;
      const a = n / (4 * Z * (n - 1));
      if (a < 1 && rnd() >= a) return;
      const u = rnd(), K = u * KMAX < 1 ? KMAX : (1 / u) | 0;
      if (!cert(c, K)) return;
      g[c] &= 254; cells[i] = cells[--n]; removed++;
    }
  }

  // begin refining `coarse` (G x G bytes) about the coarse site (fx, fy)
  function start(coarse, fx, fy) {
    g.fill(0); n = 0; sweepsDone = 0; movesDone = 0; added = removed = 0;
    const h = HALF >> 1, sx0 = fx - h, sy0 = fy - h;
    for (let cy = 0; cy < HALF; cy++) {
      const yy = sy0 + cy; if (yy < 0 || yy >= G) continue;
      for (let cx = 0; cx < HALF; cx++) {
        const xx = sx0 + cx; if (xx < 0 || xx >= G || !coarse[xx + yy * G]) continue;
        const p = 2 * cx + 2 * cy * W;
        g[p] = g[p + 1] = g[p + W] = g[p + W + 1] = 1;
      }
    }
    for (let y = COL; y < W - COL; y++) for (let x = COL; x < W - COL; x++) if (g[x + y * W]) {
      if (n === cells.length) { const c2 = new Int32Array(n * 2); c2.set(cells); cells = c2; }
      cells[n++] = x + y * W;
    }
  }
  // run for about `ms` milliseconds; true when the level is finished
  function step(ms) {
    const t0 = performance.now();
    while (sweepsDone < T) {
      if (n < 1) { sweepsDone = T; break; }
      const batch = Math.max(1, n >> 3);
      for (let k = 0; k < batch; k++) move();
      movesDone += batch; sweepsDone += 0.125;
      if (performance.now() - t0 >= ms) break;
    }
    return sweepsDone >= T;
  }
  function result() {
    const out = new Uint8Array(G * G), o = PAD + COL;
    for (let y = 0; y < G; y++) { const r = (y + o) * W + o; for (let x = 0; x < G; x++) out[x + y * G] = g[r + x] & 1; }
    return out;
  }
  return { start, step, result, G, W, PAD, COL, T, get n() { return n; }, get progress() { return sweepsDone / T; }, get moves() { return movesDone; },
           grid: g, limit: (G >> 1) - (HALF >> 1) };
}

// the occupied site nearest to (x, y) in a G x G picture, within radius rmax; null if none
function nearestCell(pic, G, x, y, rmax) {
  let best = null, bd = 1e18;
  for (let r = 0; r <= rmax; r++) {
    if (best && r * r > bd) break;
    for (let dy = -r; dy <= r; dy++) for (let dx = -r; dx <= r; dx++) {
      if (Math.max(Math.abs(dx), Math.abs(dy)) !== r) continue;
      const xx = x + dx, yy = y + dy; if (xx < 0 || yy < 0 || xx >= G || yy >= G || !pic[xx + yy * G]) continue;
      const d = dx * dx + dy * dy; if (d < bd) { bd = d; best = [xx, yy]; }
    }
  }
  return best;
}
if (typeof module !== 'undefined') module.exports = { Refiner, nearestCell };
