// In-browser port of the two kernels (plain versions; see src/manyomino.cpp for the fast one).
// Polyomino lives on a W x W torus that is much larger than its extent.
function Sampler(cells, lgW, seed) {
  const W = 1 << lgW, mask = W - 1, n = cells.length;
  const grid = new Int32Array(W * W);           // 0 = empty, else index + 1
  const X = new Int32Array(n), Y = new Int32Array(n);
  const mark = new Int32Array(n), lab = new Int8Array(n);
  const DX = [1, 0, -1, 0], DY = [0, 1, 0, -1];
  const G = [[1,0,0,1],[0,-1,1,0],[-1,0,0,-1],[0,1,-1,0],[1,0,0,-1],[0,1,1,0],[-1,0,0,1],[0,-1,-1,0]];
  let stamp = 0;
  // xorshift32-ish PRNG (mulberry32)
  let s = seed >>> 0;
  function rnd() { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }
  const below = k => (rnd() * k) | 0;
  const at = (x, y) => ((y & mask) << lgW) | (x & mask);
  for (let i = 0; i < n; i++) { X[i] = cells[i][0]; Y[i] = cells[i][1]; grid[at(X[i], Y[i])] = i + 1; }
  const q = [[], [], [], []], qh = [0, 0, 0, 0], par = [0, 0, 0, 0];
  const find = a => { while (par[a] !== a) a = par[a]; return a; };
  const stats = { rel: 0, relAcc: 0, cp: 0, cpAcc: 0, moved: 0 };
  const flash = new Float64Array(n);            // time of each cell's last ride on a transplanted branch
  let clock = 0;

  // is P \ c connected?  lockstep BFS from the neighbours of c
  function removable(ic) {
    const cx = X[ic], cy = Y[ic];
    stamp++; mark[ic] = stamp; lab[ic] = -1;
    let ng = 0;
    for (let d = 0; d < 4; d++) {
      const o = grid[at(cx + DX[d], cy + DY[d])]; if (!o) continue;
      q[ng].length = 0; q[ng].push(o - 1); qh[ng] = 0; par[ng] = ng; mark[o - 1] = stamp; lab[o - 1] = ng; ng++;
    }
    if (ng <= 1) return true;
    let ncls = ng;
    for (;;) {
      for (let g = 0; g < ng; g++) {
        if (qh[g] >= q[g].length) continue;
        const j = q[g][qh[g]++], x = X[j], y = Y[j];
        for (let d = 0; d < 4; d++) {
          const o = grid[at(x + DX[d], y + DY[d])]; if (!o) continue;
          const k = o - 1;
          if (mark[k] !== stamp) { mark[k] = stamp; lab[k] = g; q[g].push(k); }
          else if (lab[k] >= 0) { const a = find(lab[k]), b = find(g); if (a !== b) { par[a] = b; if (--ncls === 1) return true; } }
        }
        if (qh[g] >= q[g].length) {
          const r = find(g); let ex = true;
          for (let h = 0; h < ng; h++) if (find(h) === r && qh[h] < q[h].length) { ex = false; break; }
          if (ex) return false;
        }
      }
    }
  }

  function relocate() {
    stats.rel++;
    const ic = below(n), iu = below(n), d = below(4);
    if (iu === ic) return false;
    const sx = X[iu] + DX[d], sy = Y[iu] + DY[d];
    if (grid[at(sx, sy)]) return false;
    let k = 0;
    for (let e = 0; e < 4; e++) { const a = grid[at(sx + DX[e], sy + DY[e])]; if (a && a - 1 !== ic) k++; }
    if (k > 1 && below(k) !== 0) return false;
    if (!removable(ic)) return false;
    grid[at(X[ic], Y[ic])] = 0; X[ic] = sx; Y[ic] = sy; grid[at(sx, sy)] = ic + 1;
    stats.relAcc++;
    return true;
  }

  const B = [];
  function transplant() {
    stats.cp++;
    const ic = below(n), d = below(4), ia = below(n), g = below(8);
    const cx = X[ic], cy = Y[ic];
    const o = grid[at(cx + DX[d], cy + DY[d])];
    if (!o) return false;
    // B = component of P \ c containing c+d, abandoned once it exceeds (n-1)/2
    const cap = (n - 1) >> 1;
    stamp++; mark[ic] = stamp; lab[ic] = -1;
    B.length = 0; B.push(o - 1); mark[o - 1] = stamp; lab[o - 1] = 1;
    for (let h = 0; h < B.length; h++) {
      const j = B[h], x = X[j], y = Y[j];
      for (let e = 0; e < 4; e++) {
        const a = grid[at(x + DX[e], y + DY[e])]; if (!a) continue;
        const k = a - 1;
        if (mark[k] !== stamp) { mark[k] = stamp; lab[k] = 1; B.push(k); if (B.length > cap) return false; }
      }
    }
    if (mark[ia] === stamp && lab[ia] === 1) return false;       // c' inside B
    if (ia === ic && g === 0) return false;
    const ax = X[ia], ay = Y[ia], M = G[g], m = B.length;
    const nx = new Int32Array(m), ny = new Int32Array(m);
    for (let i = 0; i < m; i++) {
      const j = B[i], dx = X[j] - cx, dy = Y[j] - cy;
      const x = ax + M[0] * dx + M[1] * dy, y = ay + M[2] * dx + M[3] * dy;
      nx[i] = x; ny[i] = y;
      const a = grid[at(x, y)];
      if (a && !(mark[a - 1] === stamp && lab[a - 1] === 1)) return false;      // overlaps A
      for (let e = 0; e < 4; e++) {
        const b = grid[at(x + DX[e], y + DY[e])]; if (!b) continue;
        const k = b - 1;
        if (k === ia || (mark[k] === stamp && lab[k] === 1)) continue;
        return false;                                                          // touches A away from c'
      }
    }
    for (let i = 0; i < m; i++) grid[at(X[B[i]], Y[B[i]])] = 0;
    for (let i = 0; i < m; i++) { const j = B[i]; X[j] = nx[i]; Y[j] = ny[i]; grid[at(nx[i], ny[i])] = j + 1; }
    stats.cpAcc++; stats.moved += m;
    if (m >= 6) for (let i = 0; i < m; i++) flash[B[i]] = clock;
    return true;
  }

  function rg2() {
    let sx = 0, sy = 0, sxx = 0;
    const x0 = X[0], y0 = Y[0];
    for (let i = 0; i < n; i++) { const x = X[i] - x0, y = Y[i] - y0; sx += x; sy += y; sxx += x * x + y * y; }
    return sxx / n - (sx * sx + sy * sy) / (n * n);
  }
  // graph distance from a middle cell (double sweep), for colouring
  function depths() {
    const dist = new Int32Array(n), prev = new Int32Array(n), queue = new Int32Array(n);
    function sweep(src) {
      dist.fill(-1); dist[src] = 0; prev[src] = -1; let h = 0, t = 0; queue[t++] = src;
      while (h < t) { const v = queue[h++]; for (let d = 0; d < 4; d++) { const o = grid[at(X[v] + DX[d], Y[v] + DY[d])]; if (o && dist[o - 1] < 0) { dist[o - 1] = dist[v] + 1; prev[o - 1] = v; queue[t++] = o - 1; } } }
      return queue[t - 1];
    }
    const a = sweep(0), b = sweep(a);
    let v = b; for (let i = (dist[b] >> 1); i > 0; i--) v = prev[v];
    sweep(v);
    return dist;
  }
  return { n, X, Y, relocate, transplant, rg2, depths, stats, rnd, flash, branch: B, setClock: t => { clock = t; } };
}
if (typeof module !== 'undefined') module.exports = { Sampler };
