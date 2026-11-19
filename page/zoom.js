// In-browser port of src/zoom.cpp, working on a window.
//
// A level is a G x G picture (one byte per site, 255 where there is a cell) of the part of the
// shape around the point being zoomed into.  Refining it: the middle half of the picture is
// doubled (every cell becomes a 2x2 block), giving G x G sites again, and relaxed by the local
// walk that adds and removes single cells at the critical weight per cell.
//
// The walk is not run evenly over the window.  In the middle (out to R0 of the way to the edge)
// every site gets the full T sweeps.  From there the number of sweeps falls off, halving OCT
// times before R1, and beyond R1 nothing moves at all.  This is done in stages: the walk runs
// first on everything inside R1 for a moment, then on a slightly smaller square for twice as
// long, and so on inward, with whatever lies outside the current square held fixed.  Each stage
// is the same reversible walk on a smaller region, so the law it aims at is unchanged; the
// effect is only that the outer part is less far along.  At its edge the new level is therefore
// exactly the old one drawn with four cells for one, and it turns into fully refined shape
// gradually on the way in.  A picture assembled from nested levels then has no seams: where one
// level hands over to the next coarser they are the same picture.
//
// The next level must be centred within LIMIT cells of the middle, so that all of it lies in
// the fully relaxed part of this one.  Nothing here depends on how deep the zoom is.
// where the sweeps fall off: full out to R0 of the way from the middle to the edge, halving OCT
// times from there to R1, none beyond
const ZOOM = { R0: 0.72, R1: 0.94, OCT: 12 };
// sweeps received at relative distance r from the middle (0 middle, 1 edge), as a share of the full number
function zoomRate(r) { const m = Math.ceil((r - ZOOM.R0) / (ZOOM.R1 - ZOOM.R0) * ZOOM.OCT - 1e-9); return m <= 0 ? 1 : m > ZOOM.OCT ? 0 : Math.pow(2, -m); }
// how far from the middle of a level (of G sites) the next may be centred
function zoomLimit(G) { return Math.floor((ZOOM.R0 - 0.5) * (G >> 1)) - 2; }

function Refiner(G, seed, opt) {
  opt = opt || {};
  const W = G, C = G >> 1, Z = 1 / 4.0625696, KMAX = opt.kmax || 48, T = opt.sweeps || 32;
  const R0 = ZOOM.R0, R1 = ZOOM.R1, OCT = ZOOM.OCT;
  const KA = opt.ka || 0.5;                        // the certificate box has radius K >= r with probability KA / r (r >= 2)
  const g = new Uint8Array(W * W);                 // bit 0: occupied; bits 1-3: search label + 1; bit 4: held fixed and off the list;
                                                   // bit 5: a listed cell that may not be removed (the rim of a pin)
  // stage m works on the square rings out to STAGE[m][0], for STAGE[m][1] sweeps; m runs downward
  const STAGE = []; let TOTAL = 0;
  for (let m = 0; m <= OCT; m++) STAGE.push([Math.floor(C * (R0 + (R1 - R0) * m / OCT)) - 1, T * (m < OCT ? Math.pow(2, -m - 1) : Math.pow(2, -OCT))]);
  let cells = new Int32Array(1 << 18), n = 0;
  const QCAP = (2 * KMAX + 1) * (2 * KMAX + 1) + 8;
  const q = [new Int32Array(QCAP), new Int32Array(QCAP), new Int32Array(QCAP), new Int32Array(QCAP)];
  const qn = [0, 0, 0, 0], qh = [0, 0, 0, 0], par = [0, 0, 0, 0];
  const DP = [1, W, -1, -W];
  const RP = [1, 1 + W, W, W - 1, -1, -1 - W, -W, 1 - W];
  let s = seed >>> 0;
  function rnd() { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }
  const find = a => { while (par[a] !== a) a = par[a]; return a; };
  let stage = -1, stageSweeps = 0, movesDone = 0, workDone = 0, workAll = 1;

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
    for (let a = 0; a < ng; a++) { const qa = q[a]; for (let i = 0, e = qn[a]; i < e; i++) g[qa[i]] &= 49; }
    if (was) g[p] |= 1;
    return ok;
  }

  // One proposal.  Pick a listed cell and one of 8 choices: add at one of its four neighbouring
  // sites (thinned so that every empty site next to listed cells is proposed equally often), or
  // remove the cell.  Sites outside the current square carry bit 4 and are never touched.
  function move() {
    const i = (rnd() * n) | 0, c = cells[i], r = (rnd() * 8) | 0;
    if (r < 4) {
      const t = c + DP[r];
      if (g[t] & 17) return;
      let k = 0; for (let d = 0; d < 4; d++) if ((g[t + DP[d]] & 17) === 1) k++;
      if (k > 1 && ((rnd() * k) | 0) !== 0) return;
      const a = 4 * Z * n / (n + 1);
      if (a < 1 && rnd() >= a) return;
      const u = rnd(), K = u >= KA ? 1 : u * KMAX < KA ? KMAX : (KA / u) | 0;
      if (!cert(t, K)) return;
      g[t] |= 1; if (n === cells.length) { const c2 = new Int32Array(n * 2); c2.set(cells); cells = c2; }
      cells[n++] = t;
    } else {
      if (n <= 1 || (g[c] & 32)) return;
      let k = 0; for (let d = 0; d < 4; d++) if ((g[c + DP[d]] & 17) === 1) k++;
      if (!k) return;
      const a = n / (4 * Z * (n - 1));
      if (a < 1 && rnd() >= a) return;
      const u = rnd(), K = u >= KA ? 1 : u * KMAX < KA ? KMAX : (KA / u) | 0;
      if (!cert(c, K)) return;
      g[c] &= 254; cells[i] = cells[--n];
    }
  }

  // narrow the working square to stage m: everything outside it is marked fixed
  function enter(m) {
    stage = m; stageSweeps = 0; const rmax = STAGE[m][0], lo = C - 1 - rmax, hi = C + rmax;      // sites with lo <= x, y <= hi
    if (m === OCT) {
      n = 0;
      for (let y = 0; y < W; y++) for (let x = 0; x < W; x++) {
        const p = x + y * W, inside = x >= lo && x <= hi && y >= lo && y <= hi;
        if (!inside) { g[p] |= 16; continue; }
        if ((g[p] & 17) === 1) { if (n === cells.length) { const c2 = new Int32Array(n * 2); c2.set(cells); cells = c2; } cells[n++] = p; }
      }
    } else {
      const plo = C - 1 - STAGE[m + 1][0], phi = C + STAGE[m + 1][0];
      for (let y = plo; y <= phi; y++) for (let x = plo; x <= phi; x++) if (x < lo || x > hi || y < lo || y > hi) g[x + y * W] |= 16;
      let k = 0; for (let i = 0; i < n; i++) if (!(g[cells[i]] & 16)) cells[k++] = cells[i];
      n = k;
    }
  }
  // The shape is connected, so no part of it can sit inside a window without reaching the
  // window's edge.  Setting a pin into a doubled level can leave such islands: a piece of the pin
  // with nothing outside to go on into, or a piece of the surroundings whose only link ran
  // through the middle.  Each is given the shortest run of cells, through empty sites outside
  // the pin, to the nearest part that does reach the edge, if there is one within REACH sites,
  // and as many cells are taken away again from the tips nearest the link.  Coming down from
  // the far end this happens at the coarsest level where the island first appears, where the
  // gap is a few sites; every finer level then inherits the link and the walk roughens it like
  // any other strand.  (A link is new shape, and one made several levels out grows into a
  // whole strand by the time it is seen, so the levels of a batch are kept few: see Outward.)
  let joined = 0, joinCells = 0, missed = 0; const REACH = 48;
  function ground() {
    const st = [], flood = () => { while (st.length) { const p = st.pop(); for (let d = 0; d < 4; d++) { const t = p + DP[d]; if (t >= 0 && t < W * W && (d & 1 || ((t / W) | 0) === ((p / W) | 0)) && (g[t] & 65) === 1) { g[t] |= 64; st.push(t); } } } };
    for (let k = 0; k < W; k++) for (const p of [k, k + (W - 1) * W, k * W, k * W + W - 1]) if ((g[p] & 65) === 1) { g[p] |= 64; st.push(p); }
    flood();
    const prev = new Map();
    for (let p0 = 0; p0 < W * W; p0++) {
      if ((g[p0] & 193) !== 1) continue;
      // an island: collect it, then search breadth-first from its empty neighbours for the grounded part
      const isle = [p0]; g[p0] |= 128; st.push(p0);
      while (st.length) { const p = st.pop(); for (let d = 0; d < 4; d++) { const t = p + DP[d]; if ((g[t] & 193) === 1) { g[t] |= 128; st.push(t); isle.push(t); } } }
      prev.clear(); let fr = [], goal = -1;
      for (const p of isle) for (let d = 0; d < 4; d++) { const t = p + DP[d]; if (g[t] === 0 && !prev.has(t)) { prev.set(t, -1); fr.push(t); } }
      const reach = isle.length >= 256 ? 6 * REACH : REACH;                 // a large piece is worth a long link
      for (let dist = 0; dist < reach && fr.length && goal < 0; dist++) {
        const nx = [];
        for (const p of fr) {
          const x = p % W, y = (p / W) | 0; if (x < 1 || y < 1 || x > W - 2 || y > W - 2) continue;
          for (let d = 0; d < 4; d++) { const t = p + DP[d]; if (g[t] & 64) { goal = p; break; } if (g[t] === 0 && !prev.has(t)) { prev.set(t, p); nx.push(t); } }
          if (goal >= 0) break;
        }
        fr = nx;
      }
      if (goal < 0) { missed++; continue; }
      for (const p of isle) g[p] |= 64;
      joined++; const path = new Set(); for (let p = goal; p >= 0; p = prev.get(p)) { g[p] = 65; joinCells++; path.add(p); }
      // the link must not add to the shape: take as many cells away again from the tips nearest to it
      const gx = goal % W, gy = (goal / W) | 0; let need = path.size;
      for (let rad = 4; need > 0 && rad <= 128; rad *= 2) {
        for (let pass = 0; pass < 3 && need > 0; pass++)
          for (let y = Math.max(1, gy - rad); y <= Math.min(W - 2, gy + rad) && need > 0; y++) for (let x = Math.max(1, gx - rad); x <= Math.min(W - 2, gx + rad) && need > 0; x++) {
            const p = x + y * W; if ((g[p] & 49) !== 1 || path.has(p)) continue;
            // a cell with one neighbour, or the outer corner of a block (two neighbours joined to each other diagonally)
            const e = g[p + 1] & 1, w = g[p - 1] & 1, u = g[p + W] & 1, v = g[p - W] & 1, deg = e + w + u + v;
            if (deg === 1 || (deg === 2 && ((e && u && (g[p + 1 + W] & 1)) || (e && v && (g[p + 1 - W] & 1)) || (w && u && (g[p - 1 + W] & 1)) || (w && v && (g[p - 1 - W] & 1))))) { g[p] = 0; need--; }
          }
      }
    }
    for (let p = 0; p < W * W; p++) g[p] &= 63;
  }
  // begin refining `coarse` (G x G bytes) about the coarse site (fx, fy).  With a pin (S x S
  // bytes), the middle S x S sites of the new level are set to it and held: its cells stay, its
  // empty sites stay empty, and the walk goes on around it.  Cells on the pin's outer ring are
  // kept on the list where `rim` says so (4 S flags: top, bottom, left, right), so that the shape
  // outside can grow from them as from any other cell; the rest of the pin is inert.
  function start(coarse, fx, fy, pin, S, rim) {
    g.fill(0); movesDone = 0; workDone = 0;
    const h = C >> 1, sx0 = fx - h, sy0 = fy - h; let all = 0;
    for (let cy = 0; cy < C; cy++) {
      const yy = sy0 + cy; if (yy < 0 || yy >= G) continue;
      for (let cx = 0; cx < C; cx++) {
        const xx = sx0 + cx; if (xx < 0 || xx >= G || !coarse[xx + yy * G]) continue;
        const p = 2 * cx + 2 * cy * W; all += 4;
        g[p] = g[p + 1] = g[p + W] = g[p + W + 1] = 1;
      }
    }
    if (pin) {
      const o = C - (S >> 1);
      for (let y = 0; y < S; y++) for (let x = 0; x < S; x++)
        g[(x + o) + (y + o) * W] = !pin[x + y * S] ? 16 : ((y === 0 && rim[x]) || (y === S - 1 && rim[S + x]) || (x === 0 && rim[2 * S + y]) || (x === S - 1 && rim[3 * S + y])) ? 33 : 17;
      ground();
    }
    enter(OCT); workAll = Math.max(1, T * all * R0 * R0);           // rough, for the progress figure only
  }
  // run for about `ms` milliseconds; true when the level is finished
  function step(ms) {
    const t0 = performance.now();
    while (stage >= 0) {
      if (stageSweeps >= STAGE[stage][1] || n < 1) { if (stage === 0) { stage = -1; break; } enter(stage - 1); continue; }
      const left = STAGE[stage][1] - stageSweeps, part = Math.min(left, 0.125), batch = Math.max(1, Math.round(n * part));
      for (let k = 0; k < batch; k++) move();
      movesDone += batch; workDone += batch; stageSweeps += part;
      if (performance.now() - t0 >= ms) break;
    }
    return stage < 0;
  }
  function result() { const out = new Uint8Array(G * G); for (let p = 0; p < G * G; p++) if (g[p] & 1) out[p] = 255; return out; }
  return { start, step, result, G, T, R0, R1, OCT, limit: zoomLimit(G), rate: zoomRate,
           get joined() { return [joined, joinCells, missed]; }, get n() { return n; }, get progress() { return Math.min(1, workDone / workAll); }, get moves() { return movesDone; } };
}

// a picture at half the resolution: a site is occupied if any of the four it stands for is
function orHalf(pic, S) {
  const h = S >> 1, out = new Uint8Array(h * h);
  for (let y = 0; y < h; y++) for (let x = 0; x < h; x++) { const p = 2 * x + 2 * y * S; if (pic[p] | pic[p + 1] | pic[p + S] | pic[p + S + 1]) out[x + y * h] = 255; }
  return out;
}

// Going the other way: given a level, make the coarser levels around it.
//
// What lies outside a window cannot be had by running a local walk outwards; there is nothing
// there to start from.  It is built from the far end instead.  A far picture, some other
// typical piece of shape, stands for the surroundings at a coarser scale: the scale at which
// the given level is FW = 64 of its sites across.  Where in it should the given level go?  Where
// it looks most like the given level seen that coarsely.  Both are reduced to an 8 x 8 pattern
// of which blocks hold anything, and the far picture is searched, over positions and its eight
// orientations, for the 64 x 64 patch whose pattern differs from the given level's in the
// fewest blocks, with a charge as well for what the far picture has just outside each stretch
// of the patch's edge: heavy where the given level runs off its window there and the far
// picture has nothing for it to run on into, lighter where it has something and the given
// level does not reach the edge.  The far picture is then refined towards that patch J times as for zooming in,
// and at every step the middle of the new level is held to what is known there: the given
// level at that resolution (orHalf, applied the right number of times).  The walk runs around
// the held part and joins on to it, and stray islands are linked (Refiner, ground).  The
// result is J new levels, the outermost of which becomes the given level of the next batch.
//
// The far picture comes from a chain of refinements of its own, one more step each batch at a
// randomly chosen place, so it is never the same.
function Outward(G, seed, opt) {
  const R = Refiner(G, seed, opt), C = G >> 1, FW = 64, NB = 8, J = Math.round(Math.log2(G / FW)) - 1, LIM = zoomLimit(G);
  let s = (seed ^ 0x5bd1e995) >>> 0;
  function rnd() { s = (s + 0x6D2B79F5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }
  // In real samples the block just outside a stretch of edge holds something 99% of the time if the
  // window has cells on that stretch and 23% of the time if it has none; a far patch is charged
  // the improbability (minus the log) of what it has there, an inside block that differs being charged 1.
  const edge = (cells, outside) => cells ? (outside ? 0.01 : 4.5) : (outside ? 1.47 : 0.26);
  let exits = null;
  let far = null, pins = null, rims = null, want = null, jobs = [], cur = null, at = null, out = null, running = false, fit = [0, 0];
  function place(pic) {                                // an occupied site within reach of the middle, at random
    for (let t = 0; t < 2000; t++) { const x = C + (((rnd() * 2 - 1) * LIM) | 0), y = C + (((rnd() * 2 - 1) * LIM) | 0); if (pic[x + y * G]) return [x, y]; }
    return nearestCell(pic, G, C, C, LIM) || [C, C];
  }
  function turned(pic, k) {                            // one of the eight orientations
    const o = new Uint8Array(G * G);
    for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) {
      let u = k & 1 ? G - 1 - x : x, v = k & 2 ? G - 1 - y : y; if (k & 4) { const t = u; u = v; v = t; }
      o[x + y * G] = pic[u + v * G];
    }
    return o;
  }
  // the orientation and place in the far picture whose 8 x 8 pattern is nearest the wanted one
  function match() {
    const bs = FW / NB, P = G + 1, sat = new Int32Array(P * P); let best = 1e9, pick = null, ties = 0;
    for (let k = 0; k < 8; k++) {
      const pic = turned(far, k);
      for (let y = 0; y < G; y++) { let row = 0; for (let x = 0; x < G; x++) { row += pic[x + y * G] ? 1 : 0; sat[(x + 1) + (y + 1) * P] = sat[(x + 1) + y * P] + row; } }
      for (let fy = C - LIM; fy <= C + LIM; fy += 4) for (let fx = C - LIM; fx <= C + LIM; fx += 4) {
        let d = 0; const X0 = fx - FW / 2, Y0 = fy - FW / 2;
        const occ = (bx, by) => { const x0 = X0 + bx * bs, y0 = Y0 + by * bs; return sat[(x0 + bs) + (y0 + bs) * P] - sat[x0 + (y0 + bs) * P] - sat[(x0 + bs) + y0 * P] + sat[x0 + y0 * P] > 0; };
        for (let by = 0; by < NB && d <= best; by++) for (let bx = 0; bx < NB; bx++) if (occ(bx, by) !== !!want[bx + by * NB]) d++;
        // where the given level runs off its edge, the far picture should have something just outside to run on into
        for (let k2 = 0; k2 < NB && d <= best; k2++)
          d += edge(exits[k2], occ(k2, -1)) + edge(exits[NB + k2], occ(k2, NB)) + edge(exits[2 * NB + k2], occ(-1, k2)) + edge(exits[3 * NB + k2], occ(NB, k2));
        if (d < best) { best = d; ties = 1; pick = [k, fx, fy]; } else if (d === best && rnd() * ++ties < 1) pick = [k, fx, fy];
      }
    }
    fit = [best, NB * NB]; return pick;
  }
  // begin a batch outward from the level F (G x G bytes)
  function begin(F) {
    // the given level at each coarser resolution, and where along its edge it really has cells on its outermost row or column:
    // a coarse site on the pin's rim may be occupied only because of cells further in, and nothing should grow from it then
    pins = []; rims = []; let p = F, S = G, rm = new Uint8Array(4 * G);
    for (let t = 0; t < G; t++) { rm[t] = F[t]; rm[G + t] = F[t + (G - 1) * G]; rm[2 * G + t] = F[t * G]; rm[3 * G + t] = F[t * G + G - 1]; }
    for (let i = 1; i <= J; i++) {
      p = orHalf(p, S); const h = S >> 1, r2 = new Uint8Array(4 * h); for (let q = 0; q < 4; q++) for (let t = 0; t < h; t++) r2[q * h + t] = rm[q * S + 2 * t] | rm[q * S + 2 * t + 1];
      S = h; rm = r2; pins.push(p); rims.push(rm);
    }
    while (S > NB) { p = orHalf(p, S); S >>= 1; } want = p;
    // the stretches of each side (top, bottom, left, right; NB to a side) where F has cells on its outermost row or column
    exits = new Uint8Array(4 * NB); const seg = G / NB;
    for (let t = 0; t < G; t++) { const k2 = (t / seg) | 0; if (F[t]) exits[k2] = 1; if (F[t + (G - 1) * G]) exits[NB + k2] = 1; if (F[t * G]) exits[2 * NB + k2] = 1; if (F[t * G + G - 1]) exits[3 * NB + k2] = 1; }
    if (!far) far = F;
    jobs = ['far', 'match']; for (let i = J; i >= 1; i--) jobs.push(i);
    out = []; cur = null; running = false;
  }
  // run for about `ms` milliseconds; the J new levels, nearest first, once they are all made
  function step(ms) {
    const t0 = performance.now();
    while (jobs.length) {
      const job = jobs[0];
      if (job === 'match') { const m = match(); cur = turned(far, m[0]); at = [m[1], m[2]]; jobs.shift(); continue; }
      if (!running) {
        if (job === 'far') { const f = place(far); R.start(far, f[0], f[1]); } else { R.start(cur, at[0], at[1], pins[job - 1], G >> job, rims[job - 1]); at = [C, C]; }
        running = true;
      }
      if (!R.step(Math.max(1, ms - (performance.now() - t0)))) return null;
      running = false; jobs.shift();
      if (job === 'far') far = R.result(); else { cur = R.result(); out.unshift(cur); }
      if (performance.now() - t0 >= ms) break;
    }
    return jobs.length ? null : out;
  }
  return { begin, step, J, G, get joined() { return R.joined; }, get fit() { return fit; } };
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
if (typeof module !== 'undefined') module.exports = { Refiner, Outward, orHalf, nearestCell, ZOOM, zoomRate, zoomLimit };
