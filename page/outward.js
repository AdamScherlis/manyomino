// The exact zoom out.  A level is a "world": a whole polyomino of N cells that stands for the shape
// around the root at one scale, twice as coarse as the level before.  The next world is made from
// the last in two moves (src/outward.c, compiled to WebAssembly):
//   1. halve: every 2x2 block holding a cell becomes a cell, and a local walk thins the result
//      without leaving the finer shape: no cell further than one site from where the finer shape
//      has cells, and none of those places left without a cell on it or next to it;
//   2. extend: keep the box of half-width B about the root, and run the real sampler on all the
//      rest.  That draws the rest of an N-cell polyomino given the box, uniformly at random.
//      Strands that touch in the box after halving but were not joined there in the finer world
//      are not taken for joined: the sampler has to join them outside, as the finer world did.
// A world's cells are big on the screen, so its ring (from the box out to twice the box) is refined
// J times (a cell becomes four, then the local walk), each time with the finer level's picture held
// fixed inside the box and joined to it.
// Level 0 is real cells; no picture is ever refined below real cells (there a cell just becomes
// four).  Coordinates put the root at the origin; r is max(|x|, |y|) with -1 folded onto 0, so
// r < R is the square of side 2R about the root.

// B     half-width of the held box, in cells of the world          N     cells in a world
// S     moves of the sampler per cell                              TC    sweeps of the walk that thins a halved world
// RULE  3: that walk adds cells only on or next to the halved shape   COVER 1: and leaves none of its sites bare
// W     weight of a cell in that walk, times 1/lambda             PARTS 5: parts of the finer world are kept apart in the box (0: not)
// TINY  a part of this many cells or fewer is merged with what it touches
// T     sweeps of the walk at each stage of a picture; T1 at the first stage, halving from there down to T
// COLLAR  depth of the strip inside the seam that the walk may rework
// NEAR  a join this short is made wherever ring and inside have come apart; REACH the longest join, made only if needed
// CARRY 1: a strand grown in over the seam stays with the ring if it is all that joins two things
const OUT = { B: 8, N: 16000, S: 100, TC: 16, RULE: 3, COVER: 1, W: 0.65, PARTS: 5, TINY: 3, T: 16, T1: 32, COLLAR: 4, NEAR: 2, REACH: 32, CARRY: 1 };
const outRad = (x, y) => { const a = x < 0 ? -x - 1 : x, b = y < 0 ? -y - 1 : y; return a > b ? a : b; };

// E: the engine's exports (after E.init); o: { J, seed, and any of OUT's fields }
function Outgoer(E, o) {
  const P = Object.assign({}, OUT, o), B = P.B, J = P.J, G = 4 * B << J;
  let world = null, prev = null, level = 0;
  E.seed(P.seed >>> 0, 0x5eed);
  const put = (at, x, y) => { new Int32Array(E.memory.buffer, E.xs() + 4 * at, x.length).set(x); new Int32Array(E.memory.buffer, E.ys() + 4 * at, y.length).set(y); };
  const get = (from, to) => ({ x: new Int32Array(E.memory.buffer, E.xs() + 4 * from, to - from).slice(), y: new Int32Array(E.memory.buffer, E.ys() + 4 * from, to - from).slice() });
  const keep = (c, f) => { const n = c.x.length, x = new Int32Array(n), y = new Int32Array(n); let m = 0; for (let i = 0; i < n; i++) if (f(outRad(c.x[i], c.y[i]))) { x[m] = c.x[i]; y[m++] = c.y[i]; } return { x: x.subarray(0, m), y: y.subarray(0, m) }; };
  const keep2 = (c, f) => { const n = c.x.length, x = new Int32Array(n), y = new Int32Array(n); let m = 0; for (let i = 0; i < n; i++) if (f(c.x[i], c.y[i])) { x[m] = c.x[i]; y[m++] = c.y[i]; } return { x: x.subarray(0, m), y: y.subarray(0, m) }; };
  const twice = c => { const n = c.x.length, x = new Int32Array(4 * n), y = new Int32Array(4 * n); for (let i = 0; i < n; i++) for (let k = 0; k < 4; k++) { x[4 * i + k] = 2 * c.x[i] + (k & 1); y[4 * i + k] = 2 * c.y[i] + (k >> 1); } return { x, y }; };
  const join = (a, b) => { const x = new Int32Array(a.x.length + b.x.length), y = new Int32Array(x.length); x.set(a.x); x.set(b.x, a.x.length); y.set(a.y); y.set(b.y, a.y.length); return { x, y }; };
  let rs = (P.seed ^ 0x9e3779b9) >>> 0; const rnd = () => { rs ^= rs << 13; rs >>>= 0; rs ^= rs >>> 17; rs ^= rs << 5; rs >>>= 0; return rs / 4294967296; };
  const KEY = 1 << 15, key = (x, y) => (y + KEY) * 2 * KEY + (x + KEY);

  // level 0: real cells, root at the origin
  function begin(cells) {
    world = { x: Int32Array.from(cells.x), y: Int32Array.from(cells.y), lab: null }; level = 0;
    prev = [keep(world, r => r < 2 * B)]; for (let j = 1; j <= J; j++) prev.push(twice(prev[j - 1]));
    return picture(prev[J]);
  }
  function picture(c) { const pic = new Uint8Array(G * G), h = G >> 1; for (let i = 0; i < c.x.length; i++) pic[(c.y[i] + h) * G + c.x[i] + h] = 255; return pic; }

  // Which cells of the world, out to twice the box, are joined to which.  Cells of the box that touch are joined
  // only if they carry the same label (see labels()); all other touching cells are joined.
  const NB = [[1, 0], [0, 1], [-1, 0], [0, -1]];
  function fineParts() {
    const pos = new Map(), idx = []; for (let i = 0; i < world.x.length; i++) if (outRad(world.x[i], world.y[i]) < 2 * B) { pos.set(key(world.x[i], world.y[i]), idx.length); idx.push(i); }
    const part = new Int32Array(idx.length).fill(-1), L = world.lab; let np = 0;
    for (let k = 0; k < idx.length; k++) {
      if (part[k] >= 0) continue; part[k] = np; const st = [k];
      while (st.length) { const a = st.pop(), i = idx[a]; for (const [dx, dy] of NB) { const b = pos.get(key(world.x[i] + dx, world.y[i] + dy)); if (b === undefined || part[b] >= 0) continue; const j = idx[b]; if (L && L[i] && L[j] && L[i] !== L[j]) continue; part[b] = np; st.push(b); } }
      np++;
    }
    return { pos, part, np };
  }
  // Labels for the cells of the halved world that lie in the box: the part of the finer shape each stands for.
  // Two parts that touch after halving without being joined in the finer shape keep different labels, and neither
  // the thinning walk nor the sampler takes them for joined; so the sampler has to join them outside the box, as
  // the finer shape did.  Where one cell stands for two parts they become one.
  function labels(Y, fp) {
    const m = Y.x.length, lab = new Int32Array(m), up = Array.from({ length: fp.np }, (_, a) => a);
    const find = a => { while (up[a] !== a) a = up[a] = up[up[a]]; return a; };
    for (let i = 0; i < m; i++) { if (outRad(Y.x[i], Y.y[i]) >= B) continue; lab[i] = -1;
      for (let q = 0; q < 4; q++) { const k = fp.pos.get(key(2 * Y.x[i] + (q & 1), 2 * Y.y[i] + (q >> 1))); if (k === undefined) continue; if (lab[i] < 0) lab[i] = fp.part[k]; else up[find(fp.part[k])] = find(lab[i]); } }
    // a part of only a cell or two is not kept apart from what it touches
    if (P.TINY) { const size = new Int32Array(fp.np), at = new Map(); for (let k = 0; k < fp.part.length; k++) size[fp.part[k]]++;
      for (let i = 0; i < m; i++) if (outRad(Y.x[i], Y.y[i]) < B) at.set(key(Y.x[i], Y.y[i]), i);
      const sz = a => { let t = 0; for (let q = 0; q < fp.np; q++) if (find(q) === find(a)) t += size[q]; return t; };
      for (const [k, i] of at) { if (sz(lab[i]) > P.TINY) continue; for (const [dx, dy] of NB) { const j = at.get(key(Y.x[i] + dx, Y.y[i] + dy)); if (j !== undefined && find(lab[j]) !== find(lab[i])) { up[find(lab[i])] = find(lab[j]); break; } } } }
    for (let i = 0; i < m; i++) if (outRad(Y.x[i], Y.y[i]) < B) lab[i] = find(lab[i]) + 1;
    let kinds = 0; for (let a = 0; a < fp.np; a++) if (find(a) === a) kinds++;
    return { lab, of: p_ => find(p_) + 1, kinds };
  }
  // Halve the world, in the engine: every 2x2 block holding a cell becomes a cell, and a local walk thins the result
  // while staying on or next to the finer shape.  Leaves the halved world in the engine, labelled, the box held.
  function halve(info) {
    const fp = fineParts();
    E.clear(); put(0, world.x, world.y); E.place(world.x.length);
    info.halved = E.coarsen(); const lb = labels(get(0, E.count()), fp); info.partsHalved = lb.kinds;
    E.hold(-B, -B, B, B, 1);
    if (P.PARTS) {
      new Int32Array(E.memory.buffer, E.partlabels(), lb.lab.length).set(lb.lab); E.parts(P.PARTS);
      const ol = new Int32Array(E.memory.buffer, E.orlabels(), 4 * B * B).fill(0), Y = get(0, E.count());
      for (let i = 0; i < lb.lab.length; i++) if (lb.lab[i]) ol[(Y.y[i] + B) * 2 * B + Y.x[i] + B] = lb.lab[i];
    }
    E.tube(P.RULE); E.shadow(P.COVER, P.W); E.relax(4 * P.TC); E.tube(0); E.shadow(0, 1); if (P.PARTS) E.parts(1); info.thinned = E.count(); info.fineParts = fp.np;
    info.partsThinned = new Set(new Int32Array(E.memory.buffer, E.partlabels(), E.count()).filter(v => v)).size;
    return Int32Array.from(fp.part, lb.of);            // for each cell of the finer world out to twice the box, the label of its part
  }

  // The connected pieces of a set of cells, all at r < R: for each cell, the number (from 1) of its piece.
  function pieces(c, R) {
    const W = 2 * R, g = new Int32Array(W * W), n = c.x.length, id = new Int32Array(n), st = new Int32Array(n); let np = 0;
    for (let i = 0; i < n; i++) g[(c.y[i] + R) * W + c.x[i] + R] = i + 1;
    for (let i = 0; i < n; i++) {
      if (id[i]) continue; id[i] = ++np; let top = 0; st[top++] = i;
      while (top) { const a = st[--top], x = c.x[a] + R, y = c.y[a] + R;
        for (let d = 0; d < 4; d++) { const u = x + NB[d][0], v = y + NB[d][1]; if (u < 0 || v < 0 || u >= W || v >= W) continue; const k = g[v * W + u] - 1; if (k >= 0 && !id[k]) { id[k] = np; st[top++] = k; } } }
    }
    return { id, np, at: (x, y) => (x < -R || y < -R || x >= R || y >= R) ? -1 : g[(y + R) * W + x + R] - 1 };
  }

  // The next level out, as a piece of work that can be done a little at a time (it yields between sweeps); its
  // value at the end is the picture (G x G bytes, the root at the middle), with .info.
  function* work() {
    const info = { level: level + 1 }, t0 = Date.now(), pinLab = halve(info);
    E.grow(P.N); E.pad(P.N); yield;
    for (let s = 0; s < P.S; s++) { E.run(E.count()); yield; }
    world = get(0, E.count()); world.lab = P.PARTS ? new Int32Array(E.memory.buffer, E.partlabels(), world.x.length).slice() : null; E.parts(0); level++; info.msWorld = Date.now() - t0;
    info.parts = world.lab ? new Set(world.lab.filter(v => v)).size : 0;
    // The picture.  Stage 0 is the world's own cells out to twice the box; each stage after is twice as fine.
    // The ring (box to twice the box) is refined as a closed region: nothing beyond it counts, and where a strand
    // leaves it one cell on the outer edge stays put.  Inside it is the finer level's picture, held fixed except
    // for a collar a few cells deep just inside the seam, which the walk may rework along with the ring.
    const st = [keep(world, r => r < 2 * B)], all = new Set(); for (let i = 0; i < world.x.length; i++) all.add(key(world.x[i], world.y[i]));
    let ring = keep(world, r => r >= B && r < 2 * B), inner = keep(world, r => r < B), anchors = [];
    for (let i = 0; i < ring.x.length; i++) { const x = ring.x[i], y = ring.y[i]; if (outRad(x, y) === 2 * B - 1) for (const [dx, dy] of NB) if (outRad(x + dx, y + dy) === 2 * B && all.has(key(x + dx, y + dy))) anchors.push(x, y, dx, dy); }
    // labels: of the inside's cells at the last stage (to start with, the box's own cells), and of the finer picture's cells at the stage before
    let innerLab = new Int32Array(inner.x.length), baseLab = world.lab ? pinLab : new Int32Array(prev[0].x.length).fill(1);
    { let k = 0; for (let i = 0; i < world.x.length; i++) if (outRad(world.x[i], world.y[i]) < B) innerLab[k++] = world.lab ? world.lab[i] : 1; }
    info.links = 0; info.unjoined = 0;
    for (let j = 1; j <= J; j++) {
      const b = B << j, fine = prev[j - 1];
      if (j > level) { st.push(twice(st[j - 1])); continue; }      // finer than real cells: a cell just becomes four
      // Where the ring touched the inside at the last stage it must touch it at this one, and the same piece of it.
      // Every cell of the inside carries the label of its piece; the finer picture's pieces take the labels of
      // the pieces they refine (found by looking each cell up in the picture it came from).
      let fineLab = baseLab, trans = null;
      if (j > 1) {
        const pc = pieces(prev[j - 2], b >> 1), pf = pieces(fine, b), vote = new Map(), best = new Int32Array(pf.np + 1), score = new Int32Array(pf.np + 1);
        for (let i = 0; i < fine.x.length; i++) { const k = pc.at(fine.x[i] >> 1, fine.y[i] >> 1); if (k < 0 || !baseLab[k]) continue; const v = pf.id[i] * 1048576 + baseLab[k], c = (vote.get(v) || 0) + 1; vote.set(v, c); if (c > score[pf.id[i]]) { score[pf.id[i]] = c; best[pf.id[i]] = baseLab[k]; } }
        for (let i = 0; i < fine.x.length; i++) {                 // a small piece may have no cell over where it came from: look next door
          if (best[pf.id[i]]) continue;
          for (let d = 0; d < 9 && !best[pf.id[i]]; d++) { const k = pc.at((fine.x[i] >> 1) + d % 3 - 1, (fine.y[i] >> 1) + (d / 3 | 0) - 1); if (k >= 0 && baseLab[k]) best[pf.id[i]] = baseLab[k]; }
        }
        let fresh = 0; for (let k = 0; k < baseLab.length; k++) if (baseLab[k] > fresh) fresh = baseLab[k];
        for (let q = 1; q <= pf.np; q++) if (!best[q]) best[q] = ++fresh;
        fineLab = Int32Array.from(pf.id, q => best[q]);
        // two pieces that only touched before may be one piece now, under one label: an old label leads to the piece that holds most of it
        const own = new Map(); for (const [v, c] of vote) { const lab = v % 1048576, o = own.get(lab); if (!o || c > o[1]) own.set(lab, [best[Math.floor(v / 1048576)], c]); }
        for (let q = 1; q <= pf.np; q++) if (!own.has(best[q])) own.set(best[q], [best[q], 0]);
        trans = lab => { const o = own.get(lab); return o ? o[0] : 0; };
      }
      const known = new Set(fineLab), touch = [], ig = pieces(inner, b >> 1);
      for (let i = 0; i < ring.x.length; i++) { const x = ring.x[i], y = ring.y[i]; for (const [dx, dy] of NB) { const k = ig.at(x + dx, y + dy); if (k >= 0) { const w = trans ? trans(innerLab[k]) : innerLab[k]; touch.push(x, y, dx, dy, known.has(w) ? w : 0); } } }
      // on the outer edge, one of the two cells facing out of each cell that led out stays where it is
      const held = new Set(), next = [];
      for (let k = 0; k < anchors.length; k += 4) {
        const dx = anchors[k + 2], dy = anchors[k + 3], t = rnd() < 0.5 ? 0 : 1, x = dx ? 2 * anchors[k] + (dx > 0 ? 1 : 0) : 2 * anchors[k] + t, y = dy ? 2 * anchors[k + 1] + (dy > 0 ? 1 : 0) : 2 * anchors[k + 1] + t;
        next.push(x, y, dx, dy); held.add(key(x, y));
      }
      // the cells in order: the finer picture short of the collar, the held cells of the ring, the collar, the rest of the ring
      // (a strand of the ring that has grown in over the seam is still the ring's; where its cells would land on the finer picture they are dropped)
      const taken = new Uint8Array(4 * b * b); for (let i = 0; i < fine.x.length; i++) taken[(fine.y[i] + b) * 2 * b + fine.x[i] + b] = 1;
      const r2 = keep2(twice(ring), (x, y) => outRad(x, y) >= b || !taken[(y + b) * 2 * b + x + b]), X_ = [], Y_ = [], L_ = []; let nfix = 0, maxLab = 0;
      for (let i = 0; i < fine.x.length; i++) if (outRad(fine.x[i], fine.y[i]) < b - P.COLLAR) { X_.push(fine.x[i]); Y_.push(fine.y[i]); L_.push(fineLab[i]); }
      const ncore = X_.length;
      for (let i = 0; i < r2.x.length; i++) if (held.has(key(r2.x[i], r2.y[i]))) { X_.push(r2.x[i]); Y_.push(r2.y[i]); L_.push(0); }
      nfix = X_.length;
      for (let i = 0; i < fine.x.length; i++) if (outRad(fine.x[i], fine.y[i]) >= b - P.COLLAR) { X_.push(fine.x[i]); Y_.push(fine.y[i]); L_.push(fineLab[i]); }
      for (let i = 0; i < r2.x.length; i++) if (!held.has(key(r2.x[i], r2.y[i]))) { X_.push(r2.x[i]); Y_.push(r2.y[i]); L_.push(0); }
      E.clear(); put(0, X_, Y_); E.place(X_.length);
      for (let i = 0; i < L_.length; i++) if (L_[i] > maxLab) maxLab = L_[i];
      new Int32Array(E.memory.buffer, E.partlabels(), L_.length).set(L_); E.hold(-b, -b, b, b, 1); E.parts(4);
      // Joins.  Pieces of the finer picture that carry one label (they stand for one part of the world) are joined to
      // each other by the shortest way through empty sites.  Then every place where the ring touched the inside gets
      // a short join if it has come apart; a longer one is made only where the ring's cell and the piece it touched
      // would otherwise not be joined at all.
      { const pa = pieces({ x: X_, y: Y_ }, 2 * b), up = Array.from({ length: pa.np + 1 }, (_, a) => a), find = a => { while (up[a] !== a) a = up[a] = up[up[a]]; return a; }, of = new Map();
        const size = new Int32Array(pa.np + 1), lists = new Map();
        for (let i = 0; i < L_.length; i++) { if (!L_[i]) continue; const q = pa.id[i]; size[q]++; let l = lists.get(L_[i]); if (!l) lists.set(L_[i], l = []); if (!l.includes(q)) l.push(q); }
        for (const [lab, l] of lists) {
          l.sort((u, v) => size[v] - size[u]); of.set(lab, l[0]);
          for (let t = 1; t < l.length; t++) {
            // breadth first from the smaller piece through empty sites inside, to any cell joined to the largest
            const q = l[t], seen = new Map(), qx = [], qy = []; let h = 0, hit = -1;
            if (find(q) === find(l[0])) continue;
            for (let k = 0; k < L_.length; k++) if (pa.id[k] === q) { qx.push(X_[k]); qy.push(Y_[k]); seen.set(key(X_[k], Y_[k]), -1); }
            while (h < qx.length && hit < 0 && qx.length < 60000) { const x = qx[h], y = qy[h];
              for (const [dx, dy] of NB) { const u = x + dx, v = y + dy, kk = key(u, v); if (seen.has(kk) || outRad(u, v) >= b) continue; const c = pa.at(u, v);
                if (c >= 0) { if (find(pa.id[c]) === find(l[0])) { hit = h; break; } continue; }
                seen.set(kk, h); qx.push(u); qy.push(v); }
              h++; }
            if (hit >= 0) { for (let w = hit; w >= 0 && seen.get(key(qx[w], qy[w])) !== -1; w = seen.get(key(qx[w], qy[w]))) info.links += E.add(qx[w], qy[w], lab); up[find(q)] = find(l[0]); }
          }
        }
        for (const [reach, all] of [[P.NEAR, true], [6, false], [P.REACH, false]]) for (let k = 0; k < touch.length; k += 5) {
          const want = touch[k + 4]; if (!want || !of.has(want)) continue;
          const cR = pa.id[pa.at(2 * touch[k], 2 * touch[k + 1])], cI = of.get(want); if (!all && find(cR) === find(cI)) continue;
          const a = E.link(touch[k], touch[k + 1], touch[k + 2], touch[k + 3], b, reach, want);
          if (a >= 0) { info.links += a; up[find(cR)] = find(cI); } else if (reach === P.REACH) info.unjoined++;
        } }
      E.pin(nfix, b - P.COLLAR, 2 * b); for (let q = 4 * Math.max(P.T, P.T1 >> (j - 1)); q > 0; q--) { E.relax(1); yield; }
      const got = get(0, E.count()), gl = new Int32Array(E.memory.buffer, E.partlabels(), E.count()); st.push(got); E.parts(0);
      // The ring from here on is what lies outside the seam.  A strand of it that grew in over the seam (its label is the
      // new one) is left behind with the inside, unless it is all that joins two things: two parts of the ring, or the
      // ring and a piece of the inside that it does not touch elsewhere.  Such a strand stays with the ring.
      { const n = got.x.length, kind = new Uint8Array(n), ga = pieces(got, 2 * b);                  // 1 ring, 2 inside, 3 a strand grown in
        for (let i = 0; i < n; i++) kind[i] = outRad(got.x[i], got.y[i]) >= b ? 1 : gl[i] > 0 && gl[i] <= maxLab ? 2 : 3;
        const comp = new Int32Array(n), stack = new Int32Array(n); let nc = 0;                      // pieces of the ring alone, and of the strands alone
        for (let i = 0; i < n; i++) { if (comp[i] || kind[i] === 2) continue; comp[i] = ++nc; let top = 0; stack[top++] = i;
          while (top) { const a = stack[--top]; for (const [dx, dy] of NB) { const k = ga.at(got.x[a] + dx, got.y[a] + dy); if (k >= 0 && !comp[k] && kind[k] === kind[i]) { comp[k] = nc; stack[top++] = k; } } } }
        const up = Array.from({ length: nc + maxLab + 2 }, (_, a) => a), find = a => { while (up[a] !== a) a = up[a] = up[up[a]]; return a; }, node = k => kind[k] === 2 ? nc + gl[k] : comp[k];
        for (let i = 0; i < n; i++) if (kind[i] === 1) for (const [dx, dy] of NB) { const k = ga.at(got.x[i] + dx, got.y[i] + dy); if (k >= 0 && kind[k] === 2) up[find(comp[i])] = find(nc + gl[k]); }
        const touches = new Map();                                                                   // for each strand, what it touches
        for (let i = 0; i < n; i++) if (kind[i] === 3) for (const [dx, dy] of NB) { const k = ga.at(got.x[i] + dx, got.y[i] + dy); if (k >= 0 && kind[k] !== 3) { let t = touches.get(comp[i]); if (!t) touches.set(comp[i], t = new Set()); t.add(node(k)); } }
        const stays = new Uint8Array(nc + 1);
        if (P.CARRY) for (const [c, t] of touches) { const cls = new Set([...t].map(find)); if (cls.size < 2) continue; stays[c] = 1; const f = [...cls]; for (let q = 1; q < f.length; q++) up[find(f[q])] = find(f[0]); }
        const rx = [], ry = [], ix = [], iy = [], il = [];
        for (let i = 0; i < n; i++) if (kind[i] === 1 || (kind[i] === 3 && stays[comp[i]])) { rx.push(got.x[i]); ry.push(got.y[i]); } else if (kind[i] === 2) { ix.push(got.x[i]); iy.push(got.y[i]); il.push(gl[i]); }
        ring = { x: Int32Array.from(rx), y: Int32Array.from(ry) }; inner = { x: Int32Array.from(ix), y: Int32Array.from(iy) }; innerLab = Int32Array.from(il); }
      anchors = next; baseLab = fineLab;
    }
    prev = st; info.ms = Date.now() - t0;
    const pic = picture(st[J]); pic.info = info; return pic;
  }
  function next() { const g = work(); for (;;) { const r = g.next(); if (r.done) return r.value; } }
  // the same, working for about `ms` milliseconds at a call: null until the level is finished, then its picture
  let going = null;
  function step(ms) { const t = Date.now(); if (!going) going = work(); for (;;) { const r = going.next(); if (r.done) { going = null; return r.value; } if (Date.now() - t >= ms) return null; } }
  return { begin, next, step, halve, G, get level() { return level; }, get world() { return world; }, get stages() { return prev; } };
}

if (typeof module !== 'undefined') module.exports = { OUT, outRad, Outgoer };
