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
// fixed inside the box, and together with a strip of the same world beyond the ring, so that the
// ring's outer rows are never refined against an edge.  The strip is not kept, but its first rows
// are what the next level's ring starts from beside this picture, at every stage.
// Level 0 is real cells; no picture is ever refined below real cells (there a cell just becomes
// four).  Coordinates put the root at the origin; r is max(|x|, |y|) with -1 folded onto 0, so
// r < R is the square of side 2R about the root.

// B     half-width of the held box, in cells of the world          N     cells in a world
// S     moves of the sampler per cell                              TC    sweeps of the walk that thins a halved world
// RULE  3: that walk adds cells only on or next to the halved shape   COVER 1: and leaves none of its sites bare
// W     weight of a cell in that walk, times 1/lambda             PARTS 5: parts of the finer world are kept apart in the box (0: not)
// TINY  a part of this many cells or fewer is merged with what it touches
// HOLD  layers outside the box that the sampler holds as well
// T     sweeps of the walk at each stage of a picture; T1 at the first stage, halving from there down to T; TL at the last stage
// WR    weight of a cell in the walk at a stage before the last, times 1/lambda; WL at the last stage
// COLLAR  depth of the strip inside the seam that the walk may rework at the last stage (see outCollar)
// MARGIN  layers of the world beyond twice the box that are refined along with the ring
// BAND  rows outside the seam that are taken, at each stage, from the finer level's strip
const OUT = { B: 8, N: 16000, S: 100, TC: 16, RULE: 3, COVER: 1, W: 0.72, HOLD: 2, TINY: 3, T: 16, T1: 32, TL: 64, WR: 0.985, WL: 0.9985, COLLAR: 12, MARGIN: 2, BAND: 2 };
// the depth, in cells of stage j, of the strip inside the seam that the walk may rework: COLLAR, but never more than a quarter of the way in
const outCollar = (P, j) => Math.min(P.COLLAR, (P.B << j) >> 2);
const outRad = (x, y) => { const a = x < 0 ? -x - 1 : x, b = y < 0 ? -y - 1 : y; return a > b ? a : b; };

// E: the engine's exports (after E.init); o: { J, seed, and any of OUT's fields }
function Outgoer(E, o) {
  const P = Object.assign({}, OUT, o), B = P.B, J = P.J, G = 4 * B << J;
  let world = null, prev = null, prevLab = null, prevBand = null, crossAt = null, level = 0;      // crossAt[j]: where the picture leaves its outer edge at stage j
  E.seed(P.seed >>> 0, 0x5eed);
  const put = (at, x, y) => { new Int32Array(E.memory.buffer, E.xs() + 4 * at, x.length).set(x); new Int32Array(E.memory.buffer, E.ys() + 4 * at, y.length).set(y); };
  const get = (from, to) => ({ x: new Int32Array(E.memory.buffer, E.xs() + 4 * from, to - from).slice(), y: new Int32Array(E.memory.buffer, E.ys() + 4 * from, to - from).slice() });
  const keep = (c, f) => { const n = c.x.length, x = new Int32Array(n), y = new Int32Array(n); let m = 0; for (let i = 0; i < n; i++) if (f(outRad(c.x[i], c.y[i]))) { x[m] = c.x[i]; y[m++] = c.y[i]; } return { x: x.subarray(0, m), y: y.subarray(0, m) }; };
  const keep2 = (c, f) => { const n = c.x.length, x = new Int32Array(n), y = new Int32Array(n); let m = 0; for (let i = 0; i < n; i++) if (f(c.x[i], c.y[i])) { x[m] = c.x[i]; y[m++] = c.y[i]; } return { x: x.subarray(0, m), y: y.subarray(0, m) }; };
  const twice = c => { const n = c.x.length, x = new Int32Array(4 * n), y = new Int32Array(4 * n); for (let i = 0; i < n; i++) for (let k = 0; k < 4; k++) { x[4 * i + k] = 2 * c.x[i] + (k & 1); y[4 * i + k] = 2 * c.y[i] + (k >> 1); } return { x, y }; };
  let rs = (P.seed ^ 0x9e3779b9) >>> 0; const rnd = () => { rs ^= rs << 13; rs >>>= 0; rs ^= rs >>> 17; rs ^= rs << 5; rs >>>= 0; return rs / 4294967296; };
  const KEY = 1 << 15, key = (x, y) => (y + KEY) * 2 * KEY + (x + KEY);

  // level 0: real cells, root at the origin
  function begin(cells) {
    world = { x: Int32Array.from(cells.x), y: Int32Array.from(cells.y), lab: null }; level = 0;
    prev = [keep(world, r => r < 2 * B)]; for (let j = 1; j <= J; j++) prev.push(twice(prev[j - 1]));
    crossAt = [crossings(world)];
    prevLab = [pieces(prev[0], 2 * B, null).id]; prevBand = [keep(world, r => r >= 2 * B && r < 2 * B + P.BAND)];
    for (let j = 1; j <= J; j++) { prevLab.push(Int32Array.from({ length: 4 * prevLab[j - 1].length }, (_, i) => prevLab[j - 1][i >> 2])); prevBand.push(keep(twice(prevBand[j - 1]), r => r < (2 * B << j) + P.BAND)); }
    return picture(prev[J]);
  }
  // Where the world crosses the edge of its picture (r = 2B): the cell just inside, and the way out to the cell just outside.
  function crossings(w) { const all = new Set(), out = []; for (let i = 0; i < w.x.length; i++) all.add(key(w.x[i], w.y[i]));
    for (let i = 0; i < w.x.length; i++) { const x = w.x[i], y = w.y[i]; if (outRad(x, y) === 2 * B - 1) for (const [dx, dy] of NB) if (outRad(x + dx, y + dy) === 2 * B && all.has(key(x + dx, y + dy))) out.push(x, y, dx, dy); }
    return out; }
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
  let marks = [];                                     // the marked sites of the present world, x and y in turn
  function halve(info) {
    const fp = fineParts();
    E.clear(); put(0, world.x, world.y); E.place(world.x.length);
    info.halved = E.coarsen(); const lb = labels(get(0, E.count()), fp); info.partsHalved = lb.kinds;
    // A cell of the box and a cell outside it are joined only where the finer world crosses the edge of the box
    // between them: those sites outside are marked.  Anywhere else the two only touch, as two parts in the box may.
    E.unmark(); marks = []; for (let k = 0; k < crossAt[0].length; k += 4) { const u = (crossAt[0][k] + crossAt[0][k + 2]) >> 1, v = (crossAt[0][k + 1] + crossAt[0][k + 3]) >> 1; E.flag(u, v); marks.push(u, v); }
    E.hold(-B, -B, B, B, 1);
    new Int32Array(E.memory.buffer, E.partlabels(), lb.lab.length).set(lb.lab); E.parts(5); E.strict(1);
    const ol = new Int32Array(E.memory.buffer, E.orlabels(), 4 * B * B).fill(0), Y = get(0, E.count());
    for (let i = 0; i < lb.lab.length; i++) if (lb.lab[i]) ol[(Y.y[i] + B) * 2 * B + Y.x[i] + B] = lb.lab[i];
    E.tube(P.RULE); E.shadow(P.COVER, P.W); E.relax(4 * P.TC); E.tube(0); E.shadow(0, 1); E.parts(1); info.thinned = E.count(); info.fineParts = fp.np;
    return Int32Array.from(fp.part, lb.of);            // for each cell of the finer world out to twice the box, the label of its part
  }

  // The connected pieces of a set of cells, all at r < R: for each cell, the number (from 1) of its piece.
  // Given labels, two touching cells that both have one are joined only if the labels are the same.
  function pieces(c, R, lab) {
    const W = 2 * R, g = new Int32Array(W * W), n = c.x.length, id = new Int32Array(n), st = new Int32Array(n); let np = 0;
    for (let i = 0; i < n; i++) g[(c.y[i] + R) * W + c.x[i] + R] = i + 1;
    for (let i = 0; i < n; i++) {
      if (id[i]) continue; id[i] = ++np; let top = 0; st[top++] = i;
      while (top) { const a = st[--top], x = c.x[a] + R, y = c.y[a] + R;
        for (let d = 0; d < 4; d++) { const u = x + NB[d][0], v = y + NB[d][1]; if (u < 0 || v < 0 || u >= W || v >= W) continue; const k = g[v * W + u] - 1; if (k >= 0 && !id[k] && !(lab && lab[a] && lab[k] && lab[a] !== lab[k])) { id[k] = np; st[top++] = k; } } }
    }
    return { id, np, at: (x, y) => (x < -R || y < -R || x >= R || y >= R) ? -1 : g[(y + R) * W + x + R] - 1 };
  }

  // The next level out, as a piece of work that can be done a little at a time (it yields between sweeps); its
  // value at the end is the picture (G x G bytes, the root at the middle), with .info.
  function* work() {
    const info = { level: level + 1 }, t0 = Date.now(), pinLab = halve(info);
    if (P.onHalved) P.onHalved(get(0, E.count()));     // (for tests: the halved world before anything is drawn afresh)
    // The sampler holds the box and the first layers outside it: those are the finer world's own cells, halved, and
    // say where that world crosses the edge of the box, which the box alone does not.
    E.hold(-B - P.HOLD, -B - P.HOLD, B + P.HOLD, B + P.HOLD, 1);
    E.grow(P.N); E.pad(P.N); yield;
    for (let s = 0; s < P.S; s++) { E.run(E.count()); yield; }
    world = get(0, E.count()); world.lab = new Int32Array(E.memory.buffer, E.partlabels(), world.x.length).slice(); E.parts(0); level++; info.msWorld = Date.now() - t0;
    info.parts = new Set(world.lab.filter(v => v)).size;
    // The picture.  Stage 0 is the world's own cells; each stage after is twice as fine.  The ring (box to twice
    // the box) is refined together with a strip of the same world beyond it (MARGIN layers), so that its outer rows
    // have their own neighbours beside them at every stage.  Only what lies inside twice the box is kept.  Inside
    // the ring is the finer level's picture, held; and the ring's first rows (BAND of them) are, at each stage, the
    // very cells that picture had beside it there: the first rows of the finer level's strip at that stage.
    const M0 = P.MARGIN, BAND = P.BAND, st = [], stLab = [], band = [], cross = [crossings(world)];
    const withLab = (c, lab, f) => { const n = c.x.length, x = [], y = [], l = []; for (let i = 0; i < n; i++) if (f(outRad(c.x[i], c.y[i]))) { x.push(c.x[i]); y.push(c.y[i]); l.push(lab ? lab[i] : 0); } return { x: Int32Array.from(x), y: Int32Array.from(y), lab: Int32Array.from(l) }; };
    { const k0 = withLab(world, world.lab, r => r < 2 * B); st.push({ x: k0.x, y: k0.y }); stLab.push(pieces(k0, 2 * B, k0.lab).id); band.push(keep(world, r => r >= 2 * B && r < 2 * B + BAND)); }
    let ring = keep(world, r => r >= B && r < 2 * B + M0);
    // where the world leaves the strip, one cell on its outer edge stays put at every stage (an anchor), so that a
    // strand which is joined to the rest only beyond the strip cannot draw back out of it and be left on its own
    let anchors = []; { const e = 2 * B + M0, w0 = withLab(world, world.lab, r => r < e), pw = pieces(w0, e, w0.lab), all = new Set(); for (let i = 0; i < world.x.length; i++) all.add(key(world.x[i], world.y[i]));
      for (let i = 0; i < w0.x.length; i++) { const x = w0.x[i], y = w0.y[i]; if (outRad(x, y) === e - 1) for (const [dx, dy] of NB) if (outRad(x + dx, y + dy) === e && all.has(key(x + dx, y + dy))) anchors.push(x, y, dx, dy, pw.id[i]); } }
    info.links = 0; info.loose = 0; info.dropped = 0; info.unjoined = 0; info.long = 0; info.taken = 0;
    for (let j = 1; j <= J; j++) {
      const b = B << j, R = 2 * b + (M0 << j), fine = prev[j - 1], fineLab = prevLab[j - 1], last = j === Math.min(level, J), C = last ? outCollar(P, j) : 0;
      if (j > level) { st.push(twice(st[j - 1])); stLab.push(Int32Array.from({ length: 4 * stLab[j - 1].length }, (_, i) => stLab[j - 1][i >> 2])); band.push(keep(twice(band[j - 1]), r => r < 2 * b + BAND)); continue; }      // finer than real cells: a cell just becomes four
      const own = keep(twice(ring), r => r >= b + BAND), bandIn = prevBand[j - 1], held = new Map(), nextA = [];
      for (let k = 0; k < anchors.length; k += 5) { const dx = anchors[k + 2], dy = anchors[k + 3], t = rnd() < 0.5 ? 0 : 1, x = dx ? 2 * anchors[k] + (dx > 0 ? 1 : 0) : 2 * anchors[k] + t, y = dy ? 2 * anchors[k + 1] + (dy > 0 ? 1 : 0) : 2 * anchors[k + 1] + t;
        if (!held.has(key(x, y))) { nextA.push(x, y, dx, dy, anchors[k + 4]); held.set(key(x, y), anchors[k + 4]); } }
      anchors = nextA;
      // the cells in order: the finer picture short of the collar, and the anchors (fixed); the collar (the last stage only); the band; the rest
      const X_ = [], Y_ = [], L_ = [];
      for (let i = 0; i < fine.x.length; i++) if (outRad(fine.x[i], fine.y[i]) < b - C - 1) { X_.push(fine.x[i]); Y_.push(fine.y[i]); L_.push(fineLab[i]); }
      const ndeep = X_.length;                                   // these have nothing beside them that the walk could change
      for (let i = 0; i < fine.x.length; i++) if (outRad(fine.x[i], fine.y[i]) === b - C - 1) { X_.push(fine.x[i]); Y_.push(fine.y[i]); L_.push(fineLab[i]); }
      const nK = X_.length; for (let i = 0; i < own.x.length; i++) if (held.has(key(own.x[i], own.y[i]))) { X_.push(own.x[i]); Y_.push(own.y[i]); L_.push(0); }
      const nfix = X_.length;
      for (let i = 0; i < fine.x.length; i++) if (outRad(fine.x[i], fine.y[i]) >= b - C) { X_.push(fine.x[i]); Y_.push(fine.y[i]); L_.push(fineLab[i]); }
      const nIn = X_.length;
      for (let i = 0; i < bandIn.x.length; i++) { X_.push(bandIn.x[i]); Y_.push(bandIn.y[i]); L_.push(0); }
      for (let i = 0; i < own.x.length; i++) if (!held.has(key(own.x[i], own.y[i]))) { X_.push(own.x[i]); Y_.push(own.y[i]); L_.push(0); }
      // Behind the band the ring's own cells (the last stage's, each become four) may have lost what joined them to
      // it.  What is joined to what is put back as it was in the world: all pieces that hold anchors of one piece of
      // the world (within the strip) are joined to each other, after every piece that holds no anchor has been
      // joined to whatever is nearest.  A join is the shortest way through empty sites outside the seam; smallest pieces first.
      { const n0 = X_.length, pa = pieces({ x: X_, y: Y_ }, R, last ? null : L_), up = Array.from({ length: pa.np + 1 }, (_, a) => a), find = a => { while (up[a] !== a) a = up[a] = up[up[a]]; return a; };
        const size = new Int32Array(pa.np + 1), inner = new Uint8Array(pa.np + 1), mem = Array.from({ length: pa.np + 1 }, () => []), cls = new Int32Array(pa.np + 1);   // cls: the piece of the world a group belongs to (0: not known)
        const cu = new Map(), cfind = a => { let r = a; while (cu.has(r) && cu.get(r) !== r) r = cu.get(r); return r; };      // two pieces of the world that have come to touch count as one
        for (let i = 0; i < n0; i++) { const q = pa.id[i]; size[q]++; mem[q].push(i); if (i < nK || (i >= nfix && i < nIn)) inner[q] = 1;
          if (i >= nK && i < nfix) { const c = cfind(held.get(key(X_[i], Y_[i]))); if (!cls[q]) cls[q] = c; else if (cfind(cls[q]) !== c) cu.set(c, cfind(cls[q])); } }
        const classOf = g => cls[g] ? cfind(cls[g]) : 0;
        const added = new Map(), dead = new Uint8Array(n0);
        const compAt = (u, v) => { const c = pa.at(u, v); if (c >= 0 && !dead[c]) return pa.id[c]; const a = added.get(key(u, v)); return a === undefined ? 0 : a; };   // (a join may pass over a site whose cell was taken)
        const merge = (a, b_) => { a = find(a); b_ = find(b_); if (a === b_) return; const ca = classOf(a), cb = classOf(b_); up[a] = b_; if (!cb) cls[b_] = ca; else if (ca && ca !== cb) cu.set(ca, cb); };
        // join group q to the nearest cell of another group for which ok(group) holds
        const join = (q, ok, most) => {
          const seen = new Map(), qx = [], qy = []; let h = 0, hit = -1, other = 0; const root = find(q);
          for (let g = 1; g <= pa.np; g++) if (find(g) === root) for (const i of mem[g]) if (!dead[i] && outRad(X_[i], Y_[i]) >= b - 1) { qx.push(X_[i]); qy.push(Y_[i]); seen.set(key(X_[i], Y_[i]), -1); }
          const nsrc = qx.length;
          while (h < qx.length && hit < 0 && qx.length < nsrc + 2000000) { const x = qx[h], y = qy[h];
            for (const [dx, dy] of NB) { const u = x + dx, v = y + dy, kk = key(u, v), c = compAt(u, v);
              if (c) { if (h >= nsrc && find(c) !== root && ok(find(c))) { hit = h; other = c; break; } continue; }
              const r = outRad(u, v); if (seen.has(kk) || r < b || r >= R) continue;
              seen.set(kk, h); qx.push(u); qy.push(v); }
            h++; }
          if (hit < 0) return 0;
          { let len = 0; for (let w = hit; w >= nsrc; w = seen.get(key(qx[w], qy[w]))) len++; if (len > most) return -1; info.loose++; if (len >= 9) info.long++; }
          for (let w = hit; w >= nsrc; w = seen.get(key(qx[w], qy[w]))) { added.set(key(qx[w], qy[w]), q); X_.push(qx[w]); Y_.push(qy[w]); L_.push(0); info.links++;
            for (const [dx, dy] of NB) { const c = compAt(qx[w] + dx, qy[w] + dy); if (c) merge(c, q); } }
          merge(other, q);
          // The join has added cells; as many are taken from loose ends close by, so that joins leave the number of cells
          // as it was (the walk thins what it is given by a factor, to no target, so cells added here would stay).
          // (for each cell of the join, a cell at the same distance from the root, give or take a row, so that the count
          // in each row stays what it was; the nearest such cell)
          for (let w = hit; w >= nsrc; w = seen.get(key(qx[w], qy[w]))) { const px = qx[w], py = qy[w], pr = outRad(px, py); let got = false;
            for (let slack = 1; slack <= 4 && !got; slack *= 2) for (let D = 2; D <= 64 && !got; D *= 2)
              for (let v = py - D; v <= py + D && !got; v++) for (let u = px - D; u <= px + D && !got; u++) { const c = pa.at(u, v); if (c < nIn || dead[c] || Math.abs(outRad(u, v) - pr) > slack) continue;
                // a cell with one neighbour, or the corner of a full block of four with nothing else beside it, can go (not one that a join ends at)
                let deg = 0, sib = 0; for (const [dx, dy] of NB) if (compAt(u + dx, v + dy)) { if (added.has(key(u + dx, v + dy))) deg = 9; deg++; if (((u + dx) >> 1) === (u >> 1) && ((v + dy) >> 1) === (v >> 1)) sib++; }
                if (deg === 1 || (deg === 2 && sib === 2 && compAt(u ^ 1, v ^ 1))) { dead[c] = 1; info.taken++; got = true; } } }
          return 1; };
        const bySize = []; for (let q = 1; q <= pa.np; q++) bySize.push(q); bySize.sort((u, v) => size[u] - size[v]);
        // 1: the pieces with no anchor, to whatever is nearest.  (A small piece of the ring's own cells that would
        // need a long join, and that nothing has been joined to, is dropped instead.)
        const alone = q => { const r = find(q); for (let g = 1; g <= pa.np; g++) if (g !== q && find(g) === r) return false; return !added.size || ![...added.values()].some(a => find(a) === r); };
        for (const q of bySize) { if (classOf(find(q))) continue;
          let live = 0, all = 0; { const r = find(q); for (let g = 1; g <= pa.np; g++) if (find(g) === r) for (const i of mem[g]) if (!dead[i]) { all++; if (g === q) live++; } } if (!all) continue;
          const small = !inner[q] && alone(q), r = join(q, () => true, small ? Math.max(8, size[q] >> 1) : 1e9);
          if (r > 0) continue;
          if (small) { for (const i of mem[q]) dead[i] = 1; info.dropped += live; } else info.unjoined++; }
        // 2: then the pieces of one piece of the world, to each other: all but the largest, smallest first
        for (const tried = new Set(); ;) { const sz = new Map(), byClass = new Map(); let did = false;
          for (let g = 1; g <= pa.np; g++) { const r = find(g); sz.set(r, (sz.get(r) || 0) + size[g]); }
          for (const [r, n] of sz) { const c = classOf(r); if (!c) continue; if (!byClass.has(c)) byClass.set(c, []); byClass.get(c).push([n, r]); }
          for (const [, list] of byClass) { if (list.length < 2) continue; list.sort((u, v) => u[0] - v[0]);
            for (let k = 0; k < list.length - 1; k++) { const r = list[k][1]; if (find(r) !== r || tried.has(r)) continue; tried.add(r); const c = classOf(r);
              if (join(r, g => classOf(g) === c, 1e9) > 0) did = true; else info.unjoined++; } }
          if (!did) break; }
        if (info.dropped || info.taken) { let m = 0; for (let i = 0; i < X_.length; i++) if (i >= n0 || !dead[i]) { X_[m] = X_[i]; Y_[m] = Y_[i]; L_[m++] = L_[i]; } X_.length = Y_.length = L_.length = m; } }
      if (P.onStage) P.onStage(level, j, 'start', { x: X_, y: Y_ }, { x: own.x, y: own.y }, bandIn);
      E.clear(); put(0, X_, Y_); E.place(X_.length);
      // The walk.  Before the last stage the inside is fixed and keeps the labels of its pieces: two of them that
      // only touch are not joined, so the ring cannot come to rely on that touch; the ring is joined to whatever it
      // touches.  The last stage that is not finer than real cells sets the fine texture and gets more sweeps; there
      // a collar of the inside is free as well and nothing has a label: whatever touches at the end is joined.
      if (!last) { new Int32Array(E.memory.buffer, E.partlabels(), L_.length).set(L_); E.hold(-b, -b, b, b, 1); E.parts(4); }
      E.shadow(0, last ? P.WL : P.WR);                             // (the weight of a cell, in units of 1/lambda)
      E.pin(nfix, b - C, R); E.idle(ndeep);
      for (let q = 4 * (last ? Math.max(P.TL, P.T1 >> (j - 1)) : Math.max(P.T, P.T1 >> (j - 1))); q > 0; q--) { E.relax(1); yield; }
      const got = get(0, E.count()), gl = last ? null : new Int32Array(E.memory.buffer, E.partlabels(), E.count()).slice(); E.parts(0);
      if (P.onStage) P.onStage(level, j, 'end', got);
      const k1 = withLab(got, gl, r => r < 2 * b);
      st.push({ x: k1.x, y: k1.y }); stLab.push(pieces(k1, 2 * b, last ? null : k1.lab).id); band.push(keep(got, r => r >= 2 * b && r < 2 * b + BAND));
      ring = keep(got, r => r >= b);
    }
    prev = st; prevLab = stLab; prevBand = band; crossAt = cross; info.ms = Date.now() - t0;
    const pic = picture(st[J]); pic.info = info; return pic;
  }
  function next() { const g = work(); for (;;) { const r = g.next(); if (r.done) return r.value; } }
  // the same, working for about `ms` milliseconds at a call: null until the level is finished, then its picture
  let going = null;
  function step(ms) { const t = Date.now(); if (!going) going = work(); for (;;) { const r = going.next(); if (r.done) { going = null; return r.value; } if (Date.now() - t >= ms) return null; } }
  return { begin, next, step, halve, G, get level() { return level; }, get world() { return world; }, get stages() { return prev; }, get marks() { return marks; } };
}

if (typeof module !== 'undefined') module.exports = { OUT, outRad, outCollar, Outgoer };
