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

const OUT = { B: 8, N: 16000, S: 100, TC: 16, RULE: 3, COVER: 1, W: 0.7, PARTS: 1, T: 16, REACH: 8, ANCHOR: 2 };
const outRad = (x, y) => { const a = x < 0 ? -x - 1 : x, b = y < 0 ? -y - 1 : y; return a > b ? a : b; };

// E: the engine's exports (after E.init); o: { J, seed, and any of OUT's fields }
function Outgoer(E, o) {
  const P = Object.assign({}, OUT, o), B = P.B, J = P.J, G = 4 * B << J;
  let world = null, prev = null, level = 0;
  E.seed(P.seed >>> 0, 0x5eed);
  const put = (at, x, y) => { new Int32Array(E.memory.buffer, E.xs() + 4 * at, x.length).set(x); new Int32Array(E.memory.buffer, E.ys() + 4 * at, y.length).set(y); };
  const get = (from, to) => ({ x: new Int32Array(E.memory.buffer, E.xs() + 4 * from, to - from).slice(), y: new Int32Array(E.memory.buffer, E.ys() + 4 * from, to - from).slice() });
  const keep = (c, f) => { const n = c.x.length, x = new Int32Array(n), y = new Int32Array(n); let m = 0; for (let i = 0; i < n; i++) if (f(outRad(c.x[i], c.y[i]))) { x[m] = c.x[i]; y[m++] = c.y[i]; } return { x: x.subarray(0, m), y: y.subarray(0, m) }; };
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
      new Int32Array(E.memory.buffer, E.partlabels(), lb.lab.length).set(lb.lab); E.parts(1);
      const ol = new Int32Array(E.memory.buffer, E.orlabels(), 4 * B * B).fill(0), Y = get(0, E.count());
      for (let i = 0; i < lb.lab.length; i++) if (lb.lab[i]) ol[(Y.y[i] + B) * 2 * B + Y.x[i] + B] = lb.lab[i];
    }
    E.tube(P.RULE); E.shadow(P.COVER, P.W); E.relax(4 * P.TC); E.tube(0); E.shadow(0, 1); info.thinned = E.count(); info.fineParts = fp.np;
    info.partsThinned = new Set(new Int32Array(E.memory.buffer, E.partlabels(), E.count()).filter(v => v)).size;
    return Int32Array.from(fp.part, lb.of);            // for each cell of the finer world out to twice the box, the label of its part
  }

  // the next level out; returns its picture (G x G bytes, the root at the middle), with .info
  function next() {
    const info = { level: level + 1 }, t0 = Date.now(), pinLab = halve(info), known = new Set(pinLab);
    E.grow(P.N); E.pad(P.N);
    for (let s = 0; s < P.S; s++) E.run(E.count());
    world = get(0, E.count()); world.lab = P.PARTS ? new Int32Array(E.memory.buffer, E.partlabels(), world.x.length).slice() : null; E.parts(0); level++; info.msWorld = Date.now() - t0;
    info.parts = world.lab ? new Set(world.lab.filter(v => v)).size : 0;
    // The picture.  Stage 0 is the world's own cells out to twice the box; each stage after is twice as fine.
    // The ring (box to twice the box) is refined as a closed region: the finer level's picture is held fixed inside
    // it, nothing beyond it counts, and its cells on the outer edge that had a cell beyond them stay where they are.
    const st = [keep(world, r => r < 2 * B)], all = new Set(); for (let i = 0; i < world.x.length; i++) all.add(key(world.x[i], world.y[i]));
    let free = keep(world, r => r >= B && r < 2 * B), inside = keep(world, r => r < B), anchors = [];
    for (let i = 0; i < free.x.length; i++) { const x = free.x[i], y = free.y[i]; if (outRad(x, y) === 2 * B - 1) for (const [dx, dy] of NB) if (outRad(x + dx, y + dy) === 2 * B && all.has(key(x + dx, y + dy))) anchors.push(x, y, dx, dy); }
    const edgeLab = new Map(); if (world.lab) for (let i = 0; i < world.x.length; i++) if (world.lab[i]) edgeLab.set(key(world.x[i], world.y[i]), world.lab[i]);
    info.links = 0; info.unjoined = 0;
    for (let j = 1; j <= J; j++) {
      const b = B << j, pinned = prev[j - 1];
      if (j > level) { free = twice(free); st.push(join(pinned, free)); continue; }
      // where the ring touched the inside at the last stage, it must touch it at this one
      const was = new Set(); for (let i = 0; i < inside.x.length; i++) was.add(key(inside.x[i], inside.y[i]));
      const touch = [];
      for (let i = 0; i < free.x.length; i++) {
        const x = free.x[i], y = free.y[i];
        for (const [dx, dy] of NB) { const k = key(x + dx, y + dy); if (was.has(k)) touch.push(x, y, dx, dy, j === 1 && known.has(edgeLab.get(k)) ? edgeLab.get(k) : 0); }
      }
      const held = new Set(), next = [];
      for (let k = 0; k < anchors.length; k += 4) for (let t = P.ANCHOR === 1 ? (rnd() < 0.5 ? 0 : 1) : 0, t1 = P.ANCHOR === 1 ? t + 1 : P.ANCHOR === 0 ? 0 : 2; t < t1; t++) {
        const dx = anchors[k + 2], dy = anchors[k + 3], x = dx ? 2 * anchors[k] + (dx > 0 ? 1 : 0) : 2 * anchors[k] + t, y = dy ? 2 * anchors[k + 1] + (dy > 0 ? 1 : 0) : 2 * anchors[k + 1] + t;
        next.push(x, y, dx, dy); held.add(key(x, y));
      }
      const f2 = twice(free), fx = [], fy = [], rx = [], ry = [];
      for (let i = 0; i < f2.x.length; i++) if (held.has(key(f2.x[i], f2.y[i]))) { fx.push(f2.x[i]); fy.push(f2.y[i]); } else { rx.push(f2.x[i]); ry.push(f2.y[i]); }
      const np = pinned.x.length;
      E.clear(); put(0, pinned.x, pinned.y); put(np, fx, fy); put(np + fx.length, rx, ry); E.place(np + fx.length + rx.length);
      if (j === 1 && world.lab) new Int32Array(E.memory.buffer, E.partlabels(), pinLab.length).set(pinLab);
      for (let k = 0; k < touch.length; k += 5) { const a = E.link(touch[k], touch[k + 1], touch[k + 2], touch[k + 3], b, np, P.REACH, touch[k + 4]); if (a < 0) info.unjoined++; else info.links += a; }
      E.pin(np + fx.length, b, 2 * b); E.relax(4 * P.T);
      free = get(np, E.count()); anchors = next; inside = pinned;
      st.push(join(pinned, free));
    }
    prev = st; info.ms = Date.now() - t0;
    const pic = picture(st[J]); pic.info = info; return pic;
  }
  return { begin, next, halve, G, get level() { return level; }, get world() { return world; }, get stages() { return prev; } };
}

if (typeof module !== 'undefined') module.exports = { OUT, outRad, Outgoer };
