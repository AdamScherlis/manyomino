//! Fast fixed-n single-cell MCMC sampler for uniform random polyominoes.
//!
//! No external crates (network-restricted build environment): own RNG
//! (xoshiro256++), grid-backed indexed sets, ring-test + BFS connectivity.
//!
//! The chain (see handoff doc): pick cell c uniformly; if A-{c} disconnects,
//! stay put (still a step); else pick s uniformly from the perimeter *sites*
//! of A-{c} (c itself included, so s=c is legal) and move c -> s.  The kernel
//! is symmetric, so uniform over fixed-n polyominoes is stationary.

use std::collections::HashMap;
use std::io::Write;

// ---------------------------------------------------------------- RNG

/// xoshiro256++, seeded via splitmix64.
struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        let mut x = seed;
        let mut next = || {
            x = x.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = x;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }

    #[inline(always)]
    fn next_u64(&mut self) -> u64 {
        let r = self.s[0]
            .wrapping_add(self.s[3])
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        r
    }

    /// Uniform in [0, n), unbiased (Lemire).
    #[inline(always)]
    fn below(&mut self, n: u64) -> u64 {
        let mut x = self.next_u64();
        let mut m = (x as u128) * (n as u128);
        let mut l = m as u64;
        if l < n {
            let t = n.wrapping_neg() % n;
            while l < t {
                x = self.next_u64();
                m = (x as u128) * (n as u128);
                l = m as u64;
            }
        }
        (m >> 64) as u64
    }
}

// ---------------------------------------------------------------- state

/// Polyomino on a bounded grid with a safety margin; rebuilt (recentred,
/// resized) whenever a cell lands too close to the border, so neighbor
/// arithmetic never wraps.  Absolute position is meaningless (the chain's
/// centroid random-walks); only shape matters.
struct State {
    w: usize,
    h: usize,
    /// occupancy bitmap, bit per grid position (cache-friendly hot path)
    occ_bits: Vec<u64>,
    /// index+1 into `cells`, 0 = empty (cold path: indexed-set bookkeeping)
    cell_idx: Vec<u32>,
    /// index+1 into `per`, 0 = not a perimeter site.
    per_idx: Vec<u32>,
    cells: Vec<u32>,
    per: Vec<u32>,
    visited: Vec<u32>,
    stamp: u32,
    /// 4 neighbor cell-indices per cell (NONE = empty side), kept in sync
    /// with the cell set; lets the BFS run over a compact index space
    nbr: Vec<[u32; 4]>,
    /// per-front visited bitplanes over cell *indices* (~n/8 bytes each, so
    /// they live in L1); bits are undone via the queues after every call
    fb: Vec<Vec<u64>>,
    queue: Vec<u32>,
    fqueues: Vec<Vec<u32>>,
    // coordinate sums for O(1) radius of gyration
    sx: i64,
    sy: i64,
    sx2: i64,
    sy2: i64,
    // move statistics
    n_cut_reject: u64,
    n_noop: u64,
    n_moved: u64,
    n_ring_pass: u64,
    n_leaf_pass: u64,
    n_bfs: u64,
    n_bfs_nodes: u64,
    n_rebuild: u64,
}

impl State {
    fn new(n: usize, init: &str) -> Self {
        assert!(n >= 2);
        let pts: Vec<(usize, usize)> = match init {
            "bar" => (0..n).map(|i| (i, 0)).collect(),
            "rect" => {
                let w0 = ((n as f64).sqrt().round() as usize).max(1);
                (0..n).map(|i| (i % w0, i / w0)).collect()
            }
            _ => panic!("unknown init {init}"),
        };
        let bw = pts.iter().map(|p| p.0).max().unwrap() + 1;
        let bh = pts.iter().map(|p| p.1).max().unwrap() + 1;
        let mx = (bw / 4).max(16);
        let my = (bh / 4).max(16);
        let (w, h) = (bw + 2 * mx, bh + 2 * my);
        let mut st = State {
            w,
            h,
            occ_bits: vec![0; (w * h + 63) / 64],
            cell_idx: vec![0; w * h],
            per_idx: vec![0; w * h],
            cells: Vec::with_capacity(n),
            per: Vec::new(),
            visited: vec![0; w * h],
            stamp: 0,
            nbr: Vec::with_capacity(n),
            fb: vec![vec![0; (n + 63) / 64]; 4],
            queue: Vec::new(),
            fqueues: vec![Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            sx: 0,
            sy: 0,
            sx2: 0,
            sy2: 0,
            n_cut_reject: 0,
            n_noop: 0,
            n_moved: 0,
            n_ring_pass: 0,
            n_leaf_pass: 0,
            n_bfs: 0,
            n_bfs_nodes: 0,
            n_rebuild: 0,
        };
        for (x, y) in pts {
            st.add_cell(((y + my) * w + x + mx) as u32);
        }
        st.rebuild_perimeter();
        st
    }

    #[inline(always)]
    fn occ(&self, pos: u32) -> bool {
        (self.occ_bits[(pos >> 6) as usize] >> (pos & 63)) & 1 != 0
    }

    const NONE: u32 = u32::MAX;

    #[inline(always)]
    fn dirs(&self) -> [i64; 4] {
        let w = self.w as i64;
        [1, -1, w, -w] // opposite of slot d is slot d ^ 1
    }

    #[inline(always)]
    fn add_cell(&mut self, pos: u32) {
        let i = self.cells.len() as u32;
        self.cells.push(pos);
        self.cell_idx[pos as usize] = i + 1;
        let mut links = [Self::NONE; 4];
        for (d, dd) in self.dirs().into_iter().enumerate() {
            let q = (pos as i64 + dd) as u32;
            if self.occ(q) {
                let j = self.cell_idx[q as usize] - 1;
                links[d] = j;
                self.nbr[j as usize][d ^ 1] = i;
            }
        }
        self.nbr.push(links);
        self.occ_bits[(pos >> 6) as usize] |= 1u64 << (pos & 63);
        let (x, y) = ((pos as usize % self.w) as i64, (pos as usize / self.w) as i64);
        self.sx += x;
        self.sy += y;
        self.sx2 += x * x;
        self.sy2 += y * y;
    }

    /// Swap-remove the cell at index `ci` (position `pos`), keeping the
    /// adjacency links consistent.
    #[inline(always)]
    fn remove_cell(&mut self, ci: usize, pos: u32) {
        // unlink ci from its neighbors
        for d in 0..4 {
            let j = self.nbr[ci][d];
            if j != Self::NONE {
                self.nbr[j as usize][d ^ 1] = Self::NONE;
            }
        }
        let last_idx = self.cells.len() - 1;
        if ci != last_idx {
            let last_pos = self.cells[last_idx];
            self.cells[ci] = last_pos;
            self.cell_idx[last_pos as usize] = ci as u32 + 1;
            self.nbr[ci] = self.nbr[last_idx];
            for d in 0..4 {
                let j = self.nbr[ci][d];
                if j != Self::NONE {
                    self.nbr[j as usize][d ^ 1] = ci as u32;
                }
            }
        }
        self.cells.pop();
        self.nbr.pop();
        self.cell_idx[pos as usize] = 0;
        self.occ_bits[(pos >> 6) as usize] &= !(1u64 << (pos & 63));
        let (x, y) = ((pos as usize % self.w) as i64, (pos as usize / self.w) as i64);
        self.sx -= x;
        self.sy -= y;
        self.sx2 -= x * x;
        self.sy2 -= y * y;
    }

    #[inline(always)]
    fn per_add(&mut self, pos: u32) {
        self.per.push(pos);
        self.per_idx[pos as usize] = self.per.len() as u32;
    }

    #[inline(always)]
    fn per_remove(&mut self, pos: u32) {
        let i = self.per_idx[pos as usize] as usize - 1;
        let last = *self.per.last().unwrap();
        self.per[i] = last;
        self.per_idx[last as usize] = i as u32 + 1;
        self.per.pop();
        self.per_idx[pos as usize] = 0;
    }

    fn rebuild_perimeter(&mut self) {
        for &p in &self.per {
            self.per_idx[p as usize] = 0;
        }
        self.per.clear();
        let w = self.w as i64;
        for i in 0..self.cells.len() {
            let c = self.cells[i] as i64;
            for d in [1, -1, w, -w] {
                let q = (c + d) as u32;
                if !self.occ(q) && self.per_idx[q as usize] == 0 {
                    self.per_add(q);
                }
            }
        }
    }

    #[inline(always)]
    fn has_occupied_neighbor(&self, pos: u32) -> bool {
        let w = self.w as i64;
        let p = pos as i64;
        self.occ((p + 1) as u32)
            || self.occ((p - 1) as u32)
            || self.occ((p + w) as u32)
            || self.occ((p - w) as u32)
    }

    /// True iff A - {cells[ci]} stays connected.  Cheap tests first (single
    /// occupied neighbor; 3x3 ring arc test), then BFS with early exit.
    #[inline(always)]
    fn removable(&mut self, ci: usize) -> bool {
        let c = self.cells[ci];
        let w = self.w as i64;
        let p = c as i64;
        // ring in cyclic order: E, NE, N, NW, W, SW, S, SE
        let offs = [1, 1 - w, -w, -w - 1, -1, w - 1, w, w + 1];
        let mut b = [false; 8];
        let mut orth = 0;
        for i in 0..8 {
            b[i] = self.occ((p + offs[i]) as u32);
            if i % 2 == 0 && b[i] {
                orth += 1;
            }
        }
        if orth == 1 {
            self.n_leaf_pass += 1;
            return true;
        }
        // count 0->1 transitions around the ring; a single occupied arc
        // contains all occupied orthogonal neighbors and is itself
        // 4-connected (consecutive ring cells are 4-adjacent), so removal
        // is safe.  Sufficient, not necessary.
        let mut trans = 0;
        for i in 0..8 {
            if !b[i] && b[(i + 1) % 8] {
                trans += 1;
            }
        }
        if trans <= 1 {
            self.n_ring_pass += 1;
            return true;
        }
        // Partition the occupied ring cells into maximal cyclic arcs.
        // Consecutive ring cells are 4-adjacent (orthogonal/diagonal
        // alternation), so each arc is 4-connected within A - {c}.  Only
        // arcs containing an orthogonal neighbor of c matter: c's removal
        // is safe iff those arcs end up in one component of A - {c}.
        let start = (0..8).find(|&i| !b[i]).unwrap(); // trans >= 2 => an empty exists
        let mut arcs = [[0usize; 8]; 4];
        let mut arc_len = [0usize; 4];
        let mut n_arcs = 0usize;
        let mut in_arc = false;
        for k in 1..=8 {
            let i = (start + k) % 8;
            if b[i] {
                if !in_arc {
                    in_arc = true;
                    n_arcs += 1;
                }
                let a = n_arcs - 1;
                arcs[a][arc_len[a]] = i;
                arc_len[a] += 1;
            } else {
                in_arc = false;
            }
        }
        // fronts = arcs holding at least one orthogonal neighbor
        let mut fronts = [0usize; 4];
        let mut nf = 0usize;
        for a in 0..n_arcs {
            if (0..arc_len[a]).any(|j| arcs[a][j] % 2 == 0) {
                fronts[nf] = a;
                nf += 1;
            }
        }
        if nf == 1 {
            // all orthogonal neighbors sit in one arc (ring test missed it
            // because of an unrelated diagonal arc)
            self.n_ring_pass += 1;
            return true;
        }
        // Seed each front with the *orthogonal* members of its arc, as cell
        // indices from the adjacency links (diagonal arc members are reached
        // by the BFS itself; ring index -> direction slot: E->0 N->3 W->1
        // S->2).  The BFS then runs entirely in compact index space.
        const RING_TO_SLOT: [usize; 8] = [0, 9, 3, 9, 1, 9, 2, 9];
        let mut seeds = [[0u32; 4]; 4];
        let mut seed_cnt = [0usize; 4];
        for f in 0..nf {
            let a = fronts[f];
            for j in 0..arc_len[a] {
                let ri = arcs[a][j];
                if ri % 2 == 0 {
                    seeds[f][seed_cnt[f]] = self.nbr[ci][RING_TO_SLOT[ri]];
                    seed_cnt[f] += 1;
                }
            }
        }
        self.n_bfs += 1;
        self.multifront_removable(ci as u32, &seeds, &seed_cnt, nf)
    }

    /// Alternating multi-source BFS in A - {c}: one front per neighbor arc,
    /// advanced round-robin.  Fronts that touch merge (union-find).  Stops
    /// as soon as either all fronts merge (removable) or some merged group
    /// exhausts its queues without containing every front (cut cell) — so a
    /// cut test costs O(smallest component), not O(n).
    fn multifront_removable(
        &mut self,
        ci: u32,
        seeds: &[[u32; 4]; 4],
        seed_cnt: &[usize; 4],
        nf: usize,
    ) -> bool {
        let mut queues = std::mem::take(&mut self.fqueues);
        let mut fb = std::mem::take(&mut self.fb);
        let mut heads = [0usize; 4];
        for f in 0..nf {
            queues[f].clear();
            for j in 0..seed_cnt[f] {
                let v = seeds[f][j];
                fb[f][(v >> 6) as usize] |= 1u64 << (v & 63);
                queues[f].push(v);
            }
        }

        // tiny union-find over <= 4 fronts
        let mut parent = [0u8, 1, 2, 3];
        fn find(parent: &mut [u8; 4], mut i: u8) -> u8 {
            while parent[i as usize] != i {
                parent[i as usize] = parent[parent[i as usize] as usize];
                i = parent[i as usize];
            }
            i
        }
        let mut roots = nf as u32;

        let result = 'outer: loop {
            let mut any = false;
            for f in 0..nf {
                if heads[f] < queues[f].len() {
                    any = true;
                    let u = queues[f][heads[f]] as usize;
                    heads[f] += 1;
                    self.n_bfs_nodes += 1;
                    let links = self.nbr[u];
                    for d in 0..4 {
                        let v = links[d];
                        if v == Self::NONE || v == ci {
                            continue;
                        }
                        let (wi, bi) = ((v >> 6) as usize, v & 63);
                        if fb[f][wi] >> bi & 1 != 0 {
                            continue; // already visited by own front
                        }
                        let mut owner = usize::MAX;
                        for g in 0..nf {
                            if g != f && fb[g][wi] >> bi & 1 != 0 {
                                owner = g;
                                break;
                            }
                        }
                        if owner != usize::MAX {
                            let (rf, rg) =
                                (find(&mut parent, f as u8), find(&mut parent, owner as u8));
                            if rf != rg {
                                parent[rf as usize] = rg;
                                roots -= 1;
                                if roots == 1 {
                                    break 'outer true;
                                }
                            }
                        } else {
                            fb[f][wi] |= 1u64 << bi;
                            queues[f].push(v);
                        }
                    }
                } else {
                    // front f's queue is spent; if every front in its merged
                    // group is spent, that component is fully explored
                    let rf = find(&mut parent, f as u8);
                    let group_done = (0..nf)
                        .filter(|&g| find(&mut parent, g as u8) == rf)
                        .all(|g| heads[g] >= queues[g].len());
                    if group_done {
                        break 'outer roots == 1;
                    }
                }
            }
            if !any {
                break roots == 1;
            }
        };
        // undo the visited bits (every marked cell is in exactly one queue)
        for f in 0..nf {
            for &v in queues[f].iter() {
                fb[f][(v >> 6) as usize] &= !(1u64 << (v & 63));
            }
            queues[f].clear();
        }
        self.fqueues = queues;
        self.fb = fb;
        result
    }

    /// Obviously-correct connectivity reference: full BFS over A - {c}.
    /// Used by selftest to differential-test `removable`.
    fn removable_ref(&mut self, c: u32) -> bool {
        let n = self.cells.len();
        if self.stamp >= u32::MAX - 1 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.stamp = 0;
        }
        self.stamp += 1;
        let stamp = self.stamp;
        let w = self.w as i64;
        let start = if self.cells[0] == c { self.cells[1] } else { self.cells[0] };
        self.queue.clear();
        self.queue.push(start);
        self.visited[start as usize] = stamp;
        let mut head = 0;
        let mut cnt = 1;
        while head < self.queue.len() {
            let u = self.queue[head] as i64;
            head += 1;
            for d in [1, -1, w, -w] {
                let v = (u + d) as u32;
                if v != c && self.occ(v) && self.visited[v as usize] != stamp {
                    self.visited[v as usize] = stamp;
                    cnt += 1;
                    self.queue.push(v);
                }
            }
        }
        cnt == n - 1
    }

    /// One step of the chain.  Returns true iff the state changed.
    #[inline(always)]
    fn step(&mut self, rng: &mut Rng) -> bool {
        let n = self.cells.len();
        let ci = rng.below(n as u64) as usize;
        let c = self.cells[ci];
        if !self.removable(ci) {
            self.n_cut_reject += 1;
            return false;
        }
        let w = self.w as i64;
        // remove c; perimeter becomes that of A - {c}: c itself is now a
        // perimeter site, and former sites that only touched c drop out.
        self.remove_cell(ci, c);
        self.per_add(c);
        for d in [1, -1, w, -w] {
            let q = (c as i64 + d) as u32;
            if !self.occ(q) && self.per_idx[q as usize] != 0 && !self.has_occupied_neighbor(q) {
                self.per_remove(q);
            }
        }
        // sample the new site uniformly over perimeter *sites* (s = c legal)
        let s = self.per[rng.below(self.per.len() as u64) as usize];
        self.per_remove(s);
        self.add_cell(s);
        for d in [1, -1, w, -w] {
            let q = (s as i64 + d) as u32;
            if !self.occ(q) && self.per_idx[q as usize] == 0 {
                self.per_add(q);
            }
        }
        let (x, y) = (s as usize % self.w, s as usize / self.w);
        if x < 2 || x >= self.w - 2 || y < 2 || y >= self.h - 2 {
            self.rebuild();
        }
        if s != c {
            self.n_moved += 1;
            true
        } else {
            self.n_noop += 1;
            false
        }
    }

    /// Recentre and resize the grid; recompute perimeter and sums.
    fn rebuild(&mut self) {
        self.n_rebuild += 1;
        let w = self.w;
        let (mut minx, mut miny, mut maxx, mut maxy) = (usize::MAX, usize::MAX, 0, 0);
        for &p in &self.cells {
            let (x, y) = (p as usize % w, p as usize / w);
            minx = minx.min(x);
            miny = miny.min(y);
            maxx = maxx.max(x);
            maxy = maxy.max(y);
        }
        let (bw, bh) = (maxx - minx + 1, maxy - miny + 1);
        let mx = (bw / 4).max(16);
        let my = (bh / 4).max(16);
        let (nw, nh) = (bw + 2 * mx, bh + 2 * my);
        let old: Vec<u32> = std::mem::take(&mut self.cells);
        self.w = nw;
        self.h = nh;
        self.occ_bits = vec![0; (nw * nh + 63) / 64];
        self.cell_idx = vec![0; nw * nh];
        self.per_idx = vec![0; nw * nh];
        self.visited = vec![0; nw * nh];
        self.nbr.clear();
        self.stamp = 0;
        self.per.clear();
        self.sx = 0;
        self.sy = 0;
        self.sx2 = 0;
        self.sy2 = 0;
        for &p in &old {
            let (x, y) = (p as usize % w, p as usize / w);
            self.add_cell(((y - miny + my) * nw + (x - minx + mx)) as u32);
        }
        self.rebuild_perimeter();
    }

    fn rg2(&self) -> f64 {
        let n = self.cells.len() as f64;
        let mx = self.sx as f64 / n;
        let my = self.sy as f64 / n;
        (self.sx2 + self.sy2) as f64 / n - mx * mx - my * my
    }

    fn coords(&self) -> Vec<(i64, i64)> {
        self.cells
            .iter()
            .map(|&p| ((p as usize % self.w) as i64, (p as usize / self.w) as i64))
            .collect()
    }

    /// Translation-canonical packed form for n <= 8 (4 bits per coordinate).
    fn canonical_u64(&self) -> u64 {
        let cs = self.coords();
        let minx = cs.iter().map(|c| c.0).min().unwrap();
        let miny = cs.iter().map(|c| c.1).min().unwrap();
        let mut bytes: Vec<u8> = cs
            .iter()
            .map(|c| (((c.0 - minx) as u8) << 4) | ((c.1 - miny) as u8))
            .collect();
        bytes.sort_unstable();
        bytes.iter().fold(0u64, |acc, &b| (acc << 8) | b as u64)
    }

    /// Full from-scratch consistency check (differential test vs the
    /// incremental structures).  Panics on any mismatch.
    fn check_invariants(&mut self, n: usize) {
        assert_eq!(self.cells.len(), n, "cell count");
        // no duplicate cells; index map consistent
        for (i, &p) in self.cells.iter().enumerate() {
            assert_eq!(self.cell_idx[p as usize], i as u32 + 1, "cell_idx");
        }
        // connectivity: BFS over all cells
        if self.stamp >= u32::MAX - 1 {
            self.visited.iter_mut().for_each(|v| *v = 0);
            self.stamp = 0;
        }
        self.stamp += 1;
        let stamp = self.stamp;
        let w = self.w as i64;
        self.queue.clear();
        self.queue.push(self.cells[0]);
        self.visited[self.cells[0] as usize] = stamp;
        let mut head = 0;
        let mut cnt = 1;
        while head < self.queue.len() {
            let u = self.queue[head] as i64;
            head += 1;
            for d in [1, -1, w, -w] {
                let v = (u + d) as u32;
                if self.occ(v) && self.visited[v as usize] != stamp {
                    self.visited[v as usize] = stamp;
                    cnt += 1;
                    self.queue.push(v);
                }
            }
        }
        assert_eq!(cnt, n, "connectivity");
        // perimeter set == from-scratch recompute
        let mut fresh: Vec<u32> = Vec::new();
        self.stamp += 1;
        let stamp = self.stamp;
        for i in 0..self.cells.len() {
            let c = self.cells[i] as i64;
            for d in [1, -1, w, -w] {
                let q = (c + d) as u32;
                if !self.occ(q) && self.visited[q as usize] != stamp {
                    self.visited[q as usize] = stamp;
                    fresh.push(q);
                }
            }
        }
        assert_eq!(fresh.len(), self.per.len(), "perimeter size");
        for &q in &fresh {
            assert!(self.per_idx[q as usize] != 0, "perimeter membership");
        }
        // adjacency links == from-scratch recompute
        assert_eq!(self.nbr.len(), self.cells.len(), "nbr length");
        for i in 0..self.cells.len() {
            let pos = self.cells[i] as i64;
            for (d, dd) in self.dirs().into_iter().enumerate() {
                let q = (pos + dd) as u32;
                let want = if self.occ(q) {
                    self.cell_idx[q as usize] - 1
                } else {
                    Self::NONE
                };
                assert_eq!(self.nbr[i][d], want, "nbr link");
            }
        }
        // coordinate sums
        let (mut sx, mut sy, mut sx2, mut sy2) = (0i64, 0i64, 0i64, 0i64);
        for &p in &self.cells {
            let (x, y) = ((p as usize % self.w) as i64, (p as usize / self.w) as i64);
            sx += x;
            sy += y;
            sx2 += x * x;
            sy2 += y * y;
        }
        assert_eq!((sx, sy, sx2, sy2), (self.sx, self.sy, self.sx2, self.sy2), "sums");
    }

    fn stats_line(&self) -> String {
        let total = self.n_cut_reject + self.n_noop + self.n_moved;
        format!(
            "steps={} moved={:.4} noop={:.4} cut_reject={:.4} leaf={:.4} ring={:.4} bfs={:.4} bfs_nodes_per_step={:.2} rebuilds={}",
            total,
            self.n_moved as f64 / total as f64,
            self.n_noop as f64 / total as f64,
            self.n_cut_reject as f64 / total as f64,
            self.n_leaf_pass as f64 / total as f64,
            self.n_ring_pass as f64 / total as f64,
            self.n_bfs as f64 / total as f64,
            self.n_bfs_nodes as f64 / total as f64,
            self.n_rebuild,
        )
    }
}

// ---------------------------------------------------------------- CLI

fn parse_args() -> HashMap<String, String> {
    let mut m = HashMap::new();
    let args: Vec<String> = std::env::args().collect();
    m.insert("mode".to_string(), args.get(1).cloned().unwrap_or_default());
    let mut i = 2;
    while i + 1 < args.len() + 1 {
        if i + 1 >= args.len() {
            break;
        }
        let k = args[i].trim_start_matches("--").to_string();
        m.insert(k, args[i + 1].clone());
        i += 2;
    }
    m
}

fn get<T: std::str::FromStr>(m: &HashMap<String, String>, k: &str, default: T) -> T {
    m.get(k).and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn main() {
    let args = parse_args();
    let mode = args.get("mode").cloned().unwrap_or_default();
    let n: usize = get(&args, "n", 100);
    let seed: u64 = get(&args, "seed", 1);
    let init = args.get("init").cloned().unwrap_or_else(|| "bar".into());
    let mut rng = Rng::new(seed);

    match mode.as_str() {
        "bench" => {
            let steps: u64 = get(&args, "steps", 10_000_000);
            let mut st = State::new(n, &init);
            let t0 = std::time::Instant::now();
            for _ in 0..steps {
                st.step(&mut rng);
            }
            let dt = t0.elapsed().as_secs_f64();
            println!(
                "n={} steps={} secs={:.3} moves_per_sec={:.3e}",
                n, steps, dt, steps as f64 / dt
            );
            eprintln!("{}", st.stats_line());
        }

        "bench2" => {
            let steps: u64 = get(&args, "steps", 10_000_000);
            let mut st = State::new(n, &init);
            // equilibrate a bit first so stats reflect the same shape mix
            for _ in 0..steps / 2 {
                st.step(&mut rng);
            }
            let t0 = std::time::Instant::now();
            let mut acc = 0u64;
            for _ in 0..steps {
                let ci = rng.below(st.cells.len() as u64) as usize;
                if st.removable(ci) {
                    acc += 1;
                }
            }
            let dt = t0.elapsed().as_secs_f64();
            println!(
                "removable-only n={} steps={} secs={:.3} tests_per_sec={:.3e} frac_removable={:.4}",
                n, steps, dt, steps as f64 / dt, acc as f64 / steps as f64
            );
            eprintln!("{}", st.stats_line());
        }
        "chi" => {
            // print "canonical_u64 count" lines for chi-square in Python
            assert!(n <= 8, "canonical_u64 packing supports n <= 8");
            let obs: u64 = get(&args, "obs", 1_000_000);
            let thin: u64 = get(&args, "thin", 10);
            let burn: u64 = get(&args, "burn", 100_000);
            let check_every: u64 = get(&args, "check-every", 0);
            let mut st = State::new(n, &init);
            for _ in 0..burn {
                st.step(&mut rng);
            }
            let mut counts: HashMap<u64, u64> = HashMap::new();
            for i in 0..obs {
                for _ in 0..thin {
                    st.step(&mut rng);
                }
                *counts.entry(st.canonical_u64()).or_insert(0) += 1;
                if check_every > 0 && i % check_every == 0 {
                    st.check_invariants(n);
                }
            }
            let out = std::io::stdout();
            let mut w = std::io::BufWriter::new(out.lock());
            for (k, v) in &counts {
                writeln!(w, "{k} {v}").unwrap();
            }
            w.flush().unwrap();
            eprintln!("{}", st.stats_line());
        }
        "run" => {
            // time series of rg2 (CSV step,rg2) + optional final-shape dump
            let steps: u64 = get(&args, "steps", 1_000_000);
            let record_every: u64 = get(&args, "record-every", 10_000);
            let check_every: u64 = get(&args, "check-every", 0);
            let out_path = args.get("out").cloned();
            let dump_path = args.get("dump-final").cloned();
            let dump_every: u64 = get(&args, "dump-every", 0);
            let dump_prefix = args.get("dump-prefix").cloned();
            let mut st = State::new(n, &init);
            let mut writer = out_path.map(|p| {
                let f = std::fs::File::create(p).expect("create out");
                let mut w = std::io::BufWriter::new(f);
                writeln!(w, "step,rg2").unwrap();
                w
            });
            let t0 = std::time::Instant::now();
            let mut done: u64 = 0;
            while done < steps {
                let chunk = record_every.min(steps - done);
                for _ in 0..chunk {
                    st.step(&mut rng);
                }
                done += chunk;
                if let Some(w) = writer.as_mut() {
                    writeln!(w, "{},{:.6}", done, st.rg2()).unwrap();
                }
                if check_every > 0 && done % check_every == 0 {
                    st.check_invariants(n);
                }
                if dump_every > 0 && done % dump_every == 0 {
                    if let Some(prefix) = dump_prefix.as_ref() {
                        let f = std::fs::File::create(format!("{prefix}_{done}.txt"))
                            .expect("create dump");
                        let mut w = std::io::BufWriter::new(f);
                        for (x, y) in st.coords() {
                            writeln!(w, "{x} {y}").unwrap();
                        }
                        w.flush().unwrap();
                    }
                }
            }
            let dt = t0.elapsed().as_secs_f64();
            if let Some(mut w) = writer {
                w.flush().unwrap();
            }
            if let Some(p) = dump_path {
                let f = std::fs::File::create(p).expect("create dump");
                let mut w = std::io::BufWriter::new(f);
                for (x, y) in st.coords() {
                    writeln!(w, "{x} {y}").unwrap();
                }
                w.flush().unwrap();
            }
            eprintln!(
                "n={} init={} seed={} secs={:.1} moves_per_sec={:.3e}",
                n,
                init,
                seed,
                dt,
                steps as f64 / dt
            );
            eprintln!("{}", st.stats_line());
        }
        "selftest" => {
            // heavy invariant checking at several sizes
            for (nn, steps, every) in [(2usize, 200_000u64, 100u64), (3, 200_000, 100),
                                       (5, 200_000, 100), (10, 200_000, 200),
                                       (37, 100_000, 500), (100, 100_000, 1000),
                                       (1000, 50_000, 5000)] {
                for init in ["bar", "rect"] {
                    let mut st = State::new(nn, init);
                    let mut r = Rng::new(seed ^ nn as u64);
                    for i in 0..steps {
                        st.step(&mut r);
                        if i % every == 0 {
                            st.check_invariants(nn);
                            // differential-test the fast connectivity check
                            // against a full BFS on a batch of cells
                            let m = nn.min(64);
                            for _ in 0..m {
                                let idx = r.below(nn as u64) as usize;
                                let cpos = st.cells[idx];
                                assert_eq!(
                                    st.removable(idx),
                                    st.removable_ref(cpos),
                                    "removable mismatch at n={nn}"
                                );
                            }
                        }
                    }
                    st.check_invariants(nn);
                    println!("selftest n={nn} init={init} OK  {}", st.stats_line());
                }
            }
        }
        _ => {
            eprintln!("modes: bench | chi | run | selftest");
            eprintln!("  common: --n N --seed S --init bar|rect");
            eprintln!("  bench: --steps T");
            eprintln!("  chi:   --obs M --thin K --burn B [--check-every X]");
            eprintln!("  run:   --steps T --record-every R [--out ts.csv] [--dump-final cells.txt] [--check-every X]");
            std::process::exit(2);
        }
    }
}
