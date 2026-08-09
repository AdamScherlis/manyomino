//! Fast fixed-n single-cell MCMC sampler for uniform random polyominoes.
//!
//! No external crates (network-restricted build environment): own RNG
//! (xoshiro256++), grid-backed indexed sets, ring-test + BFS connectivity.
//!
//! The chain (see handoff doc): pick cell c uniformly; if A-{c} disconnects,
//! stay put (still a step); else pick s uniformly from the perimeter *sites*
//! of A-{c} (c itself included, so s=c is legal) and move c -> s.  The kernel
//! is symmetric, so uniform over fixed-n polyominoes is stationary.

mod perm;

use std::collections::HashMap;
use std::io::Write;

// ---------------------------------------------------------------- RNG

/// xoshiro256++, seeded via splitmix64.
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    pub fn new(seed: u64) -> Self {
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
    pub fn next_u64(&mut self) -> u64 {
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
    n_cp_nopair: u64,
    n_cp_notbridge: u64,
    n_cp_toobig: u64,
    n_cp_invalid: u64,
    n_cp_accept: u64,
    /// cumulative coordinate shift from grid rebuilds (for external watchers)
    shift_x: i64,
    shift_y: i64,
    /// number of adjacent occupied pairs (edges); cycles = E - n + 1
    e_edges: u64,
    /// bridge-test node cap per front (0 = uncapped); give-up is symmetric:
    /// on any acceptable move the forward and reverse tests resolve at the
    /// same per-front pop counts (min(|B|,|R|)), so capping cannot break
    /// detailed balance; capped give-ups are rejections (self-loops)
    cp_cap: u64,
    n_cp_capped: u64,
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
            // 2 x (n/2) bar: geometrically as stringy as a 1 x n bar, but
            // with no interior cut cells, so early moves stay cheap at
            // large n (a 1 x n bar costs O(n) per cut test while unfolding)
            "bar2" => {
                let w0 = n.div_ceil(2);
                (0..n).map(|i| (i % w0, i / w0)).collect()
            }
            // width-8 bar: Rg well above equilibrium at large n, and thick
            // enough that random nicks never sever it, so cut tests stay
            // local while it crumples (bar2 develops O(n)-cost cut cells
            // as soon as opposing nicks appear)
            "bar8" => {
                let rows = 8.min(n);
                let w0 = n.div_ceil(rows);
                (0..n).map(|i| (i % w0, i / w0)).collect()
            }
            _ => panic!("unknown init {init}"),
        };
        Self::from_points(pts)
    }

    /// Resume from a cell dump ("x y" per line).
    fn from_file(path: &str) -> Self {
        let text = std::fs::read_to_string(path).expect("read init file");
        let mut pts: Vec<(i64, i64)> = Vec::new();
        for line in text.lines() {
            let mut it = line.split_whitespace();
            let x: i64 = it.next().unwrap().parse().unwrap();
            let y: i64 = it.next().unwrap().parse().unwrap();
            pts.push((x, y));
        }
        let minx = pts.iter().map(|p| p.0).min().unwrap();
        let miny = pts.iter().map(|p| p.1).min().unwrap();
        Self::from_points(
            pts.iter()
                .map(|p| ((p.0 - minx) as usize, (p.1 - miny) as usize))
                .collect(),
        )
    }

    fn from_points(pts: Vec<(usize, usize)>) -> Self {
        let n = pts.len();
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
            n_cp_nopair: 0,
            n_cp_notbridge: 0,
            n_cp_toobig: 0,
            n_cp_invalid: 0,
            n_cp_accept: 0,
            shift_x: 0,
            shift_y: 0,
            cp_cap: 0,
            n_cp_capped: 0,
            e_edges: 0,
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
        self.e_edges += links.iter().filter(|&&l| l != Self::NONE).count() as u64;
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
        self.e_edges -= self.nbr[ci].iter().filter(|&&l| l != Self::NONE).count() as u64;
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

    /// If the adjacency edge (ui, vi) is a bridge of the cell graph, return
    /// the cell *indices* of the component containing vi (must be <= n/2
    /// cells, else None).  Alternating two-front BFS excluding that edge, so
    /// a non-bridge (cycle edge) is detected as soon as the fronts meet.
    fn bridge_component(&mut self, ui: usize, vi: usize) -> Option<Vec<u32>> {
        let n = self.cells.len();
        let mut queues = std::mem::take(&mut self.fqueues);
        let mut fb = std::mem::take(&mut self.fb);
        queues[0].clear();
        queues[1].clear();
        fb[0][ui >> 6] |= 1u64 << (ui & 63);
        fb[1][vi >> 6] |= 1u64 << (vi & 63);
        queues[0].push(ui as u32);
        queues[1].push(vi as u32);
        let mut heads = [0usize; 2];
        let cap = if self.cp_cap == 0 { u64::MAX } else { self.cp_cap };
        // verdict: None undecided; Some(true) bridge; Some(false) not
        let mut verdict: Option<bool> = None;
        'outer: while verdict.is_none() {
            if heads[0] as u64 > cap && heads[1] as u64 > cap {
                // symmetric give-up: min component larger than the cap
                self.n_cp_capped += 1;
                break 'outer;
            }
            for f in 0..2 {
                if heads[f] >= queues[f].len() {
                    verdict = Some(true); // side f fully explored, never met
                    break 'outer;
                }
                let x = queues[f][heads[f]] as usize;
                heads[f] += 1;
                self.n_bfs_nodes += 1;
                let links = self.nbr[x];
                for d in 0..4 {
                    let y = links[d];
                    if y == Self::NONE {
                        continue;
                    }
                    let yu = y as usize;
                    if (x == ui && yu == vi) || (x == vi && yu == ui) {
                        continue; // the cut edge itself
                    }
                    let (wi, bi) = ((y >> 6) as usize, y & 63);
                    if fb[f][wi] >> bi & 1 != 0 {
                        continue;
                    }
                    if fb[1 - f][wi] >> bi & 1 != 0 {
                        verdict = Some(false); // fronts met: cycle edge
                        break 'outer;
                    }
                    fb[f][wi] |= 1u64 << bi;
                    queues[f].push(y);
                }
            }
        }
        let mut result = None;
        if verdict == Some(true) {
            // one side is complete; figure out sizes, enumerate the v side
            let u_done = heads[0] >= queues[0].len();
            if u_done {
                let b_size = n - queues[0].len();
                if 2 * b_size <= n {
                    // v side is exactly half; finish enumerating it
                    while heads[1] < queues[1].len() {
                        let x = queues[1][heads[1]] as usize;
                        heads[1] += 1;
                        let links = self.nbr[x];
                        for d in 0..4 {
                            let y = links[d];
                            if y == Self::NONE || (x == vi && y as usize == ui) {
                                continue;
                            }
                            let (wi, bi) = ((y >> 6) as usize, y & 63);
                            if fb[1][wi] >> bi & 1 == 0 {
                                fb[1][wi] |= 1u64 << bi;
                                queues[1].push(y);
                            }
                        }
                    }
                    result = Some(queues[1].clone());
                } else {
                    self.n_cp_toobig += 1;
                }
            } else {
                // v side complete: B enumerated in queues[1]
                if 2 * queues[1].len() <= n {
                    result = Some(queues[1].clone());
                } else {
                    self.n_cp_toobig += 1;
                }
            }
        } else {
            self.n_cp_notbridge += 1;
        }
        for f in 0..2 {
            for &v in queues[f].iter() {
                fb[f][(v >> 6) as usize] &= !(1u64 << (v & 63));
            }
            queues[f].clear();
        }
        self.fqueues = queues;
        self.fb = fb;
        result
    }

    /// Nonlocal cut-and-paste move: detach the branch hanging off a bridge
    /// edge, apply a random lattice symmetry, reattach at a uniform
    /// (cell, direction) of the remainder requiring exactly one contact.
    /// Selection probabilities are equal in both directions (a bridge is the
    /// unique edge between the sides, so the adjacency count is conserved),
    /// making the kernel symmetric with no Hastings factor.
    fn cutpaste(&mut self, rng: &mut Rng) -> bool {
        let n = self.cells.len();
        let ui = rng.below(n as u64) as usize;
        let d0 = rng.below(4) as usize;
        let vi_ = self.nbr[ui][d0];
        if vi_ == Self::NONE {
            self.n_cp_nopair += 1;
            return false;
        }
        let vi = vi_ as usize;
        let Some(bidx) = self.bridge_component(ui, vi) else {
            return false;
        };
        let w0 = self.w as i64;
        let pivot = self.cells[vi] as i64;
        let (pvx, pvy) = (pivot % w0, pivot / w0);
        let bpos: Vec<u32> = bidx.iter().map(|&i| self.cells[i as usize]).collect();
        let brel: Vec<(i64, i64)> = bpos
            .iter()
            .map(|&p| ((p as i64 % w0) - pvx, (p as i64 / w0) - pvy))
            .collect();
        let g = rng.below(8);
        let gb: Vec<(i64, i64)> = brel
            .iter()
            .map(|&(x, y)| match g {
                0 => (x, y),
                1 => (-y, x),
                2 => (-x, -y),
                3 => (y, -x),
                4 => (x, -y),
                5 => (-x, y),
                6 => (y, x),
                _ => (-y, -x),
            })
            .collect();
        // detach B
        for &p in &bpos {
            let ci = self.cell_idx[p as usize] as usize - 1;
            self.remove_cell(ci, p);
        }
        // attachment choice: uniform cell of the remainder + direction
        let upi = rng.below(self.cells.len() as u64) as usize;
        let dp = rng.below(4) as usize;
        // grid may need to grow to hold the placement; a rebuild shifts all
        // coordinates, so track the shift and re-derive positions by index
        let ext = gb
            .iter()
            .map(|&(x, y)| x.abs().max(y.abs()))
            .max()
            .unwrap() as usize;
        let mut shift = (0i64, 0i64);
        let mut wcur = self.w as i64;
        {
            let up = self.cells[upi] as i64;
            let (ux, uy) = (up % wcur, up / wcur);
            let dvec = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)][dp];
            let (tx, ty) = (ux + dvec.0, uy + dvec.1);
            let (minx, maxx) = (
                tx + gb.iter().map(|p| p.0).min().unwrap(),
                tx + gb.iter().map(|p| p.0).max().unwrap(),
            );
            let (miny, maxy) = (
                ty + gb.iter().map(|p| p.1).min().unwrap(),
                ty + gb.iter().map(|p| p.1).max().unwrap(),
            );
            if minx < 2
                || miny < 2
                || maxx >= self.w as i64 - 2
                || maxy >= self.h as i64 - 2
            {
                shift = self.rebuild_extra(ext + 4);
                wcur = self.w as i64;
            }
        }
        let up = self.cells[upi] as i64;
        let (ux, uy) = (up % wcur, up / wcur);
        let dvec = [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)][dp];
        let (tx, ty) = (ux + dvec.0, uy + dvec.1);
        // validity: every placed cell empty, and exactly one adjacency to the
        // remainder (the intended contact at u' is included in the count)
        let mut contacts = 0u32;
        let mut ok = true;
        for &(bx, by) in &gb {
            let (x, y) = (tx + bx, ty + by);
            let p = (y * wcur + x) as u32;
            if self.occ(p) {
                ok = false;
                break;
            }
            for dd in [1, -1, wcur, -wcur] {
                if self.occ((p as i64 + dd) as u32) {
                    contacts += 1;
                }
            }
        }
        if !ok || contacts != 1 {
            // reject: restore B at its original spot (shifted if rebuilt)
            self.n_cp_invalid += 1;
            let rebuilt = shift != (0, 0);
            let mut restored: Vec<u32> = Vec::with_capacity(bpos.len());
            for &p in &bpos {
                let (x, y) = ((p as i64 % w0) + shift.0, (p as i64 / w0) + shift.1);
                let np = (y * wcur + x) as u32;
                self.add_cell(np);
                restored.push(np);
            }
            if rebuilt {
                // perimeter was rebuilt for the remainder only; patch around B
                self.fix_perimeter_around(&restored);
            }
            return false;
        }
        // accept: place the transformed branch
        let mut placed: Vec<u32> = Vec::with_capacity(gb.len());
        for &(bx, by) in &gb {
            let p = ((ty + by) * wcur + (tx + bx)) as u32;
            self.add_cell(p);
            placed.push(p);
        }
        let rebuilt = shift != (0, 0);
        if !rebuilt {
            let old_sites: Vec<u32> = bpos.clone();
            self.fix_perimeter_around(&old_sites);
        }
        self.fix_perimeter_around(&placed);
        self.n_cp_accept += 1;
        true
    }

    /// Repair the perimeter indexed set around the given positions: each
    /// position and its neighbors get their membership recomputed from the
    /// final occupancy.  Idempotent.
    fn fix_perimeter_around(&mut self, sites: &[u32]) {
        let w = self.w as i64;
        for &s in sites {
            let sp = s as i64;
            for q in [sp, sp + 1, sp - 1, sp + w, sp - w] {
                let qp = q as u32;
                let should = !self.occ(qp) && self.has_occupied_neighbor(qp);
                let is = self.per_idx[qp as usize] != 0;
                if should && !is {
                    self.per_add(qp);
                } else if !should && is {
                    self.per_remove(qp);
                }
            }
        }
    }

    /// Mixed kernel: with probability 1/cp_inv try a cut-and-paste move,
    /// else a single-cell move.  cp_inv = 0 disables cut-and-paste.
    #[inline(always)]
    fn step_mixed(&mut self, rng: &mut Rng, cp_inv: u64) -> bool {
        if cp_inv > 0 && rng.below(cp_inv) == 0 {
            self.cutpaste(rng)
        } else {
            self.step(rng)
        }
    }

    /// Obviously-correct bridge reference: full BFS from u avoiding the
    /// (u, v) edge; the edge is a bridge iff v is not reached.  Also returns
    /// the v-side size when it is a bridge (n minus u-side size).
    fn bridge_ref(&mut self, ui: usize, vi: usize) -> Option<usize> {
        let n = self.cells.len();
        let mut vis = vec![0u64; (self.w * self.h + 63) / 64];
        let mark = |vis: &mut Vec<u64>, p: u32| {
            let (wi, bi) = ((p >> 6) as usize, p & 63);
            let old = vis[wi] >> bi & 1;
            vis[wi] |= 1 << bi;
            old == 0
        };
        let seen = |vis: &Vec<u64>, p: u32| vis[(p >> 6) as usize] >> (p & 63) & 1 != 0;
        let upos = self.cells[ui];
        self.queue.clear();
        self.queue.push(upos);
        mark(&mut vis, upos);
        let vpos = self.cells[vi];
        let w = self.w as i64;
        let mut head = 0;
        let mut cnt = 1usize;
        while head < self.queue.len() {
            let x = self.queue[head] as i64;
            head += 1;
            for d in [1, -1, w, -w] {
                let y = (x + d) as u32;
                if (x as u32 == upos && y == vpos) || (x as u32 == vpos && y == upos) {
                    continue;
                }
                if self.occ(y) && mark(&mut vis, y) {
                    cnt += 1;
                    self.queue.push(y);
                }
            }
        }
        if seen(&vis, vpos) {
            None
        } else {
            Some(n - cnt)
        }
    }

    /// Obviously-correct connectivity reference: full BFS over A - {c}.
    /// Used by selftest to differential-test `removable`.
    fn removable_ref(&mut self, c: u32) -> bool {
        let n = self.cells.len();
        let mut vis = vec![0u64; (self.w * self.h + 63) / 64];
        let mark = |vis: &mut Vec<u64>, p: u32| {
            let (wi, bi) = ((p >> 6) as usize, p & 63);
            let old = vis[wi] >> bi & 1;
            vis[wi] |= 1 << bi;
            old == 0
        };
        let w = self.w as i64;
        let start = if self.cells[0] == c { self.cells[1] } else { self.cells[0] };
        self.queue.clear();
        self.queue.push(start);
        mark(&mut vis, start);
        let mut head = 0;
        let mut cnt = 1;
        while head < self.queue.len() {
            let u = self.queue[head] as i64;
            head += 1;
            for d in [1, -1, w, -w] {
                let v = (u + d) as u32;
                if v != c && self.occ(v) && mark(&mut vis, v) {
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
        self.step_traced(rng).is_some()
    }

    /// One step; returns Some((c_pos, s_pos)) iff the state changed (both in
    /// the CURRENT grid frame, i.e. after any rebuild this step triggered).
    #[inline(always)]
    fn step_traced(&mut self, rng: &mut Rng) -> Option<(u32, u32)> {
        let n = self.cells.len();
        let ci = rng.below(n as u64) as usize;
        let c = self.cells[ci];
        if !self.removable(ci) {
            self.n_cut_reject += 1;
            return None;
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
        let mut cxy = ((c as usize % self.w) as i64, (c as usize / self.w) as i64);
        let mut sxy = (x as i64, y as i64);
        if x < 2 || x >= self.w - 2 || y < 2 || y >= self.h - 2 {
            let sh = self.rebuild();
            cxy = (cxy.0 + sh.0, cxy.1 + sh.1);
            sxy = (sxy.0 + sh.0, sxy.1 + sh.1);
        }
        let w = self.w as i64;
        if s != c {
            self.n_moved += 1;
            Some(((cxy.1 * w + cxy.0) as u32, (sxy.1 * w + sxy.0) as u32))
        } else {
            self.n_noop += 1;
            None
        }
    }

    /// Recentre and resize the grid; recompute perimeter and sums.
    /// Returns the (dx, dy) shift applied to every cell's (x, y).
    fn rebuild(&mut self) -> (i64, i64) {
        self.rebuild_extra(0)
    }

    fn rebuild_extra(&mut self, extra: usize) -> (i64, i64) {
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
        let mx = (bw / 4).max(16) + extra;
        let my = (bh / 4).max(16) + extra;
        let (nw, nh) = (bw + 2 * mx, bh + 2 * my);
        let old: Vec<u32> = std::mem::take(&mut self.cells);
        self.w = nw;
        self.h = nh;
        self.occ_bits = vec![0; (nw * nh + 63) / 64];
        self.cell_idx = vec![0; nw * nh];
        self.per_idx = vec![0; nw * nh];
        self.nbr.clear();
        self.e_edges = 0;
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
        let sh = (mx as i64 - minx as i64, my as i64 - miny as i64);
        self.shift_x += sh.0;
        self.shift_y += sh.1;
        sh
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

    /// Translation-canonical packed form for n <= 16 (4 bits per coordinate,
    /// one byte per cell, sorted, folded into a u128).
    fn canonical_u128(&self) -> u128 {
        let cs = self.coords();
        let minx = cs.iter().map(|c| c.0).min().unwrap();
        let miny = cs.iter().map(|c| c.1).min().unwrap();
        let mut bytes: Vec<u8> = cs
            .iter()
            .map(|c| (((c.0 - minx) as u8) << 4) | ((c.1 - miny) as u8))
            .collect();
        bytes.sort_unstable();
        bytes.iter().fold(0u128, |acc, &b| (acc << 8) | b as u128)
    }

    /// Full from-scratch consistency check (differential test vs the
    /// incremental structures).  Panics on any mismatch.
    fn check_invariants(&mut self, n: usize) {
        assert_eq!(self.cells.len(), n, "cell count");
        // no duplicate cells; index map consistent
        for (i, &p) in self.cells.iter().enumerate() {
            assert_eq!(self.cell_idx[p as usize], i as u32 + 1, "cell_idx");
        }
        // connectivity: BFS over all cells (local bitmap visited)
        let words = (self.w * self.h + 63) / 64;
        let mut vis = vec![0u64; words];
        let mark = |vis: &mut Vec<u64>, p: u32| {
            let (wi, bi) = ((p >> 6) as usize, p & 63);
            let old = vis[wi] >> bi & 1;
            vis[wi] |= 1 << bi;
            old == 0
        };
        let w = self.w as i64;
        self.queue.clear();
        self.queue.push(self.cells[0]);
        mark(&mut vis, self.cells[0]);
        let mut head = 0;
        let mut cnt = 1;
        while head < self.queue.len() {
            let u = self.queue[head] as i64;
            head += 1;
            for d in [1, -1, w, -w] {
                let v = (u + d) as u32;
                if self.occ(v) && mark(&mut vis, v) {
                    cnt += 1;
                    self.queue.push(v);
                }
            }
        }
        assert_eq!(cnt, n, "connectivity");
        // perimeter set == from-scratch recompute
        let mut fresh: Vec<u32> = Vec::new();
        let mut vis2 = vec![0u64; words];
        for i in 0..self.cells.len() {
            let c = self.cells[i] as i64;
            for d in [1, -1, w, -w] {
                let q = (c + d) as u32;
                if !self.occ(q) && mark(&mut vis2, q) {
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
        // edge count
        let mut e = 0u64;
        for i in 0..self.cells.len() {
            e += self.nbr[i].iter().filter(|&&l| l != Self::NONE).count() as u64;
        }
        // (degree sum = 2E; add_cell counts each edge once, at whichever
        // endpoint arrived second)
        assert_eq!(e % 2, 0, "odd degree sum");
        assert_eq!(self.e_edges, e / 2, "edge count");
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
        let total = self.n_cut_reject + self.n_noop + self.n_moved
            + self.n_cp_accept + self.n_cp_invalid + self.n_cp_notbridge
            + self.n_cp_toobig + self.n_cp_nopair + self.n_cp_capped;
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
        ) + &format!(
            " cp[accept={} invalid={} notbridge={} toobig={} nopair={} capped={}]",
            self.n_cp_accept, self.n_cp_invalid, self.n_cp_notbridge,
            self.n_cp_toobig, self.n_cp_nopair, self.n_cp_capped,
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
            let cp_inv: u64 = get(&args, "cp-inv", 0);
            let mut st = State::new(n, &init);
            st.cp_cap = get(&args, "cp-cap", 4000);
            let t0 = std::time::Instant::now();
            for _ in 0..steps {
                st.step_mixed(&mut rng, cp_inv);
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
            assert!(n <= 16, "canonical_u128 packing supports n <= 16");
            let obs: u64 = get(&args, "obs", 1_000_000);
            let thin: u64 = get(&args, "thin", 10);
            let burn: u64 = get(&args, "burn", 100_000);
            let check_every: u64 = get(&args, "check-every", 0);
            let cp_inv: u64 = get(&args, "cp-inv", 0);
            let mut st = State::new(n, &init);
            st.cp_cap = get(&args, "cp-cap", 4000);
            for _ in 0..burn {
                st.step_mixed(&mut rng, cp_inv);
            }
            let mut counts: HashMap<u128, u64> = HashMap::new();
            for i in 0..obs {
                for _ in 0..thin {
                    st.step_mixed(&mut rng, cp_inv);
                }
                *counts.entry(st.canonical_u128()).or_insert(0) += 1;
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
            let cp_inv: u64 = get(&args, "cp-inv", 0);
            let out_path = args.get("out").cloned();
            let dump_path = args.get("dump-final").cloned();
            let dump_every: u64 = get(&args, "dump-every", 0);
            let dump_prefix = args.get("dump-prefix").cloned();
            let mut st = match args.get("init-file") {
                Some(p) => State::from_file(p),
                None => State::new(n, &init),
            };
            st.cp_cap = get(&args, "cp-cap", 4000);
            let n = st.cells.len();
            st.check_invariants(n);
            let mut writer = out_path.map(|p| {
                let f = std::fs::File::create(p).expect("create out");
                let mut w = std::io::BufWriter::new(f);
                writeln!(w, "step,rg2,perim,cycles").unwrap();
                w
            });
            let t0 = std::time::Instant::now();
            let mut done: u64 = 0;
            while done < steps {
                let chunk = record_every.min(steps - done);
                for _ in 0..chunk {
                    st.step_mixed(&mut rng, cp_inv);
                }
                done += chunk;
                if let Some(w) = writer.as_mut() {
                    let cyc = st.e_edges as i64 - st.cells.len() as i64 + 1;
                    writeln!(w, "{},{:.6},{},{}", done, st.rg2(), st.per.len(), cyc)
                        .unwrap();
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
        "dbtest2" => {
            // End-to-end empirical measurement of the one-step transition
            // probability P(A -> A') for a specific pair, using the real
            // step semantics: each trial draws the uniform cell index; only
            // when it hits c* can the step reach A', in which case the full
            // production path (connectivity check, removal, incremental
            // perimeter, uniform site draw) runs and the destination is
            // compared against s*.  Reports hits / trials.
            let file = args.get("init-file").expect("--init-file");
            let trials: u64 = get(&args, "trials", 1_200_000_000_000);
            let parse_xy = |s: &str| -> (i64, i64) {
                let mut it = s.split(',');
                (
                    it.next().unwrap().trim().parse().unwrap(),
                    it.next().unwrap().trim().parse().unwrap(),
                )
            };
            let (cx, cy) = parse_xy(args.get("c").expect("--c"));
            let (sx, sy) = parse_xy(args.get("s").expect("--s"));
            let text = std::fs::read_to_string(file).expect("read init file");
            let mut pts: Vec<(i64, i64)> = Vec::new();
            for line in text.lines() {
                let mut it = line.split_whitespace();
                pts.push((
                    it.next().unwrap().parse().unwrap(),
                    it.next().unwrap().parse().unwrap(),
                ));
            }
            let minx = pts.iter().map(|p| p.0).min().unwrap();
            let miny = pts.iter().map(|p| p.1).min().unwrap();
            let mut st = State::from_file(file);
            let n = st.cells.len();
            st.check_invariants(n);
            let first = st.cells[0] as i64;
            let w = st.w as i64;
            let (fx, fy) = (first % w, first / w);
            let (offx, offy) = (fx - (pts[0].0 - minx), fy - (pts[0].1 - miny));
            let cpos = ((cy - miny + offy) * w + (cx - minx + offx)) as u32;
            let spos = ((sy - miny + offy) * w + (sx - minx + offx)) as u32;
            assert!(st.occ(cpos) && !st.occ(spos), "move endpoints inconsistent");
            let mut hits = 0u64;
            let mut cell_hits = 0u64;
            let nn = n as u64;
            for t in 0..trials {
                let ci = rng.below(nn) as usize;
                // a step that picks any other cell cannot produce A'
                if st.cells[ci] != cpos {
                    continue;
                }
                cell_hits += 1;
                // full production path from here
                assert!(st.removable(ci), "c* must be removable");
                st.remove_cell(ci, cpos);
                st.per_add(cpos);
                for d in [1, -1, w, -w] {
                    let q = (cpos as i64 + d) as u32;
                    if !st.occ(q) && st.per_idx[q as usize] != 0 && !st.has_occupied_neighbor(q)
                    {
                        st.per_remove(q);
                    }
                }
                let s = st.per[rng.below(st.per.len() as u64) as usize];
                if s == spos {
                    hits += 1;
                }
                // restore A via the production placement path (put c* back)
                st.per_remove(cpos);
                st.add_cell(cpos);
                for d in [1, -1, w, -w] {
                    let q = (cpos as i64 + d) as u32;
                    if !st.occ(q) && st.per_idx[q as usize] == 0 {
                        st.per_add(q);
                    }
                }
                if cell_hits % 10_000_000 == 0 {
                    st.check_invariants(n);
                }
                let _ = t;
            }
            st.check_invariants(n);
            println!(
                "trials={} cell_hits={} (expect ~{}) transition_hits={} P_hat={:.6e}",
                trials,
                cell_hits,
                trials / nn,
                hits,
                hits as f64 / trials as f64
            );
        }
        "dbtest" => {
            // Audit the sampler implementation on a specific state/transition:
            //  A. exhaustive removability: production removable() vs full-BFS
            //     reference for EVERY cell (exact, no sampling);
            //  B. after the production removal path for --c, the incremental
            //     perimeter equals the from-scratch perimeter of A - c as a
            //     set (each site exactly once) and contains --s;
            //  C. cell-pick: M raw draws through rng.below(n) over the cells
            //     array (chi-square emitted for Python);
            //  D. M forced-move trials: production removal + uniform site
            //     draw + production placement, then production-primitive
            //     undo; destination counts per site emitted for Python.
            let file = args.get("init-file").expect("--init-file");
            let trials: u64 = get(&args, "trials", 6_000_000);
            let draws: u64 = get(&args, "draws", 100_000_000);
            let out = args.get("out").expect("--out").clone();
            let parse_xy = |s: &str| -> (i64, i64) {
                let mut it = s.split(',');
                (
                    it.next().unwrap().trim().parse().unwrap(),
                    it.next().unwrap().trim().parse().unwrap(),
                )
            };
            let (cx, cy) = parse_xy(args.get("c").expect("--c"));
            let (sx, sy) = parse_xy(args.get("s").expect("--s"));
            let text = std::fs::read_to_string(file).expect("read init file");
            let mut pts: Vec<(i64, i64)> = Vec::new();
            for line in text.lines() {
                let mut it = line.split_whitespace();
                pts.push((
                    it.next().unwrap().parse().unwrap(),
                    it.next().unwrap().parse().unwrap(),
                ));
            }
            let minx = pts.iter().map(|p| p.0).min().unwrap();
            let miny = pts.iter().map(|p| p.1).min().unwrap();
            let mut st = State::from_file(file);
            let n = st.cells.len();
            st.check_invariants(n);
            let first = st.cells[0] as i64;
            let w = st.w as i64;
            let (fx, fy) = (first % w, first / w);
            let (offx, offy) = (fx - (pts[0].0 - minx), fy - (pts[0].1 - miny));
            let cpos = ((cy - miny + offy) * w + (cx - minx + offx)) as u32;
            let spos = ((sy - miny + offy) * w + (sx - minx + offx)) as u32;
            assert!(st.occ(cpos) && !st.occ(spos), "move endpoints inconsistent");

            // ---- A. exhaustive removability boundary
            let mut n_removable = 0u64;
            let mut mismatches = 0u64;
            for ci in 0..n {
                let fast = st.removable(ci);
                let refr = st.removable_ref(st.cells[ci]);
                if fast != refr {
                    mismatches += 1;
                }
                if refr {
                    n_removable += 1;
                }
            }
            println!(
                "A: exhaustive removability over {n} cells: {} removable, {} mismatches",
                n_removable, mismatches
            );
            assert_eq!(mismatches, 0);

            // ---- B. perimeter set equality after production removal of c
            let ci = st.cell_idx[cpos as usize] as usize - 1;
            assert!(st.removable(ci), "--c not removable");
            st.remove_cell(ci, cpos);
            st.per_add(cpos);
            for d in [1, -1, w, -w] {
                let q = (cpos as i64 + d) as u32;
                if !st.occ(q) && st.per_idx[q as usize] != 0 && !st.has_occupied_neighbor(q) {
                    st.per_remove(q);
                }
            }
            // from-scratch perimeter of A - c
            let mut fresh: std::collections::HashSet<u32> = std::collections::HashSet::new();
            for &p in &st.cells {
                for d in [1, -1, w, -w] {
                    let q = (p as i64 + d) as u32;
                    if !st.occ(q) {
                        fresh.insert(q);
                    }
                }
            }
            let mut per_set: std::collections::HashSet<u32> = std::collections::HashSet::new();
            for &p in &st.per {
                assert!(per_set.insert(p), "duplicate perimeter entry");
            }
            assert_eq!(per_set, fresh, "incremental perimeter != from-scratch");
            assert!(per_set.contains(&spos), "--s not in perimeter of A - c");
            let plen = st.per.len();
            println!(
                "B: incremental perimeter of A-c == from-scratch ({} sites, each exactly once; s present)",
                plen
            );
            // restore c via the production placement path
            let ci_restore = st.per_idx[cpos as usize];
            assert!(ci_restore != 0);
            st.per_remove(cpos);
            st.add_cell(cpos);
            for d in [1, -1, w, -w] {
                let q = (cpos as i64 + d) as u32;
                if !st.occ(q) && st.per_idx[q as usize] == 0 {
                    st.per_add(q);
                }
            }
            st.check_invariants(n);

            // ---- C. raw cell-pick distribution
            let mut cell_hist = vec![0u64; n];
            for _ in 0..draws {
                cell_hist[rng.below(n as u64) as usize] += 1;
            }
            let e = draws as f64 / n as f64;
            let chi2: f64 = cell_hist
                .iter()
                .map(|&o| (o as f64 - e) * (o as f64 - e) / e)
                .sum();
            println!("C: cell-pick draws={} chi2={:.1} dof={}", draws, chi2, n - 1);

            // ---- D. forced-move destination distribution
            let mut counts: HashMap<u32, u64> = HashMap::new();
            let ci0 = st.cell_idx[cpos as usize] as usize - 1;
            let _ = ci0;
            for t in 0..trials {
                let ci = st.cell_idx[cpos as usize] as usize - 1;
                st.remove_cell(ci, cpos);
                st.per_add(cpos);
                for d in [1, -1, w, -w] {
                    let q = (cpos as i64 + d) as u32;
                    if !st.occ(q) && st.per_idx[q as usize] != 0 && !st.has_occupied_neighbor(q)
                    {
                        st.per_remove(q);
                    }
                }
                let s = st.per[rng.below(st.per.len() as u64) as usize];
                *counts.entry(s).or_insert(0) += 1;
                st.per_remove(s);
                st.add_cell(s);
                for d in [1, -1, w, -w] {
                    let q = (s as i64 + d) as u32;
                    if !st.occ(q) && st.per_idx[q as usize] == 0 {
                        st.per_add(q);
                    }
                }
                // undo with the same production primitives (move s -> c)
                if s != cpos {
                    let si = st.cell_idx[s as usize] as usize - 1;
                    st.remove_cell(si, s);
                    st.per_add(s);
                    for d in [1, -1, w, -w] {
                        let q = (s as i64 + d) as u32;
                        if !st.occ(q)
                            && st.per_idx[q as usize] != 0
                            && !st.has_occupied_neighbor(q)
                        {
                            st.per_remove(q);
                        }
                    }
                    st.per_remove(cpos);
                    st.add_cell(cpos);
                    for d in [1, -1, w, -w] {
                        let q = (cpos as i64 + d) as u32;
                        if !st.occ(q) && st.per_idx[q as usize] == 0 {
                            st.per_add(q);
                        }
                    }
                }
                if t % 1_000_000 == 0 {
                    st.check_invariants(n);
                }
            }
            st.check_invariants(n);
            let f = std::fs::File::create(&out).expect("create out");
            let mut wtr = std::io::BufWriter::new(f);
            writeln!(wtr, "# n {n} plen {plen} trials {trials}").unwrap();
            for (&p, &cnt) in &counts {
                let (x, y) = ((p as i64 % w) - offx + minx, (p as i64 / w) - offy + miny);
                writeln!(wtr, "{x} {y} {cnt}").unwrap();
            }
            wtr.flush().unwrap();
            println!("D: {} forced-move trials over {} sites written to {}", trials, counts.len(), out);
        }
        "cyclewatch" => {
            // Run the pure single-cell chain from an equilibrated snapshot
            // until an accepted move closes a cycle enclosing a bounded empty
            // region of >= --threshold cells; dump the post-move state and
            // the move.  Keep running until the enclosure dies (ring breach
            // or filled); dump that transition too.  Python then verifies
            // detailed balance for both transitions from scratch.
            let file = args.get("init-file").expect("--init-file");
            let threshold: usize = get(&args, "threshold", 100);
            let max_steps: u64 = get(&args, "max-steps", 4_000_000_000);
            let prefix = args.get("out-prefix").expect("--out-prefix").clone();
            let mut st = State::from_file(file);
            let n = st.cells.len();
            st.check_invariants(n);
            let cap = 20_000usize;

            // bounded flood over empty cells from `start`; None if the region
            // exceeds `cap` (i.e. it is the outside, not a hole)
            let flood = |st: &State, start: u32| -> Option<Vec<u32>> {
                if st.occ(start) {
                    return Some(Vec::new());
                }
                let w = st.w as i64;
                let h = st.h as i64;
                let mut reg: Vec<u32> = vec![start];
                let mut seen = std::collections::HashSet::new();
                seen.insert(start);
                let mut head = 0;
                while head < reg.len() {
                    let u = reg[head] as i64;
                    head += 1;
                    let (x, y) = (u % w, u / w);
                    if x <= 1 || x >= w - 2 || y <= 1 || y >= h - 2 {
                        return None; // reached the border: this is the outside
                    }
                    for d in [1, -1, w, -w] {
                        let v = (u + d) as u32;
                        if !st.occ(v) && seen.insert(v) {
                            reg.push(v);
                            if reg.len() > cap {
                                return None;
                            }
                        }
                    }
                }
                Some(reg)
            };
            let dump = |st: &State, path: &str| {
                let f = std::fs::File::create(path).expect("create dump");
                let mut w = std::io::BufWriter::new(f);
                for (x, y) in st.coords() {
                    writeln!(w, "{x} {y}").unwrap();
                }
                w.flush().unwrap();
            };
            let xy = |st: &State, p: u32| ((p as usize % st.w) as i64, (p as usize / st.w) as i64);

            // phase 1: wait for a big-cycle birth
            let mut hole: Vec<u32> = Vec::new();
            let mut birth_step = 0u64;
            let mut t = 0u64;
            'search: while t < max_steps {
                t += 1;
                let Some((c, s)) = st.step_traced(&mut rng) else {
                    continue;
                };
                // gap fill? s has >= 2 occupied neighbors
                let w = st.w as i64;
                let mut occn = 0;
                for d in [1, -1, w, -w] {
                    if st.occ((s as i64 + d) as u32) {
                        occn += 1;
                    }
                }
                if occn < 2 {
                    continue;
                }
                for d in [1, -1, w, -w] {
                    let e = (s as i64 + d) as u32;
                    if st.occ(e) {
                        continue;
                    }
                    if let Some(reg) = flood(&st, e) {
                        // genuine cycle birth only: the region must have been
                        // OPEN before this move, i.e. flooding with s treated
                        // as empty (and c occupied) must escape the cap
                        let open_before = {
                            let mut reg2: Vec<u32> = vec![e];
                            let mut seen = std::collections::HashSet::new();
                            seen.insert(e);
                            let mut head = 0;
                            let mut escaped = false;
                            let h = st.h as i64;
                            while head < reg2.len() {
                                let u = reg2[head] as i64;
                                head += 1;
                                let (x, y) = (u % w, u / w);
                                if x <= 1 || x >= w - 2 || y <= 1 || y >= h - 2 {
                                    escaped = true;
                                    break;
                                }
                                for dd in [1, -1, w, -w] {
                                    let v = (u + dd) as u32;
                                    let empty_before =
                                        (v == s || !st.occ(v)) && v != c;
                                    if empty_before && seen.insert(v) {
                                        reg2.push(v);
                                        if reg2.len() > cap {
                                            escaped = true;
                                            break;
                                        }
                                    }
                                }
                                if escaped {
                                    break;
                                }
                            }
                            escaped
                        };
                        if open_before && reg.len() >= threshold {
                            let (cx, cy) = xy(&st, c);
                            let (sx2, sy2) = xy(&st, s);
                            dump(&st, &format!("{prefix}_birth_after.txt"));
                            let mut mf = std::fs::File::create(format!("{prefix}_birth_meta.txt")).unwrap();
                            writeln!(mf, "step {t}\nc {cx} {cy}\ns {sx2} {sy2}\narea {}", reg.len()).unwrap();
                            eprintln!("birth at step {t}: hole area {} (c=({cx},{cy}) -> s=({sx2},{sy2}))", reg.len());
                            hole = reg;
                            birth_step = t;
                            break 'search;
                        }
                    }
                }
            }
            assert!(!hole.is_empty(), "no big cycle within max-steps");

            // phase 2: watch until the enclosure dies.  Track the hole as
            // (x, y) pairs so grid rebuilds (which shift all coordinates by a
            // known delta) are easy to follow.
            let mut hole_xy: Vec<(i64, i64)> = {
                let w = st.w as i64;
                hole.iter().map(|&p| (p as i64 % w, p as i64 / w)).collect()
            };
            let mut shift0 = (st.shift_x, st.shift_y);
            loop {
                t += 1;
                let Some((c, s)) = st.step_traced(&mut rng) else {
                    continue;
                };
                if (st.shift_x, st.shift_y) != shift0 {
                    let (dx, dy) = (st.shift_x - shift0.0, st.shift_y - shift0.1);
                    for p in hole_xy.iter_mut() {
                        p.0 += dx;
                        p.1 += dy;
                    }
                    shift0 = (st.shift_x, st.shift_y);
                    eprintln!("rebuild during watch at step {t}; hole translated by ({dx},{dy})");
                }
                let w = st.w as i64;
                // refresh the hole from any still-empty member
                let seed = hole_xy
                    .iter()
                    .map(|&(x, y)| (y * w + x) as u32)
                    .find(|&p| !st.occ(p));
                let Some(seed) = seed else {
                    let (cx, cy) = xy(&st, c);
                    let (sx2, sy2) = xy(&st, s);
                    dump(&st, &format!("{prefix}_death_after.txt"));
                    let mut mf =
                        std::fs::File::create(format!("{prefix}_death_meta.txt")).unwrap();
                    writeln!(mf, "step {t}\nc {cx} {cy}\ns {sx2} {sy2}\ncause filled").unwrap();
                    eprintln!(
                        "death (filled) at step {t} after {} steps of life",
                        t - birth_step
                    );
                    break;
                };
                match flood(&st, seed) {
                    Some(reg) => {
                        hole_xy = reg
                            .iter()
                            .map(|&p| (p as i64 % w, p as i64 / w))
                            .collect();
                    }
                    None => {
                        let (cx, cy) = xy(&st, c);
                        let (sx2, sy2) = xy(&st, s);
                        dump(&st, &format!("{prefix}_death_after.txt"));
                        let mut mf =
                            std::fs::File::create(format!("{prefix}_death_meta.txt")).unwrap();
                        writeln!(mf, "step {t}\nc {cx} {cy}\ns {sx2} {sy2}\ncause breach").unwrap();
                        eprintln!(
                            "death (breach) at step {t} after {} steps of life",
                            t - birth_step
                        );
                        break;
                    }
                }
            }
            st.check_invariants(n);
            eprintln!("done; {}", st.stats_line());
        }
        "probe" => {
            // Execute one specific single-cell transition on a loaded shape,
            // using the exact production move machinery: move the cell at
            // --leaf "x,y" to the empty site --target "x,y" (raw dump
            // coordinates).  Verifies: leaf passes the connectivity check,
            // target is a selectable perimeter site of A - leaf, and the
            // final state passes all invariants.  Prints the per-step
            // proposal probability of exactly this transition.
            let file = args.get("init-file").expect("--init-file");
            let parse_xy = |s: &str| -> (i64, i64) {
                let mut it = s.split(',');
                (
                    it.next().unwrap().trim().parse().unwrap(),
                    it.next().unwrap().trim().parse().unwrap(),
                )
            };
            let (lx, ly) = parse_xy(args.get("leaf").expect("--leaf"));
            let (tx, ty) = parse_xy(args.get("target").expect("--target"));
            // replicate from_file's translation so raw coords line up
            let text = std::fs::read_to_string(file).expect("read init file");
            let mut pts: Vec<(i64, i64)> = Vec::new();
            for line in text.lines() {
                let mut it = line.split_whitespace();
                pts.push((
                    it.next().unwrap().parse().unwrap(),
                    it.next().unwrap().parse().unwrap(),
                ));
            }
            let minx = pts.iter().map(|p| p.0).min().unwrap();
            let miny = pts.iter().map(|p| p.1).min().unwrap();
            let mut st = State::from_file(file);
            let n = st.cells.len();
            st.check_invariants(n);
            // from_file translates by (-minx,-miny); State::new adds margins
            // (mx,my) = position of cell (0,0).  Recover the offset from any
            // known cell: use the first point in the file.
            let p0 = st.cells[st.cell_idx
                .iter()
                .position(|&v| v == 1)
                .map(|_| 0usize)
                .unwrap_or(0)];
            let _ = p0;
            // offset = grid position of file point pts[0] minus its shifted xy
            let first = st.cells[0] as i64; // cells pushed in file order
            let (fx, fy) = (first % st.w as i64, first / st.w as i64);
            let (offx, offy) = (fx - (pts[0].0 - minx), fy - (pts[0].1 - miny));
            let w = st.w as i64;
            let lpos = ((ly - miny + offy) * w + (lx - minx + offx)) as u32;
            let tpos = ((ty - miny + offy) * w + (tx - minx + offx)) as u32;
            assert!(st.occ(lpos), "leaf position not occupied");
            assert!(!st.occ(tpos), "target position not empty");
            let ci = st.cell_idx[lpos as usize] as usize - 1;
            // 1. production connectivity check
            let removable = st.removable(ci);
            println!("leaf ({lx},{ly}) removable: {removable}");
            assert!(removable, "leaf failed the connectivity check");
            // 2. production removal + perimeter update (exactly as step())
            st.remove_cell(ci, lpos);
            st.per_add(lpos);
            for d in [1, -1, w, -w] {
                let q = (lpos as i64 + d) as u32;
                if !st.occ(q) && st.per_idx[q as usize] != 0 && !st.has_occupied_neighbor(q) {
                    st.per_remove(q);
                }
            }
            let per_len = st.per.len();
            let in_per = st.per_idx[tpos as usize] != 0;
            println!(
                "target ({tx},{ty}) in perimeter of A-leaf: {in_per}  (|P| = {per_len}, n = {n})"
            );
            assert!(in_per, "target is not a selectable perimeter site");
            // 3. complete the transition exactly as step() would
            st.per_remove(tpos);
            st.add_cell(tpos);
            for d in [1, -1, w, -w] {
                let q = (tpos as i64 + d) as u32;
                if !st.occ(q) && st.per_idx[q as usize] == 0 {
                    st.per_add(q);
                }
            }
            st.check_invariants(n);
            println!(
                "transition executed; all invariants pass; proposal probability = 1/{} * 1/{} = {:.3e}",
                n,
                per_len,
                1.0 / (n as f64) / (per_len as f64)
            );
        }
        "perm" => {
            // Independent PERM sampler: estimates a_k (animal counts) and
            // the weighted <Rg^2> at the target size.
            let tours: u64 = get(&args, "tours", 1_000_000);
            let p: f64 = get(&args, "p", 0.59);
            let mut pm = perm::Perm::new(n, p, seed, 16);
            let t0 = std::time::Instant::now();
            for _ in 0..tours {
                pm.tour(tours);
            }
            let dt = t0.elapsed().as_secs_f64();
            println!(
                "perm n={} p={} tours={} samples_at_n={} secs={:.1}",
                n, p, tours, pm.samples, dt
            );
            for k in 1..=n.min(10) {
                println!("a_{k} estimate: {:.4e}", pm.log_ak(k).exp());
            }
            if let Some(path) = args.get("ak-out") {
                let f = std::fs::File::create(path).expect("create ak-out");
                let mut wtr = std::io::BufWriter::new(f);
                writeln!(wtr, "k log_ak").unwrap();
                for k in 1..=n {
                    writeln!(wtr, "{k} {:.8}", pm.log_ak(k)).unwrap();
                }
                wtr.flush().unwrap();
            }
            if n > 2 {
                let lam = (pm.log_ak(n) - pm.log_ak(n - 1)
                    + ((n as f64) / (n as f64 - 1.0)).ln() * 1.0)
                    .exp()
                    * (n as f64 / (n as f64 - 1.0)).powf(0.0);
                println!("lambda_hat (a_n/a_(n-1)): {:.4}", lam);
            }
            let (rg2, se) = pm.rg2_mean_stderr();
            println!(
                "<Rg^2>_{n} = {:.4} +- {:.4}   weight_ESS = {:.1} (of {} samples)",
                rg2, se, pm.weight_ess(), pm.samples
            );
        }
        "selftest" => {
            // heavy invariant checking at several sizes
            for (nn, steps, every) in [(2usize, 200_000u64, 100u64), (3, 200_000, 100),
                                       (5, 200_000, 100), (10, 200_000, 200),
                                       (37, 100_000, 500), (100, 100_000, 1000),
                                       (1000, 50_000, 5000)] {
                for init in ["bar", "rect", "bar2", "bar8"] {
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
                // mixed kernel (25% cut-and-paste) + bridge differential
                for init in ["bar", "rect"] {
                    let mut st = State::new(nn, init);
                    let mut r = Rng::new(seed ^ (nn as u64) ^ 0xC0FFEE);
                    for i in 0..steps / 2 {
                        st.step_mixed(&mut r, 4);
                        if i % every == 0 {
                            st.check_invariants(nn);
                            for _ in 0..nn.min(32) {
                                let ui = r.below(nn as u64) as usize;
                                let d0 = r.below(4) as usize;
                                let vi = st.nbr[ui][d0];
                                if vi == State::NONE {
                                    continue;
                                }
                                let fast = st.bridge_component(ui, vi as usize);
                                let refr = st.bridge_ref(ui, vi as usize);
                                match refr {
                                    None => assert!(fast.is_none(), "bridge false pos n={nn}"),
                                    Some(bsz) => {
                                        if 2 * bsz <= nn {
                                            let f = fast.expect("bridge missed");
                                            assert_eq!(f.len(), bsz, "bridge size n={nn}");
                                        } else {
                                            assert!(fast.is_none(), "toobig not rejected");
                                        }
                                    }
                                }
                            }
                        }
                    }
                    st.check_invariants(nn);
                    println!("selftest-mixed n={nn} init={init} OK  {}", st.stats_line());
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
