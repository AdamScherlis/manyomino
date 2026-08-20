# manyomino

Uniformly sample extremely large polyominoes.

Samples polyominoes (lattice animals on Z², 4-connectivity) **uniformly at
random for fixed size n**, from n ~ 100 up to n = 100,000+, renders them, and
proves — via a validation suite — that the samples really are uniform and the
chain really has mixed.  Target distribution is uniform over **fixed**
polyominoes (translation classes, OEIS
[A001168](https://oeis.org/A001168)).  See `docs/HANDOFF.md` for the full
specification this implements.

Gallery of converged samples:
<https://adamscherlis.github.io/static/manyomino/>

## The chain

State: a 4-connected set A ⊂ Z², |A| = n.  One move:

1. Pick a cell c uniformly from A.
2. If A − {c} is disconnected: stay put (still counts as a step).
3. Else pick s uniformly from the perimeter *sites* of A − {c} (c itself is
   always one of them, so s = c is a legal no-op).
4. New state: (A − {c}) ∪ {s}.

Forward and reverse proposal probabilities are equal (A − c = A′ − s), so the
kernel is symmetric and the uniform distribution is stationary; every legal
proposal is accepted.

## Layout

- `python/polyref.py` — slow, obviously-correct reference chain
  (from-scratch connectivity + perimeter every move).
- `python/enumerate_fixed.py` — brute-force enumerator of fixed polyominoes;
  matches A001168 for n = 1..10 (`python3 enumerate_fixed.py 10`).
- `python/validate_reference.py` — chi-square uniformity of the reference
  chain against the exact shape list (`python3 validate_reference.py 6
  1000000`).
- `rust/` — the fast sampler, std-only (no crates; the build environment has
  no network).  `cargo build --release`; binary `manyomino` with modes:
  - `selftest` — invariant checks (cell count, connectivity, perimeter set
    and adjacency links vs from-scratch recomputes) plus a differential test
    of the fast connectivity check against a full BFS, across sizes and both
    seeds.
  - `chi --n 5|6|8 --obs M --thin K --burn B` — emit canonical-shape counts
    for `python/validate_fast.py`.
  - `run --n N --init bar|rect --steps T --record-every R --out ts.csv
    [--dump-every D --dump-prefix P] [--check-every X]` — production runs:
    R_g² time series + periodic shape snapshots.
  - `bench` — moves/second.
- `python/analyze.py` — τ (Sokal-windowed, binned), ESS, means ± stderr for a
  time series; the stderr is the more conservative of the Sokal estimate and
  coarse batch means, since windowing misses slow modes.
- `python/twoseed.py` — two-seed convergence test (bar vs rect chain).
- `python/nufit.py` — weighted log-log fit of ⟨R_g⟩ vs n.
- `python/render.py` — stdlib PNG renderer (`bw` and graph-distance `dist`
  modes).
- `rust/src/perm.rs` — PERM (pruned-enriched Rosenbluth) sampler: Leath-style
  FIFO growth, log-space weights, adaptive pruning/enrichment; independent
  cross-check of the MCMC and estimator of the animal counts a_k
  (`manyomino perm --n N --tours T --p 0.4`).
- `python/make_page.py` — regenerates the gallery page from current runs;
  only publishes sizes whose two-seed test passes and whose snapshots are
  ≥ 5τ into their chain.

## Implementation notes

- Grid-backed indexed sets (dense array + position map, swap-remove) for O(1)
  uniform sampling of cells and perimeter *sites* (not (cell, direction)
  pairs — that would bias the chain).
- The perimeter used for the proposal is that of A − {c}: c is removed first,
  the perimeter is updated incrementally, then s is drawn.
- Connectivity test: leaf shortcut → 3×3 ring-arc test → alternating
  multi-front BFS over maintained adjacency links, one front per neighbor arc
  of c, merged with a tiny union-find.  Proving "cut" costs O(smallest
  component) ≈ O(√n) for branched-polymer shapes, not O(n).  At equilibrium
  ~70% of proposals are cut cells, so this is the hot path.
- The chain's centroid random-walks; the grid is rebuilt (recentred, resized)
  whenever a cell nears the border, and R_g² is tracked by O(1) coordinate
  sums.

## Validation status

| check | result |
|---|---|
| exact counts n=1..10 vs A001168 | exact match |
| reference chain chi-square, n=5 (63 shapes, 10⁶ obs) | p = 0.70 |
| reference chain chi-square, n=6 (216 shapes, 10⁶ obs) | p = 0.66 |
| fast sampler chi-square, n=5 (2×10⁶ obs) | p = 0.03 |
| fast sampler chi-square, n=6 (2×10⁶ obs) | p = 0.39 |
| fast sampler chi-square, n=8 (2725 shapes, 4×10⁶ obs) | p = 0.49 |
| move-type fractions, reference vs fast, n=5,6 | agree to ~10⁻³ |
| mixed-kernel chi-square, n=5 / n=6 / n=8 | p = 0.94 / 0.19 / 0.43 |
| cycle + near-miss stats vs exact census, n=12 (505,861 shapes) | all \|z\| ≤ 1.4 |
| two-seed convergence, n=1000 / 2000 / 3000 / 10000 | \|z\| = 0.06 / 0.49 / 0.09 / 0.67, all PASS |
| two-seed, n=5000, four independent chains pooled by init | \|z\| = 1.87, PASS |
| leaf-to-gap probe: largest almost-cycle (156/108 cells) in n=10⁴ snapshots | legal, probability (1/n)(1/\|P\|), invariants pass |
| detailed balance on captured cycle-seal/breach transitions (n=10⁴) | exact: \|P\| identical both ways; implementation audit: removability boundary exhaustive 0 mismatches, perimeter set == from-scratch, destination uniform (p=0.11–0.83) |
| end-to-end P(A→A′) vs P(A′→A), 1.2×10¹² unmodified-step trials each | seal: 10099 vs 10047 hits (z=+0.4); breach: 9840 vs 10118 (z=−2.0); all four consistent with exact 1/(n\|P\|) |
| ν fit, n = 100..100000 | 0.6441 ± 0.0008 stat; previously 0.6436 ± 0.0010 stat ± 0.001 sys (corrections-to-scaling scan Δ=0.5..1.5: ν=0.6436–0.6447; expect ≈ 0.6408) |
| PERM (independent non-MCMC sampler): a₁..a₁₀ vs A001168 | within ~0.5% statistical |
| PERM λ̂ | ≈ 4.01–4.08 (Klarner 4.0626) |
| PERM ⟨R_g²⟩ vs MCMC: n=10 (exact) / 100 / 300 / 1000 | 3.1794±27 vs 3.1801 exact; z = 1.0 / 0.4 / pooled 1380±37 vs 1339–1382 |
| bridge-cap symmetry (chi-square n=8, cap=2, cap fired 1.5×10⁶ times) | p = 0.885 |
| two-seed, n=30000 | \|z\| = 0.61, PASS |
| n=100,000: 12 pooled chain segments by lineage, PDG-inflated errors | \|z\| = 2.87 then 1.84 (consecutive), pooled ESS 136/122, PASS |
| n=100,000 local observables across all segments | perimeter/cell 1.1946–1.1965, cycles/cell 0.0725–0.0730, no lineage grouping |

Measured autocorrelation time, single-cell kernel: τ ≈ 0.3 · n^2.2 moves.
The mixed kernel (33% cut-and-paste branch moves, see below) collapses this:
τ(n=1000) drops from ~2×10⁶ to 2.9×10⁴ moves and τ(n=10⁴) from ~2×10⁸ to
2.2×10⁵ — a ~10³ speedup at the largest sizes.

## Cut-and-paste branch moves

Nonlocal move mixed at 1/3 probability with the single-cell move: pick a
uniform ordered adjacent pair (u, v); if that adjacency is a *bridge* of the
cell graph (alternating two-front BFS excluding the edge), detach the
component on v's side (≤ n/2 cells), apply a uniform lattice symmetry about
v, and re-attach it at a uniform (cell, direction) of the remainder,
requiring the placement to touch the remainder at exactly one contact.  A
bridge is the unique edge between the two sides, so the adjacency count is
conserved and forward/reverse selection probabilities are equal: the kernel
is symmetric with no Hastings factor.  Validated by the same chi-square
battery (see table) before production use.

(vibecoded — proceed at own risk. if it breaks you get to keep both pieces.)
