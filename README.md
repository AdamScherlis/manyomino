# manyomino

Uniformly sample extremely large polyominoes.

Samples polyominoes (lattice animals on Z², 4-connectivity) **uniformly at
random for fixed size n**, from n ~ 100 up to n ~ 10,000, renders them, and
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
  time series.
- `python/twoseed.py` — two-seed convergence test (bar vs rect chain).
- `python/nufit.py` — weighted log-log fit of ⟨R_g⟩ vs n.
- `python/render.py` — stdlib PNG renderer (`bw` and graph-distance `dist`
  modes).
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
| two-seed convergence, ν fit | in progress (production runs) |

Measured autocorrelation time: τ ≈ 0.3 · n^2.2 moves (n = 100..1000), so
n = 10⁴ needs ~2×10⁸ moves per independent sample.

(vibecoded — proceed at own risk. if it breaks you get to keep both pieces.)
