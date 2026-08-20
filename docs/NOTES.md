# Methods notes — what was built beyond the handoff spec

This documents design decisions, correctness arguments, and measured results
accumulated while implementing `docs/HANDOFF.md`, in enough detail to
reconstruct or audit any of it.  Everything here is implemented in
`rust/src/main.rs` (+ `rust/src/perm.rs`) and the `python/` analysis tools.

## 1. Single-cell kernel (the handoff baseline)

As specified: pick a uniform cell c; if A−{c} disconnects, stay put; else
move c to a uniform perimeter *site* of A−{c} (c itself included, so s=c is a
legal no-op).  Symmetric because A−c = A′−s forces |P| equal both ways.

Implementation: grid-backed indexed sets for cells and perimeter sites
(dense array + position index, swap-remove); incremental perimeter
maintenance (remove c → c becomes a site, neighbors that only touched c drop
out; add s → s leaves, its empty neighbors join); O(1) R_g², gyration
tensor, perimeter length and cycle count (E − n + 1) via incremental sums.

Connectivity test, in escalating order of cost:
1. leaf shortcut (one occupied neighbor);
2. 3×3 ring-arc test: consecutive ring cells are 4-adjacent
   (orthogonal/diagonal alternation), so a single occupied arc containing
   all orthogonal neighbors proves removability; also the one-arc-with-all-
   orthogonals case with a separate diagonal arc;
3. alternating multi-front BFS **in cell-index space** over maintained
   adjacency links (`nbr: Vec<[u32;4]>`), one front per neighbor arc, tiny
   union-find; per-front visited bitplanes over indices (~n/8 bytes, L1-
   resident).  Proving "cut" costs O(min component) ≈ O(√n) at equilibrium
   (subtree sizes ~ k^{-3/2}).

At equilibrium ~70% of proposals are cut cells, so (3) dominates the move
cost: ~40 ns/node · ~2√n nodes.  Memory-layout experiments (grid bitmaps,
u8 tags) did not beat this; it is ~5 dependent accesses per node.

## 2. Cut-and-paste branch moves (mixed kernel)

Move (mixed at probability 1/3): pick a uniform *ordered adjacent pair*
(u,v) via uniform cell + uniform direction (empty → no-op).  If the edge
(u,v) is a **bridge** of the cell-adjacency graph, detach the component B on
v's side (reject if |B| > n/2), apply a uniform D₄ element about v, and
re-attach at a uniform (cell u′, direction d) of the remainder R, requiring
the placement to be disjoint from R with **exactly one** R-contact.

No Hastings factor: a bridge is the unique edge between the sides, so
E(A′) = E(A); forward and reverse tuple probabilities are both
(1/4n)·(1/8)·(1/4|R|), and the tuple↔reverse-tuple map is a bijection.
The empty-pair/not-bridge/too-big/invalid-placement outcomes are rejections
(self-loops) and need no balancing.

Bridge test: two-front alternating BFS excluding the tested edge; fronts
meeting ⇒ cycle edge; a side exhausting ⇒ bridge, B enumerated from its
queue.  **Node cap** (`--cp-cap`, default 4000): give up (reject) once both
fronts exceed K pops.  Symmetric because on any acceptable move the forward
and reverse tests have identical component sizes and front roles
(front 0 = R-side both ways), so resolution-within-cap depends only on
min(|B|,|R|).  Validated adversarially: chi-square at n=8 with K=2 (cap
firing 1.5×10⁶ times) stays uniform, p = 0.885.

Effect: τ(n=1000) drops ~70× (2.9×10⁴ vs ~2×10⁶ moves), τ(n=10⁴) ~10³×
(2.2×10⁵ vs ~2×10⁸); scaling roughly τ ~ n^1.1 vs n^2.2.

## 3. Seeds for the two-seed test at large n

- 1×n bar: correct but pathological — every interior cell is a cut cell
  costing O(n) to prove.
- 2×(n/2) bar (`bar2`): no interior cut cells *initially*, but as soon as
  opposing nicks appear, cells under nicks are cut cells with O(n) proofs.
- **width-8 bar (`bar8`)**: R_g still ~5× equilibrium at n=10⁵, but random
  nicks cannot sever it, so tests stay local while it crumples.  7×
  throughput vs bar2 at n=10⁵.
- Compact seeds (√n×√n `rect`) at n ≳ 3×10⁵ crawl through a *coarsening*
  transient (dumbbell necks are cut cells with O(n) proofs).  Fix:
  **inflation** (`python/inflate.py`): map each cell of an equilibrated
  smaller sample to a 2×2 block.  Large-scale geometry is already
  equilibrium-like; only sub-inflation-scale (fast) modes must relax.
  Two-seed honesty is preserved by inflating *independent lineages*
  (bar-descended vs rect-descended chains).

## 4. Error methodology

Sokal-windowed integrated autocorrelation time (window W ≥ 5τ(W)) on the
binned R_g² series is biased **low** when slow modes extend past the window
— observed concretely at n=5000, where four chains were each internally
stable (half-vs-half z < 1) yet scattered more than their nominal errors.
`analyze.py` therefore reports the more conservative of the Sokal stderr and
coarse (8-batch) batch-means stderr, scaling τ/ESS consistently.  Two-seed
(and four-chain, `fourseed.py`) agreement between independently seeded
chains is the publication gate, not any single-chain diagnostic.

## 5. Validation ladder (all pass; see README table for numbers)

exact counts n≤10 → chi-square vs exact shape lists at n=5/6/8 (reference,
fast, mixed kernel) → move-statistics agreement reference-vs-fast →
cycle-count and diagonal-near-miss distributions vs the exact n=12 census
(505,861 shapes; thinning ≫ τ) → two-seed at every production size →
ν, λ, θ against literature → detailed balance audits (below).

## 6. Detailed-balance audits on real big-cycle transitions (n=10⁴)

`cyclewatch` captures the exact transitions that seal and later breach a
≥100-cell enclosure (strict seal: the region must flood to the outside in
the predecessor).  For both transitions:
- independent from-scratch verification (`verify_db.py`): A−c = A′−s,
  both directions legal, |P| identical, T = 1/(n|P|) exactly equal;
- implementation audit (`dbtest`): exhaustive removability boundary over
  all 10⁴ cells vs full BFS (0 mismatches), incremental perimeter equals
  from-scratch set exactly, destination distribution uniform over all ~12k
  sites (p = 0.11–0.83), forward/reverse destination distributions
  statistically identical;
- end-to-end (`dbtest2`): 1.2×10¹² unmodified one-step trials per
  direction; realized P(A→A′) and P(A′→A) equal to each other and to
  1/(n|P|) within Poisson errors.

## 6b. A rare perimeter-corruption bug, found and fixed via the time series

The only correctness bug found after production began, worth recording for
its detection story.  Scanning every recorded time series for out-of-band
perimeter values turned up exactly two single-record excursions to
perim/n ≈ 0.68 (n=1000 b=0.3 step 6.978M; n=3000 b=0 step 26.366M), each
with a simultaneous R_g² spike and an asphericity ≈ 0.0155 record
immediately after — too similar to be fluctuations.  Because every run
logs its seed and the chain is deterministic, the exact states could be
replayed: `--dump-final` at the anomalous step showed the *actual* shape
was a normal branched polymer whose true perimeter (1188 sites at n=1000)
disagreed with the maintained set (685) — state corruption, not physics.
A per-step verification mode (`--percheck-from`) then pinpointed the
corrupting move: a **rejected cut-and-paste whose mid-move grid rebuild
returned shift = (0,0)**.  The code inferred "a rebuild happened" from
shift ≠ (0,0), but a rebuild can legitimately land on a zero shift (new
margin equal to the old bbox offset in both axes, needs a long branch
poking the border while the remainder sits deep inside — hence the R_g²
spike); the reject path then skipped the perimeter repair around the
restored branch, leaving per = perimeter of the remainder only, ~45% of
sites missing.  The periodic recentering rebuild heals it ~6.5k steps
later, which is why the excursions were single records.  Fixed by
tracking rebuilds with an explicit flag at all four (move × accept/reject)
sites; verified by full-length replays of both affected runs (25M and
36M steps) with invariant checks every 2000 steps (perimeter set,
adjacency links, edge count, sums, connectivity): zero failures.
Event rate before the fix:
~1 per 25M steps at n=1000–3000 (needs mid-move rebuild ∧ zero shift ∧
reject), a ~2×10⁻⁴ corrupted-sample fraction — no visible effect on any
published average, but the proposal distribution was briefly wrong in
those windows.

Moral: record cheap redundant observables (perimeter, cycles) in every
production run and band-scan them; log seeds so any anomaly is exactly
replayable.

## 7. PERM (independent non-MCMC cross-check)

Leath-style growth with deterministic FIFO discovery order: each candidate
site is occupied with probability p or permanently blocked.  Given the FIFO
rule, a positioned animal containing the seed has exactly **one** decision
path, so P(A) = p^{n−1}(1−p)^{b_A} exactly, and W = 1/P is an exact
importance weight (no growth-order multiplicity to count).  Log-space
weights; adaptive pruning/enrichment (thresholds from running Ẑ_k);
E[ΣW at size k] = k·a_k gives the animal counts.

Tuning: the percolation-critical p ≈ 0.593 produces compact clusters and
weight degeneracy (weight-ESS ~1.5% at n=300, upward-biased ⟨R_g²⟩ with
deceptive error bars).  Matching the animal ensemble's perimeter/site ratio
(measured 1.195/cell → p* = 1/(1+1.195) ≈ 0.456) raises weight-ESS ~20×;
p ≈ 0.40–0.45 used in production.  Always report weight-ESS.

Results: a₁..a₁₀ match A001168 (~0.5% statistical); λ = 4.058 and θ = 0.93
from a_k ~ Cλᵏ/k^θ over k=30..180 (Klarner λ = 4.0626…, theory θ = 1);
⟨R_g²⟩ agrees with MCMC at n = 10 (exact), 100, 200, 300, 1000.

## 8. Measured constants (uniform fixed polyominoes, 2D)

- ν = 0.6436 ± 0.0010 (stat) ± ~0.001 (sys, corrections-to-scaling scan),
  n = 100..30000; literature ≈ 0.6408.
- amplitude: ⟨R_g⟩ ≈ 0.422·n^0.6436.
- λ = 4.058 (PERM; exact 4.0626…), θ ≈ 0.93 (theory 1).
- cycle density → 0.073 per cell; perimeter density → 1.195 per cell.
- asphericity ⟨(λ₁−λ₂)²/(λ₁+λ₂)²⟩ = 0.392 ± 0.022 (n=10⁵ equilibrium
  time series, 8 segments).

## 9. Tilted ensemble: pi(A) ~ exp(-beta R_g^2)

A Metropolis factor min(1, exp(-beta*dRg2)) multiplies every (symmetric)
kernel; dRg2 is O(1) from the maintained coordinate sums, and rejected
moves restore state through the same production primitives used by
acceptance (validated by chi-square against the exact reweighted
distribution at small n, `python/validate_beta.py`).  Natural coupling:
**b = beta * <Rg^2>_0(n)** (penalty in units of the typical fluctuation
scale at that size).  Sweeps at n=1000 (b=-3..1000) and n=3000
(b=-2..3000), plus renders at n=100,000; figure `results/betafig.png`,
tables via `python/beta_analysis.py`.

Three regimes:

1. **Linear response, |b| <~ 1.**  d<Rg2>/db = -Var(Rg2)/<Rg2>_0 =
   -relVar * <Rg2>_0; measured relative slope ~= -(0.27)^2 per unit b at
   both sizes (relative sigma of Rg2 is ~0.27, n-independent).  Local
   observables (perimeter/cell, cycles/cell, movable fraction) do not
   move at all at this order; the entropy cost is O(1) nats *total*
   (dS/n ~ 1/n): the tilt just re-weights the global shape mode.

2. **Universal squeeze, 1 << b << 0.18n.**  <Rg2>/<Rg2>_0 collapses onto
   a single n-independent curve ~ b^-0.2 (n=1000 and n=3000 curves lie
   on top of each other; see figure).  Asphericity falls from 0.39
   toward 0 (shapes round out before they densify).  Locals still
   nearly frozen: at b=100 the perimeter density has moved only 1.19 ->
   1.18.  Entropy cost grows ~ b^0.85 but remains subextensive until...

3. **Collapse crossover at b_c ~ 0.18 n** (i.e. beta ~ 0.18n/<Rg2>_0,
   where the penalty starts competing with the *per-cell* free energy).
   Locals finally move: at n=1000, b=1000: perimeter/cell 1.199 ->
   1.058, cycles/cell 0.073 -> 0.099 (+37%), movable fraction 0.295 ->
   0.364, asphericity -> 0.008 (near-disk).  Entropy cost becomes
   extensive: dS/n ~= -0.076 nats/cell at (n=1000, b=1000), -0.081 at
   (n=3000, b=3000) — consistent with a common collapse branch as a
   function of b/n.  The animal is turning into a dense droplet: fewer
   perimeter sites, more cycles, and (counterintuitively) *more* movable
   cells, because compact bulk means fewer cut vertices than the
   branched-polymer phase.

   Negative beta (reward large Rg2) at b=-3 stretches Rg2 by 1.33x and
   raises asphericity 0.37 -> 0.55 (elongation), locals again frozen.
   The negative branch has no thermodynamic limit (logZ diverges as
   stretched configurations proliferate), so only small |b| is
   meaningful there.

Entropy via thermodynamic integration: d logZ/d beta = -<Rg^2>_beta,
integrated over the sweep grid (trapezoid), S(beta) - S(0) = logZ +
beta*<Rg2>_beta.  Since the beta=0 entropy per cell is log(lambda) ~=
1.40 nats, the collapse at b ~ 0.18n has burned ~5% of the animal's
entropy at the deepest points sampled.

n=100,000 renders (four chains, 120M steps each from an equilibrated
beta=0 snapshot; equilibrated stats over the last 60M):

|    b   |   beta    | Rg2/Rg2_0 | perim/n | cyc/n  |  asph | movable |
|--------|-----------|-----------|---------|--------|-------|---------|
|   -3   | -6.3e-6   | 1.41±0.09 | 1.1954  | 0.0726 | 0.545 | 0.2950  |
|   30   |  6.3e-5   | 0.551     | 1.1948  | 0.0730 | 0.083 | 0.2958  |
|  1000  |  2.1e-3   | 0.273     | 1.1934  | 0.0730 | 0.006 | 0.2941  |
| 30000  |  6.3e-2   | 0.140     | 1.1785  | 0.0754 | 0.003 | 0.3031  |

b=30 reproduces the n=1000 ratio (0.5455) at 100x the size — the squeeze
curve collapses in b alone.  The locals stay frozen at b=1000 here even
though the same b collapsed them at n=1000: local densification is
controlled by b/n (b_c ~ 0.18n = 18,000 at this size), cleanly separating
the two scaling variables.  At b=30000 (b/b_c = 1.7) the locals have just
begun to move.  Renders: gallery/tilt_b{m3,30,1000,30000}.png.

beta=0 distributions (fixed n, `python/betadist.py`, runs of 700-2000
tau at n=1000/3000/10000): two clean classes.
- Global shape modes, n-independent distributions: Rg2 rel-sigma
  0.27-0.28 with skewness +1.1-1.3 (broad right tail — stretched shapes
  are entropically favored); asphericity mean 0.37, rel-sigma 0.58-0.60.
  These never self-average.
- Local densities, self-averaging ~ n^{-1/2}: perim/cell rel-sigma
  0.016 / 0.010 / 0.005 and cycles/cell 0.124 / 0.071 / 0.038 at
  n = 1000 / 3000 / 10000; skewness ~ 0.  (The movable-fraction column
  is measured by a k=128-cell probe, so its recorded spread ~0.041 is
  probe-binomial noise; the true fluctuation is <~ 0.015.)
- Means: perim/cell 1.196-1.201, cycles/cell 0.0722-0.0727, movable
  0.293-0.295, all sizes — consistent with the section-8 constants.
- Caveat: the n=3000 series contains the one known corrupted perimeter
  record (section 6b; recorded 2039, true 3633 at step 26.366M), which
  drags its perim skewness to -5.8; all other moments unaffected.

## 10. Performance summary

- mixed kernel, equilibrium: ~2.5×10⁵ moves/s (n=10³), ~5×10⁴ (n=10⁴),
  ~3×10⁴ (n=10⁵) on one core; cost ~ 40ns · 2√n per move.
- memory: grids scale with the shape's bounding box (~n^1.29); ~1GB/chain
  at n=3×10⁵, ~3.5GB at n=10⁶ after the local-bitmap refactor of the cold
  paths.
- the remaining asymptotic bottleneck is exact cut-cell testing
  (O(√n)/move); dynamic biconnectivity would remove it but is deliberately
  out of scope.
