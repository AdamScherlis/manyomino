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
- asphericity ⟨(λ₁−λ₂)²/(λ₁+λ₂)²⟩ ≈ 0.35–0.42 (few-sample; the asph
  column in newer time series will sharpen this).

## 9. Performance summary

- mixed kernel, equilibrium: ~2.5×10⁵ moves/s (n=10³), ~5×10⁴ (n=10⁴),
  ~3×10⁴ (n=10⁵) on one core; cost ~ 40ns · 2√n per move.
- memory: grids scale with the shape's bounding box (~n^1.29); ~1GB/chain
  at n=3×10⁵, ~3.5GB at n=10⁶ after the local-bitmap refactor of the cold
  paths.
- the remaining asymptotic bottleneck is exact cut-cell testing
  (O(√n)/move); dynamic biconnectivity would remove it but is deliberately
  out of scope.
