# manyomino

Uniformly sample extremely large polyominoes.

(vibecoded — proceed at own risk. if it breaks you get to keep both pieces.)

A *fixed polyomino* is a finite edge-connected set of cells of the square grid, counted up to
translation only. `manyomino` draws one uniformly at random from all fixed polyominoes with
`n` cells, for `n` up to the millions, using a Markov chain with two kinds of move:

* **relocate** — move one cell to an empty square on the boundary, if the rest stays connected;
* **transplant** — cut off a whole branch at a cut cell, rotate/reflect it, and re-attach it
  at another cell, if it fits there touching nothing else.

Each move is exactly as likely as its reverse, so the uniform distribution is stationary, and
relocation alone can reach every polyomino. The writeup, with pictures and a live demo, is at
<https://adam.scherl.is/static/manyomino/explain.html>.

## Why two moves

Relocation by itself is correct but slow: its autocorrelation time (for the squared radius of
gyration) grows like `n^2.6` moves, because mass can only diffuse through the shape one cell
at a time. With transplants it is about `15 n` moves for the radius, and about `40 n` for the
slowest property a deliberate search could find. See [Results](#results) and
[Convergence checks](#convergence-checks).

## Build and run

```sh
g++ -O3 -march=native -std=c++17 -o manyomino src/manyomino.cpp

# 10,000 cells, 20M moves, a line of observables every 30,000 moves, final shape to a file
./manyomino -n 10000 -moves 2e7 -every 30000 -out ts.txt -snap sample.bin

# continue from a saved shape; or inflate one 2x (4x the cells) as a starting point
./manyomino -load sample.bin -moves 2e7 -snap sample2.bin
./manyomino -load sample.bin -inflate 2 -n 40000 -moves 1e8 -snap big.bin
```

Useful options:

| option | meaning |
|---|---|
| `-n N` | number of cells |
| `-moves M`, `-tsec S` | stop after M proposed moves, or S seconds |
| `-pcp P` | fraction of moves that are transplants (default 0.5; the production runs use 0.4) |
| `-relK 0` | use the plain relocate kernel with an exact global connectivity test |
| `-relA A` | relocate draws its radius with `P(K ≥ r) = r^-A` (default 1; every result below was produced with 2, see Loops) |
| `-shape rect\|bar` | starting shape (a near-square block, or a 1×N bar) |
| `-load F [-inflate k]` | start from a saved shape, optionally blown up k× |
| `-every E`, `-out F` | write observables every E moves |
| `-blobs F` | also write, at each of those times, how many bridgeless components have at least 4, 5, …, 3,445 cells |
| `-snap F`, `-snapEvery E`, `-ckpt F` | save shapes: at the end, periodically, or as a rolling checkpoint |
| `-p P` | tilt by `(1-P)^perimeter`: the law of a size-N site-percolation cluster at density P |
| `-mode exact` | histogram every shape visited and chi-square it against uniform (small N) |
| `-mode means` | high-precision means of local observables |
| `-revtest 1` | check reversibility of every accepted move as the chain runs |

Saved shapes are little-endian `int32`: the cell count, then `x y` per cell. The time series
columns are `move sec rg2 asph perim cycles noncut leaves w h holes deg1 deg2 deg3 deg4`.

To start a chain at a large size, climb a ladder: equilibrate a small polyomino, inflate it
2×, equilibrate, repeat (`scripts/ladder.sh`). Inflation only supplies a starting point with
roughly the right overall shape; each rung is then run for many autocorrelation times.

## The moves, precisely

**Relocate.** Pick a cell `c`, a cell `u` and a direction `d` uniformly, and a radius `K`
with `P(K ≥ r) = r⁻¹`. Let `s = u + d` and `Q = P \ c`. If `s` is empty, continue with
probability `1/k`, where `k` is the number of neighbours of `s` in `Q`. Require that the
neighbours of `c` in `Q` are connected to each other within the radius-`K` box around `c`,
and likewise for `s`. Then move `c` to `s`. Given `K`, every valid pair `(c, s)` has
probability `1/4n²` and the conditions are symmetric in `c` and `s`. The condition at `c`
implies `Q` is connected. Most draws of `K` are small, so the test is cheap, but arbitrarily
large `K` has positive probability, which is what makes the chain irreducible. `-relK 0` sets
`K = ∞` always: the textbook kernel. `-relA A` makes the law `r^-A`; the kernel is
symmetric for any law of `K`. Every result in this README was produced with `-relA 2`, the
default until the study of loops below; `RELA=2` in front of any of the scripts runs them that way.

**Transplant.** Pick a cell `c`, a direction `d`, a cell `c'` and one of the eight lattice
symmetries `g`, all uniformly. Let `B` be the component of `P \ c` containing `c + d`. Require
`2|B| ≤ n − 1` and `c' ∉ B`. Let `T(x) = c' + g(x − c)`. Require `T(B)` to be disjoint from
`A = P \ B` and adjacent to `A` only at `c'`. Replace `B` by `T(B)`. The reverse move is
`(c', g(d), c, g⁻¹)`, which has the same probability `1/32n²`; the adjacency condition is
what guarantees `T(B)` is again a full component of `P' \ c'`.

The implementation decides a transplant with breadth-first searches run in lockstep from the
neighbours of `c` and from `c'`, testing the image of each branch cell as it is found and
stopping as soon as the move is certain to fail. That changes the cost, not the outcome:
`ref/manyomino_v1.cpp` is the naive version and produces bit-identical trajectories.

The same idea for lattice *trees* is the cut-and-paste algorithm of Janse van Rensburg and
Madras (J. Phys. A 25, 1992).

## Validation

* `-mode exact`: for `n ≤ 9` every fixed polyomino is visited (9,910 at `n = 9`) with
  frequencies consistent with uniform, for each kernel and the mixture, and for the tilted
  ensemble against its exact weights.
* `analysis/exact_small.py` computes exact averages over all polyominoes by brute force;
  `-mode means` reproduces them.
* The optimized transplant and the reference implementation agree move for move.
* `-revtest 1` replays the reverse of every accepted move and checks it is accepted and
  restores the shape: no failures in millions of moves at sizes from 60 to 64,000.
* The `-mode exact` and `-revtest` checks were repeated when the default radius law changed
  to `-relA 1` (9 cells: chi-square 9,868 for 9,909 degrees of freedom with both moves, 9,903
  relocate alone; no reversal failures at 60, 4,000 and 64,000 cells).
* Independent ladders from different seeds agree at every size.
* Counting removable cells recovers Klarner's constant (below).
* Samplers with no Markov chain agree with it on 60 observables: see
  [Convergence checks](#convergence-checks).

## Results

Two independent ladders (seeds A and B) at every other size, a third (C) at the sizes between;
one core each, about 16 core-hours in total; the largest size was then run about ten times
longer on four chains (`scripts/extend_hero.sh`, 30 core-hours). `analysis/summarize.py` and `analysis/build_page.py`
regenerate everything below from the raw time series, which are not in the repository
(`scripts/launch_all.sh` reruns them).

| cells | chains | run length / tau | <Rg^2> | <Rg> | tau / n | perimeter / n | removable / n | holes / n | lambda estimate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 125 | 1 | 81,208 | 90.948 ± 0.17 | 9.4562 ± 0.0084 | 14.5 | 1.23006 | 0.30302 | 0.02094 | 4.0603 ± 0.0032 |
| 250 | 2 | 64,669 | 222.51 ± 0.37 | 14.788 ± 0.012 | 15.6 | 1.21278 | 0.29866 | 0.02199 | 4.0610 ± 0.0029 |
| 500 | 1 | 13,709 | 543.3 ± 3 | 23.11 ± 0.061 | 15.3 | 1.20384 | 0.29614 | 0.02245 | 4.0652 ± 0.0058 |
| 1,000 | 2 | 15,379 | 1323.8 ± 4.1 | 36.07 ± 0.056 | 14.5 | 1.19938 | 0.29541 | 0.02284 | 4.0601 ± 0.0027 |
| 2,000 | 1 | 6,219 | 3250 ± 14 | 56.502 ± 0.12 | 15.1 | 1.19714 | 0.29475 | 0.02293 | 4.0619 ± 0.0007 |
| 4,000 | 2 | 9,467 | 7853.6 ± 34 | 87.858 ± 0.19 | 14.4 | 1.19612 | 0.29460 | 0.02300 | 4.0601 ± 0.0022 |
| 8,000 | 1 | 2,384 | 19186 ± 2e+02 | 137.29 ± 0.7 | 13.7 | 1.19526 | 0.29433 | 0.02309 | 4.0610 ± 0.0026 |
| 16,000 | 2 | 2,081 | 46699 ± 4.4e+02 | 214.22 ± 0.99 | 17.6 | 1.19516 | 0.29428 | 0.02307 | 4.0613 ± 0.0022 |
| 32,000 | 1 | 634 | 1.1246e+05 ± 1.8e+03 | 332.32 ± 2.5 | 14.2 | 1.19492 | 0.29446 | 0.02309 | 4.0581 ± 0.0028 |
| 64,000 | 2 | 430 | 2.7559e+05 ± 4.6e+03 | 520.5 ± 4.3 | 21.8 | 1.19513 | 0.29417 | 0.02306 | 4.0628 ± 0.0023 |
| 128,000 | 1 | 225 | 6.9162e+05 ± 2.2e+04 | 822.09 ± 12 | 12.0 | 1.19520 | 0.29410 | 0.02304 | 4.0639 ± 0.0026 |
| 256,000 | 2 | 104 | 1.6449e+06 ± 5.3e+04 | 1273.7 ± 19 | 22.5 | 1.19511 | 0.29424 | 0.02305 | 4.0616 ± 0.0017 |
| 512,000 | 1 | 130 | 3.8113e+06 ± 1.3e+05 | 1935.5 ± 33 | 12.4 | 1.19512 | 0.29400 | 0.02308 | 4.0651 ± 0.0012 |
| 1,024,000 | 4 | 626 | 9.6899e+06 ± 1.8e+05 | 3087.9 ± 28 | 8.1 | 1.19503 | 0.29423 | 0.02308 | 4.0615 ± 0.0007 |

Each run starts from a 2x-inflated sample of the size below; the first 300 n moves (and at least a fifth of the run) are discarded. tau is the integrated autocorrelation time of Rg^2 in moves, and run length counts all moves. Above about 100,000 cells the runs are only tens of tau long after burn-in, so tau itself is poorly determined there (hence the scatter in tau / n).

Relocate-only autocorrelation time (moves), n = 50: 3.84e+03, 100: 2.41e+04, 200: 1.41e+05, 400: 9.33e+05, 800: 6.09e+06, 1600: 1.54e+07; fit tau ∝ n^2.62 (±0.04).
With transplants: tau ∝ n^0.98 (±0.02), about 15 moves per cell.
Size exponent from n ≥ 4,000: nu = 0.6416 ± 0.0013 (literature 0.6412 ± 0.0005).
Klarner's constant from n ≥ 1,000: lambda = 4.0619 ± 0.0004 (accepted 4.0625696; chi-square of the sizes about their mean 10.9 for 10 degrees of freedom, error scaled accordingly).

## Convergence checks

There is no proof that this chain mixes as fast as it appears to (only the weak bound described
at the end of this section), and the radius of gyration could in principle relax quickly while
something else lags. `scripts/convergence.sh` and `analysis/convergence.py`
attack that from three sides, using 60 observables per shape (`src/features.h`: moments of the
mass distribution to order 6, radial profile, local pattern counts, holes, the cut structure
of the shape, graph distances).

**1. References with no Markov chain.** Exact enumeration of all 7,204,874 / 104,592,937
polyominoes of 14 / 16 cells (`src/enumerate.cpp`), uniform and percolation-weighted. Exact
independent percolation clusters of 250, 500 and 1,000 cells (`src/leath.cpp`), which have
the law the chain samples with `-p`. Uniform polyominoes from grown, resampled populations
(`src/grow.cpp`: each polyomino with a marked cell has exactly one growth history, so
importance weights are exact; a look-ahead tilt `kappa^(untried sites)` keeps the weights
flat). A chain run at the wrong density is the positive control.

| reference (no Markov chain) | details | comparisons | rms z | largest | beyond 2 s.e. |
|---|---|---:|---:|---:|---:|
| all 7,204,874 shapes of 14 cells, uniform | 2,000,001 chain samples; precision 7.5e-04 sd | 60 | 1.03 | 2.33 | 2 |
| all 104,592,937 shapes of 16 cells, uniform | 2,000,001 chain samples; precision 7.7e-04 sd | 60 | 0.70 | 2.23 | 1 |
| all 7,204,874 shapes of 14 cells, weighted (1-p)^perimeter, p = 0.5 | 2,000,001 chain samples; precision 7.6e-04 sd | 60 | 1.09 | 2.24 | 3 |
| all 7,204,874 shapes of 14 cells, weighted (1-p)^perimeter, p = 0.592746 | 2,000,001 chain samples; precision 8.6e-04 sd | 60 | 0.87 | 1.84 | 0 |
| 906,300 exact percolation clusters of 60 cells at p = 0.5 | 1 chain; learned test z = -0.67; quartiles of each distribution: 4 of 176 beyond 2 s.e., largest 2.43 | 60 | 0.93 | 2.17 | 2 |
| 532,808 exact percolation clusters of 250 cells at p = 0.5 | 13 chains (chi2/dof between them 1.07); learned test z = +0.80; quartiles of each distribution: 6 of 176 beyond 2 s.e., largest 2.47 | 60 | 0.72 | 2.08 | 2 |
| 20,000 exact percolation clusters of 500 cells at p = 0.55 | 2 chains (chi2/dof between them 1.14); learned test z = -0.44; quartiles of each distribution: 1 of 176 beyond 2 s.e., largest 2.63 | 60 | 0.81 | 2.00 | 0 |
| 2,585 exact percolation clusters of 1,000 cells at p = 0.592746 | 2 chains (chi2/dof between them 1.03); learned test z = -0.64; quartiles of each distribution: 5 of 176 beyond 2 s.e., largest 2.71 | 60 | 0.92 | 1.82 | 0 |
| 8 grown populations of 200,000 at 250 cells, uniform | count check -0.1% ± 0.9%; 18,453 distinct ancestors at n/2; precision 0.004 sd | 60 | 0.97 | 2.96 | 2 |
| 8 grown populations of 200,000 at 500 cells, uniform | count check +0.4% ± 1.6%; 6,636 distinct ancestors at n/2; precision 0.008 sd | 60 | 0.92 | 3.07 | 2 |
| 8 grown populations of 199,000 at 1,000 cells, uniform (fails its own count check) | count check -11.3% ± 3.5%; 2,132 distinct ancestors at n/2; precision 0.021 sd | 60 | 1.28 | 3.64 | 6 |
| 6 grown populations of 100,000 at 2,000 cells, uniform (fails its own count check) | count check -35.9% ± 9.2%; 296 distinct ancestors at n/2; precision 0.066 sd | 60 | 2.35 | 9.71 | 15 |
| 17 observables that are exactly zero by symmetry, at 6 sizes from 250 to 256,000 cells | two chains per size | 102 | 1.00 | 3.89 | 4 |
| chains started from a straight bar and a solid square of 16,000 cells, after 2,000 moves per cell, against the ladder chains | learned test z = -1.66 | 60 | 0.75 | 2.40 | 2 |
| control: chain at p = 0.51 against exact clusters at p = 0.50, 250 cells | learned test z = +6.68 | 60 | 7.26 | 18.56 | 31 |

z is (chain − reference) / standard error for each observable in `src/features.h`; chain errors are batch means. Over the 12 tests with a trustworthy reference: 762 comparisons, 20 beyond 2 standard errors (2.6%; 4.6% expected if every one were independent and exact), 3 beyond 3, largest 3.89 (m14 at 64,000 cells, mirror symmetry, 250 to 256,000). Observables are correlated, so these are not 762 independent tests.

The growth sampler estimates the number of polyominoes as it goes, and that number is known. Its estimate is the reference's own health check:

| cells | population | populations | count estimate vs known | distinct ancestors at n/2 | rms t | largest |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 199,900 | 8 | -0.1% ± 0.9% | 18,453 | 0.97 | 2.96 |
| 500 | 199,900 | 8 | +0.4% ± 1.6% | 6,636 | 0.92 | 3.07 |
| 1,000 | 198,800 | 8 | -11.3% ± 3.5% | 2,132 | 1.28 | 3.64 |
| 2,000 | 99,800 | 6 | -35.9% ± 9.2% | 296 | 2.35 | 9.71 |
| 1,000 | 12,400 | 24 | -15.5% ± 8.2% | 136 | 1.45 | 4.30 |
| 1,000 | 49,700 | 12 | -3.1% ± 9.0% | 527 | 1.35 | 3.39 |
| 2,000 | 25,000 | 16 | +52.7% ± 71.6% | 70 | 1.49 | 3.67 |

Eight independent chains at 250 cells from one starting shape (4.3e+09 moves in all): chi-square of the chain means about their mean, per degree of freedom, averaged over 60 observables: 1.06 (largest 2.26; 1 if the batch-means errors are right).

The largest symmetry violation (m12 at 64,000 cells, -0.092 ± 0.024 of its standard deviation) was measured again on two fresh chains at that size: +0.031 ± 0.033.

Autocorrelation times in moves per cell (integrated, Sokal window), two independent chains per size:

| cells | samples (one per n moves) | Rg² | slowest single observable | slowest combination found (held out) | A-vs-B learned test z |
|---:|---:|---:|---:|---:|---:|
| 250 | 857,880 | 15.0 ± 0.3 | 17.5 ± 0.4 (trunk10) | 22.5 ± 0.8 | +1.09 |
| 1,000 | 349,060 | 14.8 ± 0.5 | 25.8 ± 1.1 (trunk03) | 30.9 ± 2.0 | +0.49 |
| 4,000 | 125,954 | 14.4 ± 0.8 | 35.4 ± 2.9 (trunk01) | 43.1 ± 5.5 | +1.03 |
| 16,000 | 76,476 | 14.4 ± 1.0 | 41.3 ± 4.7 (trunk01) | 47.9 ± 8.1 | +1.40 |
| 64,000 | 30,585 | 13.6 ± 1.4 | 34.9 ± 5.8 (trunk01) | 41.5 ± 10.7 | +2.09 |
| 256,000 | 1,890 | 21.8 ± 14.2 | 46.5 ± 56.1 (trunk01) | 21.3 ± 15.7 | -2.07 |

At 256,000 cells the runs are too short to time the slow observables (errors above 50%); that size is left out of the fits and the chart.

Fits from 1,000 cells up, tau/n ∝ n^x: Rg² x = -0.016 ± 0.020; slowest single observable x = +0.131 ± 0.048; slowest combination x = +0.128 ± 0.050.

Two honest footnotes. The 250-cell cluster comparison was first run with five chains, and one
observable (a fourth-order moment) stood 3.6 standard errors from the exact value; eight
further chains from fresh exact starts, with more than twice the data, did not reproduce it,
and the table pools all thirteen. The "A-vs-B learned test" column trains a small network on
the first halves of two independent chains and scores it on their second halves; it is positive
at five of six sizes (largest +2.1), which is mild evidence that batch-means errors are
slightly small for the slowest properties at the larger sizes. The direct test of that at 250
cells (eight chains) finds none.

**2. A search for slow modes.** The autocorrelation time of every observable, and of the
slowest linear combinations of the observables and their squares (time-lagged component
analysis, fitted on one chain and timed on an independent one). The slowest thing found is
the trunk of the shape (number of cells whose removal detaches at least 1%), about three
times slower than the radius. All of this is a lower bound on the true relaxation time.

**3. Other starts and symmetry.** Chains started from a straight bar and from a solid square
agree with the ladder chains after burn-in; the 17 observables that are odd under a lattice
symmetry average to zero.

**A proof, of a weak kind.** [`proof/`](proof/README.md) shows that the relaxation time of this
chain is at most n^(O(log n)): no property of the shape can take exponentially long to
forget. The bound is far above what is measured, and the argument has not been reviewed.

**Properties nobody measured.** Two more checks address slow properties outside the 60
observables. For shapes of up to 10 cells, `-mode trans` and `analysis/spectrum.py` measure the
full table of transition rates between shapes and take its relaxation time, which bounds every
property at once: the slowest is the shape's orientation, and the 60 observables capture 98.6%
of it at 10 cells. At large sizes, `scripts/fossil.sh` and `analysis/fossil.py` follow the
largest piece of a starting shape that is still intact (cells keep their identity as they
move), the slowest kind of property there is; it is gone after a few hundred moves per cell.
`analysis/hero_slow.py` times each observable directly at 1,024,000 cells. The trunk is no
slower there than at 64,000 cells, but one observable, the size of the largest clump of
cells held together by loops (`blob_max`, the largest bridgeless component), reads about three
times slower than at 64,000 cells, with an uncertainty of more than half. The next paragraph
is about why.

**Loops.** Large clumps are single thin rings: a clump of 200 cells is a loop one cell wide,
about 43 cells across, with unit squares attached. Transplant cannot open a loop, and relocate
lifts a cell `c` only if its neighbours are connected within the radius-`K` box, which for a
cell on a ring means a box that takes in the whole ring. With `P(K ≥ r) = r⁻²`, the rule all
the results here were produced with, that is rare for a large ring, so large rings open and close slowly (their number is right, since closing is
suppressed by the same factor). `scripts/clumps.sh` and `analysis/clumps.py` record, once per
move per cell, how many clumps of at least each size the shape holds (`-blobs`), at 4,000,
16,000 and 64,000 cells. Two things come out. The number of clumps of each size per cell is
the same at every shape size, falling about 12-fold per doubling of the clump, so a larger
shape holds larger clumps only because it has more room. And the autocorrelation time of
that number is also the same at every shape size but grows with the clump, roughly in
proportion to its size from 32 cells up: 5 moves per cell at 16 cells, 20 at 64, 40 at 128,
about 80 at 363. Putting 64 unrelated stretches of the 16,000-cell record side by side as a
stand-in for one shape of 1,024,000 cells reproduces the million-cell reading (largest clump
140 ± 8 cells remembered for 76 ± 36 moves per cell; measured 142 cells and 96 ± 56). So the
relaxation time of 40 to 46 moves per cell is that of properties a typical shape has plenty
of; above it sits a ladder of rare ring properties, taller the larger the shape. They carry
almost no weight in the averages reported (a 4,000-cell shape holds a clump of 128 cells or
more 0.3% of the time).

`-relA 1` draws large boxes more often (`P(K ≥ r) = 1/r`). The clump counts are unchanged,
the exhaustive test passes, and the ladder is gone: clumps of every size from 16 to 363
cells are forgotten in under 2 moves per cell, and at 1,024,000 cells the largest clump in
1.1 ± 0.3 (two chains of about 200 moves per cell). The gain at each size is close to the
width of the rings of that size (25-fold at 128 cells, where rings are about 30 cells
across; 58-fold at 256, about 50 across), as it should be: a box wide enough for a ring `d`
across is drawn `d` times as often. The trunk is not helped (40 ± 2 against 41 ± 4 at 16,000
cells; 62 ± 16 against 39 ± 4 at 64,000, on a run half as long). It cannot be hurt either:
the connectivity test passes in a larger box whenever it passes in a smaller one, so the
`-relA 1` kernel has every off-diagonal transition probability at least as large as that of
`-relA 2`, and by Peskun's theorem (Biometrika 60, 1973) the asymptotic variance of every
average is no larger per move, and the spectral gap no smaller. The cost is about 9% more
time per move at 4,000 and 16,000 cells, and less at larger sizes, where transplants take
most of the time (searched cells per move rise by 7% at 64,000 cells and about 5% at
1,024,000). `-relA 1` is the default from here on; every result in this README was produced
with `-relA 2`.

Clumps (`analysis/clumps.py`, `analysis/clump_report.py`): number per cell and autocorrelation time, in moves per cell, by clump size.

| clump of at least | per cell | lasts with `-relA 2`, 4,000 cells | lasts with `-relA 2`, 16,000 cells | lasts with `-relA 2`, 64,000 cells | lasts with `-relA 1` |
|---:|---:|---:|---:|---:|---:|
| 4 cells | 0.0486 | 4.4 ± 0.1 | 4.4 ± 0.0 | 4.4 ± 0.2 | 4.1 ± 0.0 |
| 8 cells | 0.00663 | 3.4 ± 0.0 | 3.5 ± 0.1 | 3.1 ± 0.1 | 2.6 ± 0.0 |
| 16 cells | 0.000575 | 5.0 ± 0.1 | 5.1 ± 0.0 | 4.9 ± 0.3 | 1.7 ± 0.0 |
| 32 cells | 4.86e-05 | 11.2 ± 0.2 | 11.1 ± 0.2 | 12.7 ± 0.8 | 1.6 ± 0.0 |
| 64 cells | 4.95e-06 | 24.0 ± 1.7 | 19.9 ± 0.8 | 19.9 ± 1.1 | 1.6 ± 0.0 |
| 128 cells | 5.48e-07 | 41.5 ± 3.6 | 37.6 ± 5.1 | 32.0 ± 10.4 | 1.6 ± 0.1 |
| 182 cells | 2.08e-07 | 52.7 ± 4.8 | 53.5 ± 9.0 |  | 1.5 ± 0.1 |
| 256 cells | 8.45e-08 | 67.1 ± 6.8 | 72.6 ± 27.1 |  | 1.2 ± 0.1 |
| 363 cells | 6.28e-08 | 79.1 ± 10.0 |  |  | 1.1 ± 0.2 |

Times are integrated autocorrelation times, in moves per cell, of the number of clumps of at least that size (2,381,877 samples at 4,000 cells, 755,268 samples at 16,000 cells, 50,564 samples at 64,000 cells, one per move per cell; errors by leaving out blocks).

| cells | `-relA 2`: largest clump, mean | its time | `-relA 1`: mean | time |
|---:|---:|---:|---:|---:|
| 4,000 | 25.5 ± 0.1 | 39.1 ± 5.8 | 25.3 ± 0.0 | 1.7 ± 0.0 |
| 16,000 | 38.4 ± 0.2 | 40.9 ± 12.7 | 38.4 ± 0.2 | 1.6 ± 0.1 |
| 64,000 | 56.8 ± 1.1 | 33.1 ± 9.0 | 58.8 ± 0.5 | 1.7 ± 0.1 |
| 1,024,000 |  |  | 134.3 ± 4.6 | 1.1 ± 0.3 |

The largest clump of K independent stretches of a smaller record, taken together as one shape of K times the size (time: median over 24 relative alignments of the stretches, with the larger of their spread and the block error):

| stands for | built from | largest clump, mean | its time | range over alignments |
|---:|---|---:|---:|---|
| 256,000 | 64 × 4,000 | 98 ± 3 | 62 ± 10 | 50 to 74 |
| 1,024,000 | 256 × 4,000 | 163 ± 11 | 65 ± 27 | 49 to 147 |
| 4,096,000 | 1024 × 4,000 | 272 ± 23 | 72 ± 16 | 36 to 95 |
| 1,024,000 | 64 × 16,000 | 140 ± 8 | 76 ± 36 | 57 to 164 |
| 4,096,000 | 256 × 16,000 | 221 ± 15 | 84 ± 36 | 48 to 127 |

Other observables with `-relA 1` against `-relA 2` (`conv/` chains), moves per cell:

| observable | 16,000: `-relA 2` | 16,000: `-relA 1` | 64,000: `-relA 2` | 64,000: `-relA 1` |
|---|---:|---:|---:|---:|
| trunk01 | 40.5 ± 3.8 | 40.0 ± 2.5 | 39.0 ± 4.3 | 61.6 ± 16.2 |
| trunk03 | 34.7 ± 5.2 | 35.2 ± 1.3 | 31.7 ± 2.4 | 48.1 ± 14.8 |
| trunk10 | 21.6 ± 1.8 | 22.9 ± 1.2 | 23.8 ± 2.5 | 24.6 ± 4.4 |
| perim | 20.3 ± 1.5 | 18.2 ± 0.7 | 18.3 ± 1.8 | 16.9 ± 2.0 |
| rg2 | 14.7 ± 0.9 | 13.7 ± 0.2 | 13.3 ± 0.9 | 15.1 ± 1.5 |
| blob_max | 27.7 ± 8.3 | 1.6 ± 0.1 | 30.3 ± 9.6 | 1.7 ± 0.1 |

Speed, million moves per second (`-relA 2`, `-relA 1`): 4,000 cells: 0.88, 0.80; 16,000 cells: 0.56, 0.51.

Relaxation time of the whole chain on small shapes, from the measured table of transition rates (`analysis/spectrum.py`), in moves per cell:

| cells | shapes | slowest property of all | slowest combination of the observables | share of the slowest property inside the library |
|---|---|---|---|---|
| 4 | 19 | 3.68 | 3.68 | 100.0% |
| 5 | 63 | 4.79 | 4.79 | 100.0% |
| 6 | 216 | 6.96 | 6.92 | 99.8% |
| 7 | 760 | 7.54 | 7.43 | 99.6% |
| 8 | 2,725 | 9.68 | 9.40 | 99.2% |
| 9 | 9,910 | 9.64 | 9.32 | 98.9% |
| 10 | 36,446 | 11.61 | 11.09 | 98.6% |

Largest piece of a starting shape still intact (`analysis/fossil.py`), median over runs, in cells:

| cells | runs | after 5 n moves | after 50 n moves | after 100 n moves | after 200 n moves | after 400 n moves | after 600 n moves | every cell has moved by |
|---|---|---|---|---|---|---|---|---|
| 4,000 | 7 | 2,821 | 235 | 27 | 9 | 3 | 3 | 20 to 190 n moves |
| 16,000 | 8 | 10,856 | 796 | 134 | 16 | 3 | 3 | 40 to 200 n moves |
| 64,000 | 4 | 42,434 | 1,990 | 343 | 35 | 4 | 3 | 40 to 190 n moves |
| 256,000 | 1 | 111,280 | 5,397 | 234 | 187 | 13 | 3 | 65 n moves |
| 1,024,000 | 16 pairs of snapshots 150 n moves apart | | | | | | | largest intact piece 181 to 3,218 cells |
| 1,024,000 | 12 pairs of snapshots 300 n moves apart | | | | | | | largest intact piece 28 to 162 cells |
| 1,024,000 | 8 pairs of snapshots 450 n moves apart | | | | | | | largest intact piece 6 to 12 cells |
| 1,024,000 | 4 pairs of snapshots 600 n moves apart | | | | | | | largest intact piece 3 to 4 cells |

Autocorrelation times at 1,024,000 cells from the extended runs (`analysis/hero_slow.py`; 4 chains, 3,885 moves per cell in all), in moves per cell, with the same estimator applied to equally short stretches of the smaller chains (mean ± spread over sets of stretches) and the long-run value there:

| observable | 1,024,000 cells | 16,000: short stretches | 16,000: long run | 64,000: short stretches | 64,000: long run |
|---|---|---|---|---|---|
| rg2 | 8.3 ± 0.3 | 15.4 ± 5.2 | 14.4 | 13.2 ± 3.0 | 13.6 |
| blob_max | 96.2 ± 55.8 | 19.4 ± 12.2 | 27.1 | 21.5 ± 6.2 | 30.5 |
| trunk01 | 32.0 ± 8.7 | 34.7 ± 14.5 | 41.3 | 39.9 ± 18.3 | 34.9 |
| trunk03 | 30.1 ± 7.8 | 34.6 ± 15.2 | 36.3 | 28.2 ± 8.2 | 28.7 |
| m15 | 19.9 ± 6.8 | 14.0 ± 3.4 | 13.7 | 13.6 ± 4.0 | 13.0 |
| trunk10 | 18.4 ± 4.5 | 22.6 ± 8.3 | 21.9 | 21.4 ± 5.0 | 19.8 |
| m13 | 17.8 ± 5.1 | 15.5 ± 3.7 | 15.3 | 14.2 ± 4.5 | 14.1 |
| perim | 17.5 ± 1.4 | 21.8 ± 10.0 | 20.8 | 16.5 ± 4.6 | 18.4 |
| detach | 15.5 ± 3.2 | 17.7 ± 9.3 | 16.5 | 15.8 ± 3.4 | 16.6 |

What this does not cover: beyond 1,000 cells there is no independent sampler (the growth
sampler's own count estimate fails from about 1,000 cells at these population sizes), and a
slow mode outside the span of the 60 observables would not be seen, as the ring
properties above were not until `blob_max` was timed at 1,024,000 cells.

## Music video

`video/` renders the chain as a two-minute music video. A 50×50 square is dissolved by the two
moves at 138 BPM: relocation proposals run on every frame, transplant proposals come in bursts
on the sixteenth-note grid, so branches jump on the beat. (Any fixed sequence of the two
kernels leaves the uniform distribution stationary, so it is still the same chain.) The
soundtrack is synthesised from the event log. Each accepted relocation is a tick. The largest
branch moved in each burst is a mallet note, lower and longer the bigger the branch, and in
the drops it also fires the lead synth: one stab per sixteenth in which a branch moved,
brighter, louder and higher up the chord the bigger the branch. Chords, bass and drums are a
fixed arrangement: a mellow opening that builds into two trance drops.

```sh
node video/dry.js events.json                 # the chain's event log (deterministic)
node video/make.js out                        # out/silent.mp4, via headless Chromium + ffmpeg
python3 video/synth.py events.json music.wav
ffmpeg -i out/silent.mp4 -i music.wav -c:v copy -c:a aac -shortest video.mp4
```

## Layout

```
src/manyomino.cpp      the sampler
src/features.h         60 shape observables (geometry, local counts, cut structure, graph distances)
src/enumerate.cpp      exact averages of them over every polyomino of n cells
src/leath.cpp          exact percolation clusters of a given size
src/grow.cpp           uniform polyominoes from a resampled population of growing clusters
ref/manyomino_v1.cpp   naive reference implementation
proof/                 a quasi-polynomial bound on the relaxation time, and its machine check
scripts/               production ladders (resumable); convergence runs
analysis/              autocorrelation times, data reduction, rendering, page build
page/                  source of the explainer page
video/                 music video generator
```
