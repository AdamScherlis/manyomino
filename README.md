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
| `-freeze x0 y0 x1 y1` | cells in the box x0 ≤ x < x1, y0 ≤ y < y1 are never moved and none is put there: samples the rest of the shape given the box (with `-load`) |
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
more 0.3% of the time). Slow is not the same as wrong: the 26 saved shapes of 1,024,000 cells
hold 1,295,251 clumps of 4 or more cells against 1,293,551 expected from the long records,
and agree with them at every size up to 256 cells within counting error, except for 159
clumps of 64 or more against 132 expected, 2.4 counting errors high
(`analysis/clump_census.py`).

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
1,024,000). The two rules must give the same averages, and do (`analysis/rela_means.py`): of
the 60 observables, none differs by more than two standard errors at 16,000 cells (148,000
readings against 69,000) and 2 do at 64,000 (43,000 against 43,000; the largest 2.5). `-relA 1` is the default from here on; every result in this README was produced
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

## Endless zoom

`src/zoom.cpp` refines a polyomino to twice the resolution: every cell becomes a 2×2 block, and
the result is relaxed for a set number of sweeps by a walk that adds and removes single cells,
with the number of cells free and weight `z = 1/λ` per cell. The reason for that weight is
the usual one for a small part of a large system. In a uniformly random polyomino of `N`
cells, a region can gain a cell only if the rest of the shape gives one up, and the rest has
`A(N−1)/A(N) → 1/λ` times as many arrangements with one cell fewer. So, given its
surroundings, the fillings of a region that keep the whole connected are weighted by
`λ^(−cells inside)`. That is a heuristic, not a theorem; the tests below are the evidence.
The walk is reversible with respect to `z^n` on connected shapes (`-mode exact` checks this
over every shape of up to 8 cells), so on a doubled shape it rebuilds the fine detail and,
run briefly, leaves larger structure alone.
Connectivity is certified inside a box of random radius, as in relocate (`P(K ≥ r) = 1/r`,
capped at `-kmax`, default 256), so every move is local and the walk can be run on a window
(`-win`) with the surroundings held fixed. The cap matters: a loop wider than it is never
opened, doubling widens every loop, and with a cap of 64 four refinements in a row left the
largest clump a third too large on average.

`page/zoom.js` is the same walk in JavaScript on a 2,048×2,048 window around one point
(1,024 on a narrow screen), run in a worker, which is what the page's endless zoom runs. Each
level is the middle half of the previous one, doubled and relaxed, so the cost per level does
not depend on depth (about half a second). Three things differ from `src/zoom.cpp`, for speed
and for a picture without seams. It runs 16 sweeps, not 32 or more. Its certificate box has
radius `K ≥ r` with probability `1/(2r)`, capped at 48, which halves the searching. And the
sweeps are not spread evenly: the full number out to 72% of the way from the middle of the
window to its edge, then half as many on each of 12 successive square rings, and none beyond
94%. (This is done in stages, on nested squares with the outside held fixed, so each stage is
the same reversible walk on a smaller region.) At its edge a level is therefore exactly the
coarser level drawn with four cells for one, and the next level is always centred inside the
fully relaxed part.

The page draws a nest of these levels in one WebGL2 shader: each point on screen comes from
the finest level whose window contains it, with area-averaged coverage (mipmaps) since cells
are 0.1 to 0.65 pixels wide. No level is ever exchanged for another at a visible size; a new
one enters at the centre at a tenth of a pixel per cell. A level with cells of width `p` is
brightened by `p^-(2 − 1/ν)`, which is the factor by which finer levels have less ink per
unit area, and by the measured cell count where it has had fewer sweeps, so the brightness is
even across the nest and steady in time. Stepping the zoom by a fiftieth of a level at a time
across a level boundary, the frame-to-frame change stays between 1.4 and 2.3 grey levels out
of 255 throughout, with no step at any point.

**Outward.** The page also zooms out from single cells with no largest view
(`page/outward.js`, with `src/outward.c` compiled to `page/outward.wasm` by
`scripts/make_outward.sh`). A local walk cannot supply surroundings that are not there, but
the sampler can, exactly, at one scale: hold some cells of a polyomino fixed and run the
sampler on the rest (`-freeze x0 y0 x1 y1` in `manyomino`; the same moves in `outward.c`).
Nothing in the box moves and nothing is put there, the moves are still symmetric, so the
walk is uniform over the polyominoes of that size that agree with the box: a fair draw of
the surroundings of the box, given the box.

To go on without end the view is handed from one such polyomino to another twice as coarse.
A *world* is a polyomino of 16,000 cells standing for the shape around the starting cell
(the root, at the origin) at one scale. The next world is made in two steps.

1. *Halve.* Every 2×2 block holding a cell becomes a cell. That is too thick (2.1 cells of
   the finer world for each of the coarser, where a polyomino has 2.95), so a local walk
   thins it: the walk of the zoom in, with a weight of 0.72/λ per cell and two rules that
   keep the result lined up with the finer world. No cell is added further than one site
   from where the halved shape had cells, and none is taken away if that would leave one of
   those sites with no cell on it or round it. Sixteen sweeps; the result no longer changes
   after that. The weight is the one at which halving divides the number of cells by 2.95,
   and then the halved shape also has, layer by layer round the root, the number of cells a
   real shape has round one of its cells (`analysis/world_test.py`, below). Run freely for
   256 sweeps the same walk gives a typical polyomino too (`analysis/coarsen_test.py`), but
   it moves branch tips by many sites, and nothing drawn on that could be joined to the
   finer picture.
2. *Extend.* The box of half-width 8 about the root is held, and with it the two layers of
   the halved world just outside the box. Everything else is drawn afresh: the world is
   padded back to 16,000 cells and the sampler runs for 100 moves per cell. The padding is
   turned copies of the halved world hung on its own far tips. Padding with cells at random
   sites next to the shape gives a blob that the sampler needs 400 to 800 moves per cell to
   undo; the copies are forgotten in about 50 (`analysis/mix_test.py`).

The two layers are held because the box alone does not say where the finer world crosses
its edge, and the finer world knows: it has cells out there already. Drawn afresh, the
layer outside the box has about the right number of cells, but not where the finer picture's
strands end, and every crossing that does not line up has to be bridged in the picture.

Two strands that pass within a cell of each other without meeting touch after halving. If
they were taken for joined, the sampler would not join them anywhere else, and the picture,
which shows the finer world inside the box, would fall into pieces (with this left out, the
largest connected piece of a generated picture held 40% of its cells where a real one holds
84%). So each cell of the box carries the label of the part of the finer world it stands
for (a part is a connected piece of what lies within twice the box), and touching is not
taken for joined in two cases. Two cells of the box with different labels are not joined.
And a cell of the box is joined to a cell outside it only where the finer world really
crosses the edge of the box between the two; those sites are marked. Both hold for the
thinning walk (which also keeps each part in one piece inside the box, and never takes its
last cell) and for the sampler, which therefore has to join the parts outside the box, as
the real shape must. A part of three cells or fewer is merged with a part it touches.
`outward_native cut` checks the sampler with a pair of labelled cells against a count of
the shapes by exhaustive search (`analysis/held_count.py`), as `outward_native exact` and
`held` do without labels and with a held box.

The picture takes from each world only a ring, from its box out to twice the box; inside is
the finer level's picture and outside the next one's. A world's cells are large on the
screen, so the ring is refined six times (five on a small screen) by doubling and the local
walk, as for the zoom in, but never below real cells (finer than that, a cell just becomes
four). The walk runs for 32 sweeps at the first stage, 16 at the later ones and 64 at the
last, which is the one that sets the fine texture: with 16 or 32 there a classifier can tell
the texture from the real thing (below), with 64 it cannot.

What makes this hard is that a region refined on its own comes out thin along its edges.
Cells that should be joined to the rest through cells beyond the edge have nothing there,
and at each stage the last rows end up short: after three stages the last half layer of a
ring had 76% of the cells it should (counted as `analysis/stage_test.py` counts). An earlier version
refined each ring as a closed region and joined it to the picture inside at one held cell
for each place where the coarse world crosses the seam. That was right for the first two
seams and wrong by more the deeper the seam, because a real seam seen at the resolution of
real cells has several times as many joins as the coarse world has crossings, and because
both sides had been thinned once more at every stage. Run through five stages it left the
three rows either side of the fifth seam 14 to 18% short of cells, with 19% too few joins
across, while the ring as a whole was 5% over (the walk's weight had been raised to make up
the total). Now nothing that is kept is ever refined next to an edge:

- **The strip.** The ring is refined together with the two layers of its own world beyond
  it. Its outer rows therefore have their own neighbours beside them at every stage. The
  strip is not part of the picture, since the next world's ring belongs there; where the
  world leaves the strip one cell on its outer edge is held at every stage (an *anchor*),
  so that a strand joined to the rest only beyond the strip cannot draw back out of it.
- **The band.** Those two layers are the ones the sampler held when it drew the next world,
  so the next ring's first layers come from the same cells. At every stage the next ring's
  first two rows outside the seam are not its own: they are the first two rows of this
  level's strip at that stage, the very cells this picture was refined against. Its own
  cells start behind them. So the picture inside always meets exactly what it met when it
  was made, and the place where two separately refined things meet is two rows out in
  cells of that stage: half as far out, in real cells, at each finer stage. What lies
  within a given distance of the seam has therefore been refined on its own only since the
  stage at which that distance was a couple of cells.
- **Joins.** Behind the band the ring's own strands do not all meet it. What is joined to
  what is put back as the world has it. A piece that holds an anchor is joined to the rest
  beyond the strip, and so are all pieces holding anchors of one piece of the world (within
  the strip), which are joined to each other; a piece with no anchor is joined to whatever
  is nearest. A join is the shortest way through empty sites outside the seam, smallest
  pieces first; most are one or two cells. A small piece of the ring's own cells that
  would need a join longer than itself is dropped instead.
- **Joins add nothing.** For each cell a join adds, a cell is taken from a loose end in the
  same row or the next (a cell with one neighbour, or the spare corner of a block of four).
  The walk thins what it is given by a factor and to no target, so cells added at the seam
  stay there and are multiplied at every later stage: before this the rows outside the
  third seam ended 10% over. Taking them from the same rows matters too; taken from
  anywhere near, the first dozen rows outside the fifth seam were still 8% over.
- **The inside is held** until the last stage, every cell of it carrying the number of its
  piece in the finer level's picture at that stage. Two pieces of the inside that touch
  are not joined (they were not, or they would be one piece), so the ring cannot come to
  rely on a touch that the next stage, which puts the inside in afresh, would undo. The
  ring is joined to whatever it touches.
- **The last stage** has no labels, and a collar of the inside is free as well, 12 cells
  deep or a quarter of the way in where that is less: the picture is the one that stays,
  so whatever touches is joined, and the walk by itself never parts what is joined.
- **Weights.** A cell has weight 0.985/λ before the last stage and 0.9985/λ at it, not
  1/λ: at 1/λ a stage leaves about 1% too many cells. Both are set by the count of cells
  in rings of real shapes refined this way, and with no edges to make up for they are below
  one where the earlier version's was above.

On the page a level is drawn out to where the next level's reworked strip begins (three
quarters of the way to its edge for the first two levels, so 12 of the 16 real cells of the
first; all but 1/43 from the sixth level on), and the next level is drawn beyond.

`analysis/e2e_test.py` sets the whole of this against real samples. The real cells near a
cell of a 1,024,000-cell sample are given to the page's own code (run under Node by
`analysis/out_run.js`), six levels are made, and the first five, at the resolution of real
cells a picture 1,024 cells wide of which only the middle 32×32 was given, are compared with
the real cells around the same root. (The sixth is made so that the fifth ring has been
joined to a level outside it, as every ring on the page is.) The same script can leave the
sampler out: given the real cells out to 1,024 from the root it only halves and refines,
which puts the refining alone on trial, pairs every made picture with the real one it
should be, and shows a fault at a seam several times more sharply than the full run does.
`analysis/tell_test.py` then trains classifiers to tell the made pictures from the real
ones.

<!-- zoom-out table -->
390 roots in 26 samples, 5 levels. *real* and *made* are means over roots (per-cell rows: totals over roots divided), *s.e.* the difference in standard errors (paired by root; by resampling roots for the per-cell rows). Pieces are counted in the middle three quarters of the picture, so that its edge runs through the middle of the last ring.

| | real | made | difference | s.e. |
|---|---:|---:|---:|---:|
| cells 16-19 | 75.77 | 76.12 | +0.5% | +0.4 |
| cells 19-23 | 62.14 | 63.93 | +2.9% | +1.8 |
| cells 23-27 | 89.90 | 89.23 | -0.7% | -0.4 |
| cells 27-32 | 125 | 123 | -1.5% | -0.7 |
| cells 32-38 | 196 | 191 | -2.9% | -1.3 |
| cells 38-45 | 220 | 210 | -4.3% | -1.7 |
| cells 45-54 | 277 | 259 | -6.5% | -2.5 |
| cells 54-64 | 376 | 352 | -6.2% | -2.4 |
| cells 64-76 | 544 | 525 | -3.3% | -1.3 |
| cells 76-91 | 635 | 598 | -5.8% | -2.3 |
| cells 91-108 | 832 | 777 | -6.5% | -2.4 |
| cells 108-128 | 1,086 | 1,020 | -6.0% | -2.2 |
| cells 128-152 | 1,486 | 1,471 | -1.0% | -0.4 |
| cells 152-181 | 1,878 | 1,835 | -2.3% | -0.8 |
| cells 181-215 | 2,424 | 2,393 | -1.3% | -0.4 |
| cells 215-256 | 3,176 | 3,017 | -5.0% | -1.9 |
| cells 256-304 | 4,235 | 4,238 | +0.1% | +0.0 |
| cells 304-362 | 5,522 | 5,312 | -3.8% | -1.4 |
| cells 362-431 | 7,143 | 6,784 | -5.0% | -1.7 |
| cells 431-512 | 9,401 | 8,878 | -5.6% | -1.9 |
| L1 cells | 352 | 352 | -0.1% | -0.1 |
| L1 with 1 neighbours (per cell) | 0.1705 | 0.1695 | -0.6% | -0.8 |
| L1 with 2 neighbours (per cell) | 0.5427 | 0.5441 | +0.3% | +0.7 |
| L1 with 3 neighbours (per cell) | 0.2585 | 0.2580 | -0.2% | -0.4 |
| L1 with 4 neighbours (per cell) | 0.0283 | 0.0284 | +0.4% | +0.2 |
| L1 empty sites next to the shape (per cell) | 1.20 | 1.20 | +0.0% | +0.2 |
| L1 full 2x2 (per cell) | 0.0698 | 0.0700 | +0.3% | +0.2 |
| L1 boxes of 4 holding cells | 68.36 | 68.05 | -0.5% | -0.4 |
| L1 seam: cells within 1 inside | 17.64 | 17.40 | -1.4% | -0.8 |
| L1 seam: cells within 1 outside | 18.51 | 17.88 | -3.4% | -1.9 |
| L1 seam: joins across | 9.64 | 9.39 | -2.5% | -1.5 |
| L2 cells | 1,069 | 1,013 | -5.3% | -2.5 |
| L2 with 1 neighbours (per cell) | 0.1691 | 0.1691 | +0.0% | +0.0 |
| L2 with 2 neighbours (per cell) | 0.5436 | 0.5444 | +0.1% | +0.7 |
| L2 with 3 neighbours (per cell) | 0.2588 | 0.2582 | -0.2% | -0.6 |
| L2 with 4 neighbours (per cell) | 0.0285 | 0.0283 | -0.8% | -0.6 |
| L2 empty sites next to the shape (per cell) | 1.19 | 1.20 | +0.2% | +1.4 |
| L2 full 2x2 (per cell) | 0.0710 | 0.0704 | -0.8% | -0.9 |
| L2 boxes of 4 holding cells | 207 | 196 | -5.0% | -2.4 |
| L2 boxes of 16 holding cells | 26.69 | 25.70 | -3.7% | -2.1 |
| L2 seam: cells within 2 inside | 51.73 | 50.30 | -2.8% | -1.1 |
| L2 seam: cells within 2 outside | 53.63 | 52.82 | -1.5% | -0.6 |
| L2 seam: joins across | 14.08 | 13.83 | -1.7% | -0.7 |
| L3 cells | 3,096 | 2,921 | -5.7% | -2.6 |
| L3 with 1 neighbours (per cell) | 0.1690 | 0.1689 | -0.0% | -0.2 |
| L3 with 2 neighbours (per cell) | 0.5447 | 0.5448 | +0.0% | +0.1 |
| L3 with 3 neighbours (per cell) | 0.2582 | 0.2578 | -0.1% | -0.6 |
| L3 with 4 neighbours (per cell) | 0.0281 | 0.0285 | +1.3% | +1.6 |
| L3 empty sites next to the shape (per cell) | 1.20 | 1.20 | +0.0% | +0.4 |
| L3 full 2x2 (per cell) | 0.0704 | 0.0707 | +0.4% | +0.8 |
| L3 boxes of 4 holding cells | 599 | 566 | -5.6% | -2.6 |
| L3 boxes of 16 holding cells | 77.83 | 74.51 | -4.3% | -2.1 |
| L3 seam: cells within 2 inside | 77.56 | 73.48 | -5.3% | -1.8 |
| L3 seam: cells within 2 outside | 78.73 | 78.81 | +0.1% | +0.0 |
| L3 seam: joins across | 20.85 | 20.23 | -3.0% | -1.0 |
| L4 cells | 8,963 | 8,715 | -2.8% | -1.2 |
| L4 with 1 neighbours (per cell) | 0.1691 | 0.1693 | +0.1% | +1.0 |
| L4 with 2 neighbours (per cell) | 0.5440 | 0.5442 | +0.0% | +0.5 |
| L4 with 3 neighbours (per cell) | 0.2585 | 0.2582 | -0.1% | -1.0 |
| L4 with 4 neighbours (per cell) | 0.0284 | 0.0283 | -0.2% | -0.5 |
| L4 empty sites next to the shape (per cell) | 1.19 | 1.20 | +0.2% | +4.5 |
| L4 full 2x2 (per cell) | 0.0707 | 0.0704 | -0.3% | -1.1 |
| L4 boxes of 4 holding cells | 1,734 | 1,689 | -2.6% | -1.1 |
| L4 boxes of 16 holding cells | 226 | 223 | -1.6% | -0.7 |
| L4 seam: cells within 2 inside | 113 | 109 | -4.2% | -1.3 |
| L4 seam: cells within 2 outside | 114 | 113 | -1.1% | -0.3 |
| L4 seam: joins across | 30.33 | 29.44 | -3.0% | -0.9 |
| L5 cells | 26,181 | 25,092 | -4.2% | -1.7 |
| L5 with 1 neighbours (per cell) | 0.1693 | 0.1693 | -0.0% | -0.2 |
| L5 with 2 neighbours (per cell) | 0.5438 | 0.5445 | +0.1% | +2.6 |
| L5 with 3 neighbours (per cell) | 0.2585 | 0.2580 | -0.2% | -3.0 |
| L5 with 4 neighbours (per cell) | 0.0283 | 0.0283 | -0.1% | -0.3 |
| L5 empty sites next to the shape (per cell) | 1.19 | 1.20 | +0.1% | +6.3 |
| L5 full 2x2 (per cell) | 0.0707 | 0.0704 | -0.4% | -2.7 |
| L5 boxes of 4 holding cells | 5,079 | 4,876 | -4.0% | -1.7 |
| L5 boxes of 16 holding cells | 665 | 643 | -3.2% | -1.3 |
| L5 seam: cells within 2 inside | 167 | 159 | -4.9% | -1.6 |
| L5 seam: cells within 2 outside | 165 | 166 | +0.8% | +0.2 |
| L5 seam: joins across | 44.82 | 43.40 | -3.2% | -1.0 |
| floating pieces | 0.0000 | 0.0000 | +0.0% | +0.0 |
| floating pieces of 10 or more | 0.0000 | 0.0000 | +0.0% | +0.0 |
| cells in floating pieces | 0.0000 | 0.0000 | +0.0% | +0.0 |
| pieces | 23.09 | 21.53 | -6.7% | -2.0 |
| pieces of 2 to 9 cells | 7.27 | 6.92 | -4.8% | -1.1 |
| pieces of 10 to 99 | 6.57 | 5.85 | -10.9% | -2.7 |
| pieces of 100 or more | 4.98 | 4.64 | -6.7% | -1.8 |
| single cells | 4.27 | 4.12 | -3.6% | -0.8 |
| share in the largest piece | 0.8505 | 0.8465 | -0.5% | -0.4 |
| cells | 39,964 | 38,392 | -3.9% | -2.1 |

Leaving out the rows for floating pieces (which a real picture cannot have), 18 of 86 rows are beyond 2 s.e. and the largest is 6.3. Per level: 2.9, 3.0, 3.1, 3.4, 4.3 s; pieces joined again at the start of a stage (all stages of a level together) 1.9, 6.7, 13.1, 22.5, 35.7, with 2.7, 12.4, 26.7, 50.4, 77.0 cells; of those joins, longer than 8 cells 0.00, 0.04, 0.13, 0.26, 0.38; parts in the finer world's box 3.78, 3.62, 3.58, 3.63, 3.50, kept apart in the new world 2.88, 2.78, 2.64, 2.78, 2.69.

The refining alone (`analysis/e2e_test.py` with `N=0 S=0 ball=1024,250000` and the engine built for larger worlds, 312 roots): the real cells out to 1,024 from a root are halved and thinned 6 times and the rings refined again, with nothing drawn afresh, so that every world is the real shape at its scale and only the refining is on trial. Made against real, in per cent, with the difference in standard errors in brackets.

| ring | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| cells in the ring | -1.1% (-1.8) | +0.7% (+1.1) | +0.2% (+0.4) | -0.2% (-0.3) | -0.6% (-0.9) |
| cells in the two rows inside its inner seam | +1.7% (+0.8) | -4.1% (-2.9) | +0.7% (+0.4) | +2.9% (+1.7) | +1.1% (+0.7) |
| cells in the two rows outside it | +0.8% (+0.4) | +1.8% (+1.1) | +2.3% (+1.4) | +4.9% (+2.7) | +2.1% (+1.1) |
| joins across it | +1.2% (+0.6) | -3.0% (-1.7) | +0.1% (+0.1) | +2.8% (+1.5) | +1.4% (+0.7) |

(One row either side at the first seam.) 13 of the 86 rows of the full list are beyond 2 s.e., the largest 5.0; floating pieces a picture 0.000. The version before this one, on 104 roots: joins across the five seams -2.5, +2.2, +9.3, -8.8 and -18.7%, the two rows outside the fifth seam 18.5% short and the fifth ring 5.5% over.

Stage by stage (`analysis/stage_test.py`, 130 roots, again with nothing drawn afresh): cells in the picture of a level at the end of each stage, divided by the true cells at that resolution, in zones named by layers of the level's own world. The strip is the two layers beyond the ring that are refined with it and not kept; it is thin at its far side, which is the edge of what is refined.

| level | stage | first two layers of the ring | middle of the ring | last two layers of the ring | strip |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 (last) | 0.987 ± 0.019 | 0.999 ± 0.015 | 1.014 ± 0.020 | 0.977 ± 0.016 |
| 2 | 1 | 1.000 ± 0.016 | 0.992 ± 0.011 | 0.992 ± 0.015 | 0.983 ± 0.014 |
| 2 | 2 (last) | 1.041 ± 0.020 | 0.990 ± 0.014 | 1.002 ± 0.017 | 0.928 ± 0.015 |
| 3 | 1 | 0.984 ± 0.015 | 0.998 ± 0.012 | 0.997 ± 0.015 | 0.960 ± 0.013 |
| 3 | 2 | 1.013 ± 0.015 | 1.003 ± 0.012 | 0.995 ± 0.015 | 0.943 ± 0.014 |
| 3 | 3 (last) | 1.007 ± 0.014 | 1.000 ± 0.013 | 0.988 ± 0.016 | 0.916 ± 0.014 |

Far out (`analysis/mass_test.py`, 130 roots): real shapes counted in pixels 8 cells wide, out to 128 cells from a root, and the same shapes halved and thinned until a cell is a pixel, half a pixel or a quarter of one, counted in the same pixels. Brightness is a pixel's count over the mean count of lit pixels.

| cells to a pixel side | lit pixels | of the real number | spread of brightness (s.d. / mean) | brightness against real, pixels lit in both (r) | lit here only | lit in the real only |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 183 ± 4 | 55% | 0.00 |  | 25 | 173 |
| 2 | 258 ± 6 | 78% | 0.41 | 0.33 | 17 | 89 |
| 4 | 307 ± 7 | 93% | 0.53 | 0.77 | 11 | 34 |
| 8 (real cells) | 330 ± 7 | 100% | 0.60 |  | 0 | 0 |

Whether a table could supply what the dots lack (`analysis/mass_model.py`): for every dot and every empty site beside one, a number of real cells is drawn from what real shapes have at sites with the same dots around them (tables from 130 roots, drawn for 130 others). Lit pixels against the real cells' lit pixels, and the spread of brightness, for the draws and for the dots as they are shown now.

| a dot, in real cells | a pixel, in dots | lit: drawn | lit: dots | spread: real | drawn | dots |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 1 | +0.2% | -28.7% | 0.41 | 0.41 | 0.00 |
| 2 | 2 | +8.2% | -15.7% | 0.54 | 0.56 | 0.41 |
| 2 | 4 | +4.2% | -7.3% | 0.60 | 0.61 | 0.53 |
| 4 | 1 | +0.0% | -40.0% | 0.54 | 0.54 | 0.00 |
| 4 | 2 | +9.7% | -22.2% | 0.60 | 0.61 | 0.41 |
| 4 | 4 | +7.0% | -10.5% | 0.63 | 0.66 | 0.53 |
| 8 | 1 | -0.2% | -44.8% | 0.60 | 0.60 | 0.00 |
| 8 | 2 | +9.6% | -24.8% | 0.63 | 0.64 | 0.41 |
| 8 | 4 | +7.4% | -11.3% | 0.64 | 0.67 | 0.52 |

Rings against the average real ring (`analysis/ring_pop.py`; 2600 real roots, each ring's average known to 0.7%): the made rings of the full run -0.1% (-0.1), -3.8% (-2.1), -4.8% (-2.7), -2.6% (-1.5), -3.1% (-1.7); the real rings of the same 390 roots +0.0% (+0.0), +1.6% (+0.8), +1.0% (+0.5), +0.2% (+0.1), +1.1% (+0.6) (in brackets, standard errors).

The step from one world to the next (`analysis/world_test.py`, 520 roots): cells by layer round the root (layer r is the square ring at distance r) in real shapes at their own scale, in the same shapes halved, and in the first world made from them. Halving leaves one cell for 2.960 real ones. The sampler holds the layers below 10.

| layer | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real | 9.45 | 10.43 | 11.46 | 12.09 | 12.96 | 13.89 | 14.73 | 15.20 | 15.79 | 16.11 | 16.64 |
| halved | 9.51 | 10.54 | 11.21 | 11.89 | 12.50 | 13.65 | 14.40 | 14.97 | 16.01 | 16.87 | 17.62 |
| world | 9.51 | 10.54 | 11.21 | 11.89 | 12.50 | 13.49 | 14.39 | 15.04 | 16.02 | 16.65 | 17.28 |

In the box: real 61.0 ± 0.7, halved 60.1 ± 0.7. In the ring: real 117.4 ± 1.7, halved 117.9 ± 1.7, world 117.3 ± 1.8 (world less halved, root by root, -0.7 ± 1.2). In the first two layers drawn afresh: halved 28.1 ± 0.5, world 27.9 ± 0.5 (-0.17 ± 0.36). The finer world crosses the edge of the box at 9.1 sites of the halved one.

Classifiers trained to tell the made pictures from the real ones (`analysis/tell_test.py`; 390 of each, scored on pictures from samples they were not trained on; area under the ROC curve, 0.5 = cannot tell): a logistic regression on about 150 summary numbers per picture, 0.632 ± 0.019 (59% of pictures called correctly); the same without the numbers that describe how the picture falls into pieces, 0.631 ± 0.021 (59% of pictures called correctly); a small network looking at 9×9 windows centred on cells, 3,000 per picture, 0.674 ± 0.018 (63% of pictures called correctly). The same classifiers set to tell the real pictures of half the samples from those of the other half, where there is nothing to find, score 0.534 ± 0.029, 0.524 ± 0.031 and 0.586 ± 0.028.

Where strands rejoin (`analysis/part_test.py`, 208 roots): of the pieces into which the cells within 16 of a root fall, the number still separate within R cells, in real samples and in the first world made from them (at R/2 of its own cells, and by the world's own rule for what is joined).

| R | 16 | 18 | 20 | 24 | 32 | 48 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|
| real | 3.73 ± 0.13 | 2.19 ± 0.08 | 1.76 ± 0.06 | 1.39 ± 0.04 | 1.18 ± 0.03 | 1.05 ± 0.01 | 1.03 ± 0.01 |
| made | 2.84 ± 0.10 | 1.96 ± 0.07 | 1.67 ± 0.05 | 1.36 ± 0.04 | 1.13 ± 0.02 | 1.06 ± 0.02 | 1.03 ± 0.01 |

Run 40 levels out at the page's resolution (`analysis/out_deep.js`): first half of the levels against second half, cells of the world in its ring 129 and 131; parts kept apart 3.00 and 3.15; and, from level 7 on, where the picture is as fine as it gets, cells in the finished picture of the ring 82,418 and 78,430, empty sites next to the shape per cell 1.1950 and 1.1963, cells with one neighbour per cell 0.1691 and 0.1692. A level takes 8.1 s, 2.9 s of it the sampler.

Worlds only, 8 runs of 240 levels from cells of 8 different samples (`scripts/out_drift.sh`). Means over all levels: cells of the world in its box 60.3 ± 0.8, in its ring 116.6 ± 1.6, parts kept apart 2.68 ± 0.04; real shapes have 61.1 and 117.9 cells within 8 and from 8 to 16 cells of one of their cells (1,560 cells of 26 samples). A straight line through each run changes by +0.3 ± 1.0 cells per 100 levels in the box and +0.4 ± 2.0 in the ring. Quarter by quarter the ring holds 116, 117, 115, 118.
<!-- /zoom-out table -->

What is not right, as far as these tests show.

- *Rings beyond the first are a few per cent short in the full run.* Set against the
  average real ring over 2,600 roots they hold 3.8, 4.8, 2.6 and 3.1% too few cells (1.5 to
  2.7 s.e. each; `analysis/ring_pop.py`), more toward their outer edges, where the run that
  only refines has every ring within about 1%. The earlier version's raised weight had been
  making up for this along with its own edges; now that the refining is unbiased, it shows.
  Where it comes from is not settled. It is not the refining: counted on 200 roots, a stage
  makes as many cells of one with the sampler as without (2.92 and 2.92, then 2.97 and
  2.98, at the second level's two stages; 2.87 and 2.88, 2.95 and 2.95, 2.95 and 2.96 at
  the third's; each known to about 0.02). That leaves the worlds, and worlds drawn by the
  sampler, set layer by layer against halved real shapes for the same roots, agree only to
  within a noise of 3 to 4%, which cannot show a difference of 3.
- *Seams.* Within about 5% on every count in both runs, and no seam count of the full run is
  beyond 2 s.e.; the joins across are 2 to 3% low at all five seams there (each within
  1.5 s.e.), which goes with the rings being short. In the run that only refines, the two
  rows outside the fourth seam are 4.9% over (2.7 s.e.) and the two inside the second 4.1%
  under (2.9 s.e.), out of twenty seam counts.
- *Joins are made by hand.* At the start of every stage the ring's own strands are joined
  to the band again where they no longer meet it: about 36 pieces a picture at the fifth
  level, 77 cells in all, and for one join in a hundred the way is longer than 8 cells, a
  straight strand that the walk afterwards only partly roughens.
- *The fine texture is very slightly too open.* With the last stage's weight just under
  1/λ the first ring has 1% more empty sites beside the shape per cell than real ones and
  6% fewer cells with four neighbours (in the run that only refines); at the outer rings
  the same counts are within 0.2% and 2%, which that run still resolves (5 s.e. for the
  first).
- *Classifiers still tell the pictures from real ones more often than the control does.*
  On five levels the summary classifier calls 59% of pictures correctly (0.632 ± 0.019;
  the control 51%, 0.534 ± 0.029) and the one that looks at small windows 63% (0.674 ±
  0.018; the control 56%, 0.586 ± 0.028). The version before this one, put through the same
  five levels with the same roots, scored 65% and 61% (0.701 ± 0.018 and 0.659 ± 0.020):
  the summary numbers have come more than a third of the way to the control, the windows
  not at all. On the first four levels alone the two score 56% and 55% (0.562 ± 0.021 and
  0.589 ± 0.019; the control there 51% and 55%), so what the windows find is mostly in the
  fifth ring, the only one refined five times, and what it is has not been found. (That
  version's counts on five levels: rings +0.5, -2.3, -2.8, +2.0 and +2.9%; joins across
  the five seams 0.0, -2.1, -7.6, -11.3 and -16.2%, the last at 5.3 s.e.)
- *Far out, the picture is not what a real shape would look like from that far.* Past the
  sixth level a texel stands for a block of real cells, and the page draws the world's
  cells, which keep the look of single cells at every scale. A real shape seen from as far
  would be drawn by how many of its cells fall in each pixel, and that is fuller and less
  even: counted in pixels 8 cells wide, a shape thinned to two cells per pixel side (about
  what the page shows) lights 78% of the pixels the real cells light, and its brightness
  follows theirs with r = 0.33; at four per side, 93% and 0.77 (`analysis/mass_test.py`,
  table above). So between the sixth level and about the ninth a real zoom would go on
  changing, from the look of cells to the look of the limit, and the page's does not. To
  follow it each texel would have to carry the number of real cells it stands for. A table
  goes about halfway (`analysis/mass_model.py`, above): drawing that number for each dot
  and each empty site beside one, from what real shapes have at sites with the same dots
  around them, gets the counts dot by dot right by construction, but in pixels two dots
  wide it lights 8 to 10% too many where the dots alone light 16 to 25% too few. The draws
  are independent and the real cells beside a strand are not: they lie along some strands
  and not along others. Doing it properly means refining three or four more times below the
  pixel, 26 times the cells or more, or a model that keeps that coherence; neither is done.
- *Numbers that are set, not derived:* the weight of the thinning walk and the two weights
  of the refining walk, all by the requirement that a shape at twice the scale has 2.95
  times the cells; the size below which a part is merged; the depth of the strip and of
  the band.
- A world is finite: 16 to 32 of its cells from the root it holds 3% fewer cells than a
  real shape. Beyond five or six levels there is no real sample large enough to compare
  with, only the check that nothing drifts. And in shapes of a dozen cells the sampler with
  a held box cannot reach a few shapes in which the only path between two held cells has no
  cell to spare (11 of 1,801 at 11 cells, found by the same exhaustive count), and below 12
  cells it cannot take that path round the other side of the box; neither matters at
  16,000.

Three checks (`scripts/zoom_relax.sh`, `scripts/zoom_test.sh`; tables by
`analysis/zoom_report.py`):

<!-- zoom tables -->
| sweeps | cells ÷ original, 16,000 doubled | cells ÷ original, 64,000 doubled | all measurements pooled | 2×2 blocks | perimeter | tips | branchings | boxes changed: 2 | 8 | 32 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3.2727 ± 0.0023 | 3.2719 ± 0.0008 | 3.2720 ± 0.0008 | 2.325 ± 0.005 | 0.882 ± 0.000 | 0.834 ± 0.001 | 1.356 ± 0.002 | 13% | 3% | 1% |
| 2 | 3.1076 ± 0.0024 | 3.1087 ± 0.0010 | 3.1085 ± 0.0009 | 1.392 ± 0.004 | 0.952 ± 0.000 | 0.999 ± 0.001 | 1.127 ± 0.001 | 18% | 4% | 1% |
| 4 | 3.0201 ± 0.0020 | 3.0190 ± 0.0014 | 3.0194 ± 0.0012 | 1.063 ± 0.003 | 0.989 ± 0.000 | 1.025 ± 0.001 | 1.031 ± 0.001 | 24% | 5% | 1% |
| 8 | 2.9806 ± 0.0035 | 2.9815 ± 0.0020 | 2.9813 ± 0.0017 | 1.005 ± 0.003 | 0.998 ± 0.000 | 1.015 ± 0.001 | 1.007 ± 0.002 | 29% | 6% | 2% |
| 16 | 2.9625 ± 0.0064 | 2.9602 ± 0.0025 | 2.9597 ± 0.0005 | 0.996 ± 0.003 | 1.001 ± 0.000 | 1.006 ± 0.001 | 1.000 ± 0.001 | 36% | 8% | 2% |
| 32 | 2.9568 ± 0.0047 | 2.9500 ± 0.0030 | 2.9521 ± 0.0006 | 0.994 ± 0.004 | 1.001 ± 0.000 | 1.002 ± 0.001 | 0.999 ± 0.001 | 43% | 10% | 3% |
| 64 | 2.9435 ± 0.0058 | 2.9420 ± 0.0044 | 2.9473 ± 0.0008 | 0.996 ± 0.003 | 1.000 ± 0.001 | 1.000 ± 0.001 | 0.999 ± 0.001 | 52% | 12% | 3% |
| 128 | 2.9508 ± 0.0104 | 2.9352 ± 0.0055 | 2.9386 ± 0.0049 | 0.995 ± 0.004 | 1.000 ± 0.000 | 1.001 ± 0.001 | 0.999 ± 0.002 | 61% | 15% | 4% |
| 256 | 2.9416 ± 0.0112 | 2.9423 ± 0.0050 | 2.9435 ± 0.0011 | 1.001 ± 0.003 | 1.000 ± 0.001 | 1.000 ± 0.001 | 0.999 ± 0.001 | 72% | 19% | 5% |
| 512 | 2.9369 ± 0.0138 | 2.9493 ± 0.0069 | 2.9468 ± 0.0062 | 0.997 ± 0.003 | 1.001 ± 0.000 | 1.001 ± 0.001 | 0.999 ± 0.002 | 85% | 24% | 7% |
| 1024 | 2.9218 ± 0.0258 | 2.9505 ± 0.0094 | 2.9471 ± 0.0089 | 1.000 ± 0.004 | 1.000 ± 0.001 | 1.000 ± 0.001 | 0.999 ± 0.002 | 97% | 29% | 8% |

Prediction for the first columns: 2^(1/ν) = 2.9477 ± 0.0025 for ν = 0.6412 ± 0.0005. The local statistics are per cell, relative to sampled shapes, for the 64,000-cell samples; boxes are in cells of the original.

| sweeps | clumps ≥ 4 | clumps ≥ 8 | clumps ≥ 16 | clumps ≥ 32 | clumps ≥ 64 |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.53 | 3.36 | 7.33 | 8.29 | 5.55 |
| 4 | 1.04 | 1.11 | 1.11 | 1.14 | 1.48 |
| 16 | 1.00 | 1.00 | 1.01 | 0.90 | 0.80 |
| 64 | 1.00 | 1.00 | 1.00 | 0.90 | 0.89 |
| 256 | 1.00 | 1.01 | 1.02 | 0.97 | 1.16 |

Loop-connected clumps of at least that many cells, relative to the equilibrium number (`analysis/clumps.json`), in 12 refined samples.

| sweeps per level | shapes | cells | growth at levels 1 to 4 | per-cell counts (12): beyond 2 s.e. | largest | the other 48: beyond 2 s.e. | largest | largest differences |
|---:|---:|---:|---|---:|---:|---:|---:|---|
| 16 | 32 | 1,048,092 ± 7,077 | 2.9599, 2.9630, 2.9577, 2.9600 | 8 | 22.5 | 1 | 2.1 | leaves +22.5 (+0.66%), deg2 -20.7 (-0.38%), noncut +9.6 (+0.32%) |
| 32 | 32 | 1,037,350 ± 7,864 | 2.9537, 2.9535, 2.9505, 2.9525 | 7 | 9.8 | 3 | 2.6 | leaves +9.8 (+0.29%), deg2 -9.1 (-0.17%), noncut +5.5 (+0.18%) |
| 64 | 32 | 1,031,145 ± 11,354 | 2.9513, 2.9458, 2.9478, 2.9476 | 3 | 3.1 | 2 | 2.6 | deg2 -3.1 (-0.06%), cen_m1 -2.6 (-4.12%), deg3 +2.6 (+0.08%) |
| 256 | 32 | 1,026,739 ± 13,404 | 2.9487, 2.9422, 2.9468, 2.9424 | 0 | 1.1 | 3 | 2.9 | trunk01 +2.9 (+2.14%), trunk03 +2.9 (+3.63%), bw +2.1 (+7.46%) |

Samples of 13,650 cells refined four times, against the equilibrium runs at 1,024,000 cells (1295 readings in 24 blocks). A difference is in standard errors of the mean of that many equilibrium samples.
<!-- /zoom tables -->

Reading the tables. The cell count comes to rest between 2.94 and 2.95 times the original:
2.952, 2.947 and 2.943 after 32, 64 and 256 sweeps (each ± 0.001), against the 2^(1/ν) =
2.948 ± 0.003 that scaling predicts for a shape of twice the width, which the procedure is
never told. It is still creeping downward at the longest runs, so the agreement is to about
two parts in a thousand and no better. Doubling leaves several times too many large loops
(each survives at twice the width, and loops twice as wide should be about twenty times rarer
per cell); they are gone within 16 sweeps. After four refinements in a row, the 48 observables
that describe the shape at large cannot be told from those of sampled million-cell shapes at
any of the sweep counts tried. The largest difference among them, the trunk at 2.9 standard
errors, is inherited: the 32 starting samples themselves have trunks 1.6% (`trunk01`) and 3.4%
(`trunk03`) longer than the average at 13,650 cells, and the refined shapes 2.1% and 3.6%.
The 12 per-cell counts of small features, which a million cells fix to a few parts in ten
thousand, show how long the fine texture takes to settle: the worst is off by 0.9% at 16
sweeps, 0.5% at 32, under 0.1% at 64, and nothing is distinguishable at 256.
`analysis/zoom_deep.js` runs the browser's engine (32 sweeps a level) 400 levels deep and
finds its per-cell statistics where they started.

What this is not: a construction of the scaling limit of polyominoes, or evidence beyond four
doublings, where there are no samples to compare with. Each level is a stand-in for a sample
of that many cells seen around one point.

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
src/zoom.cpp           refinement to twice the resolution (cells added and removed at weight 1/λ)
src/outward.c          the engine of the page's zoom out: the sampler with a held box, halving, the walks on a picture (builds to page/outward.wasm)
src/ball.cpp           a connected piece of a sample around a cell, for the tests of the zoom out
src/features.h         60 shape observables (geometry, local counts, cut structure, graph distances)
src/enumerate.cpp      exact averages of them over every polyomino of n cells
src/leath.cpp          exact percolation clusters of a given size
src/grow.cpp           uniform polyominoes from a resampled population of growing clusters
ref/manyomino_v1.cpp   naive reference implementation
proof/                 a quasi-polynomial bound on the relaxation time, and its machine check
scripts/               production ladders (resumable); convergence runs
analysis/              autocorrelation times, data reduction, rendering, page build
page/                  source of the explainer page; page/zoom.js is the in-browser refinement, page/outward.js the zoom out
video/                 music video generator
```
