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
gyration) grows like `n^2.5` moves, because mass can only diffuse through the shape one cell
at a time. With transplants it is about `14 n` moves. See [Results](#results).

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
| `-shape rect\|bar` | starting shape (a near-square block, or a 1×N bar) |
| `-load F [-inflate k]` | start from a saved shape, optionally blown up k× |
| `-every E`, `-out F` | write observables every E moves |
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
with `P(K ≥ r) = r⁻²`. Let `s = u + d` and `Q = P \ c`. If `s` is empty, continue with
probability `1/k`, where `k` is the number of neighbours of `s` in `Q`. Require that the
neighbours of `c` in `Q` are connected to each other within the radius-`K` box around `c`,
and likewise for `s`. Then move `c` to `s`. Given `K`, every valid pair `(c, s)` has
probability `1/4n²` and the conditions are symmetric in `c` and `s`. The condition at `c`
implies `Q` is connected. Most draws of `K` are small, so the test is cheap, but arbitrarily
large `K` has positive probability, which is what makes the chain irreducible. `-relK 0` sets
`K = ∞` always: the textbook kernel.

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
* Independent ladders from different seeds agree at every size.
* Counting removable cells recovers Klarner's constant (below).

## Results

RESULTS_PLACEHOLDER

## Layout

```
src/manyomino.cpp      the sampler
ref/manyomino_v1.cpp   naive reference implementation
scripts/               production ladders (resumable)
analysis/              autocorrelation times, data reduction, rendering, page build
page/                  source of the explainer page
```
