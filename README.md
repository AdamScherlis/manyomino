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
at a time. With transplants it is about `15 n` moves. See [Results](#results).

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

Two independent ladders (seeds A and B) at every other size, a third (C) at the sizes between;
one core each, about 16 core-hours in total. `analysis/summarize.py` and `analysis/build_page.py`
regenerate everything below from the raw time series, which are not in the repository
(`scripts/launch_all.sh` reruns them).

| cells | chains | run length / tau | <Rg^2> | <Rg> | tau / n | perimeter / n | removable / n | holes / n | lambda estimate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 125 | 1 | 81,208 | 90.948 ± 0.17 | 9.4562 ± 0.0084 | 14.5 | 1.23006 | 0.30302 | 0.02094 | 4.0602 ± 0.0032 |
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
| 1,024,000 | 2 | 215 | 9.2941e+06 ± 2.3e+05 | 3029.3 ± 37 | 5.4 | 1.19485 | 0.29437 | 0.02311 | 4.0590 ± 0.0016 |

Each run starts from a 2x-inflated sample of the size below; the first 300 n moves (and at least a fifth of the run) are discarded. tau is the integrated autocorrelation time of Rg^2 in moves, and run length counts all moves. Above about 100,000 cells the runs are only tens of tau long after burn-in, so tau itself is poorly determined there (hence the scatter in tau / n).

Relocate-only autocorrelation time (moves), n = 50: 3.84e+03, 100: 2.41e+04, 200: 1.41e+05, 400: 9.33e+05, 800: 6.09e+06, 1600: 1.54e+07; fit tau ∝ n^2.62 (±0.04).
With transplants: tau ∝ n^0.99 (±0.03), about 15 moves per cell.
Size exponent from n ≥ 4,000: nu = 0.6402 ± 0.0014 (literature 0.6412 ± 0.0005).
Klarner's constant from n ≥ 1,000: lambda = 4.0619 ± 0.0006 (accepted 4.0625696; chi-square of the sizes about their mean 14.0 for 10 degrees of freedom, error scaled accordingly).

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
ref/manyomino_v1.cpp   naive reference implementation
scripts/               production ladders (resumable)
analysis/              autocorrelation times, data reduction, rendering, page build
page/                  source of the explainer page
video/                 music video generator
```
