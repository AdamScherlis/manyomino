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
