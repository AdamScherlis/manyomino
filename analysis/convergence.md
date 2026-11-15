| reference (no Markov chain) | details | comparisons | rms z | largest | beyond 2 s.e. |
|---|---|---:|---:|---:|---:|
| all 7,204,874 shapes of 14 cells, uniform | 2,000,001 chain samples; precision 7.5e-04 sd | 60 | 1.03 | 2.33 | 2 |
| all 104,592,937 shapes of 16 cells, uniform | 2,000,001 chain samples; precision 7.7e-04 sd | 60 | 0.70 | 2.23 | 1 |
| all 7,204,874 shapes of 14 cells, weighted (1-p)^perimeter, p = 0.5 | 2,000,001 chain samples; precision 7.6e-04 sd | 60 | 1.09 | 2.24 | 3 |
| all 7,204,874 shapes of 14 cells, weighted (1-p)^perimeter, p = 0.592746 | 2,000,001 chain samples; precision 8.6e-04 sd | 60 | 0.87 | 1.84 | 0 |
| 20,000 exact percolation clusters of 500 cells at p = 0.55 | 690,105 chain samples; learned test z = -0.44; quartiles of each distribution: 1 of 176 beyond 2 s.e., largest 2.63 | 60 | 0.81 | 2.00 | 0 |
| 8 grown populations of 200,000 at 250 cells, uniform | count check −0% ± 0.9%; 18,453 distinct ancestors at n/2; precision 0.004 sd | 60 | 0.97 | 2.96 | 2 |
| 8 grown populations of 200,000 at 500 cells, uniform | count check +0% ± 1.6%; 6,636 distinct ancestors at n/2; precision 0.008 sd | 60 | 0.92 | 3.07 | 2 |
| 8 grown populations of 199,000 at 1,000 cells, uniform (fails its own count check) | count check −11% ± 3.5%; 2,132 distinct ancestors at n/2; precision 0.021 sd | 60 | 1.28 | 3.64 | 6 |
| 6 grown populations of 100,000 at 2,000 cells, uniform (fails its own count check) | count check −36% ± 9.2%; 296 distinct ancestors at n/2; precision 0.066 sd | 60 | 2.35 | 9.71 | 15 |
| 17 observables that are exactly zero by symmetry, at 3 sizes from 1,000 to 16,000 cells | two chains per size | 51 | 0.76 | 2.09 | 1 |
| control: chain at p = 0.51 against exact clusters at p = 0.50, 250 cells | learned test z = +6.56 | 60 | 6.79 | 17.47 | 28 |

z is (chain − reference) / standard error for each observable in `src/features.h`; chain errors are batch means. Over the 8 tests with a trustworthy reference: 471 comparisons, 11 beyond 2 standard errors (2.3%; 4.6% expected if every one were independent and exact), 1 beyond 3, largest 3.07 (cen_m2, grown population, 500 cells). Observables are correlated, so these are not 471 independent tests.

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

Autocorrelation times in moves per cell (integrated, Sokal window), two independent chains per size:

| cells | samples (one per n moves) | Rg² | slowest single observable | slowest combination found (held out) | A-vs-B learned test z |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 349,060 | 14.8 ± 0.5 | 25.8 ± 1.1 (trunk03) | 30.9 ± 2.0 | +0.49 |
| 4,000 | 125,954 | 14.4 ± 0.8 | 35.4 ± 2.9 (trunk01) | 43.1 ± 5.5 | +1.03 |
| 16,000 | 76,476 | 14.4 ± 1.0 | 41.3 ± 4.7 (trunk01) | 47.9 ± 8.1 | +1.40 |

Fits from 1,000 cells up, tau/n ∝ n^x: Rg² x = -0.012 ± 0.025; slowest single observable x = +0.185 ± 0.039; slowest combination x = +0.178 ± 0.059.
