| clump of at least | per cell | lasts, 4,000 cells | lasts, 16,000 cells | lasts, 64,000 cells | lasts with `-relA 1` |
|---:|---:|---:|---:|---:|---:|
| 4 cells | 0.0486 | 4.4 ± 0.1 | 4.4 ± 0.0 | 4.4 ± 0.2 | 4.3 ± 0.2 |
| 8 cells | 0.00663 | 3.4 ± 0.0 | 3.5 ± 0.1 | 3.1 ± 0.1 | 2.6 ± 0.1 |
| 16 cells | 0.000575 | 5.0 ± 0.1 | 5.1 ± 0.0 | 4.9 ± 0.3 | 1.8 ± 0.0 |
| 32 cells | 4.86e-05 | 11.2 ± 0.2 | 11.1 ± 0.2 | 12.7 ± 0.8 | 1.6 ± 0.0 |
| 64 cells | 4.95e-06 | 24.0 ± 1.7 | 19.9 ± 0.8 | 19.9 ± 1.1 | 1.6 ± 0.1 |
| 128 cells | 5.48e-07 | 41.5 ± 3.6 | 37.6 ± 5.1 | 32.0 ± 10.4 | 1.7 ± 0.1 |
| 182 cells | 2.08e-07 | 52.7 ± 4.8 | 53.5 ± 9.0 |  | 1.7 ± 0.2 |
| 256 cells | 8.45e-08 | 67.1 ± 6.8 | 72.6 ± 27.1 |  | |
| 363 cells | 6.28e-08 | 79.1 ± 10.0 |  |  | |

Times are integrated autocorrelation times, in moves per cell, of the number of clumps of at least that size (2,381,877 samples at 4,000 cells, 755,268 samples at 16,000 cells, 50,564 samples at 64,000 cells, one per move per cell; errors by leaving out blocks).

| cells | largest clump, mean | its time | with `-relA 1`: mean | time |
|---:|---:|---:|---:|---:|
| 4,000 | 25.5 ± 0.1 | 39.1 ± 5.8 | 25.4 ± 0.1 | 1.8 ± 0.1 |
| 16,000 | 38.4 ± 0.2 | 40.9 ± 12.7 | 38.2 ± 0.2 | 1.6 ± 0.1 |
| 64,000 | 56.8 ± 1.1 | 33.1 ± 9.0 | 60.6 ± 2.7 | 2.2 ± 1.8 |

The largest clump of K independent stretches of a smaller record, taken together as one shape of K times the size (time: median over 24 relative alignments of the stretches, with the larger of their spread and the block error):

| stands for | built from | largest clump, mean | its time | range over alignments |
|---:|---|---:|---:|---|
| 256,000 | 64 × 4,000 | 98 ± 3 | 62 ± 10 | 50 to 74 |
| 1,024,000 | 256 × 4,000 | 163 ± 11 | 65 ± 27 | 49 to 147 |
| 4,096,000 | 1024 × 4,000 | 272 ± 23 | 72 ± 16 | 36 to 95 |
| 1,024,000 | 64 × 16,000 | 140 ± 8 | 76 ± 36 | 57 to 164 |
| 4,096,000 | 256 × 16,000 | 221 ± 15 | 84 ± 36 | 48 to 127 |

Other observables with `-relA 1` against the default (`conv/` chains), moves per cell:

| observable | 16,000: default | 16,000: `-relA 1` | 64,000: default | 64,000: `-relA 1` |
|---|---:|---:|---:|---:|
| trunk01 | 40.5 ± 3.8 | 30.9 ± 6.8 | 39.0 ± 4.3 | 37.1 ± 10.1 |
| trunk03 | 34.7 ± 5.2 | 24.7 ± 5.0 | 31.7 ± 2.4 | 25.0 ± 8.7 |
| trunk10 | 21.6 ± 1.8 | 17.9 ± 3.9 | 23.8 ± 2.5 | 35.0 ± 5.2 |
| perim | 20.3 ± 1.5 | 16.8 ± 1.6 | 18.3 ± 1.8 | 16.3 ± 2.5 |
| rg2 | 14.7 ± 0.9 | 14.5 ± 1.8 | 13.3 ± 0.9 | 9.8 ± 3.4 |
| blob_max | 27.7 ± 8.3 | 1.6 ± 0.1 | 30.3 ± 9.6 | 2.1 ± 1.7 |
