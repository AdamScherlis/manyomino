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
