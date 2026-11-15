#!/bin/bash
# Runs behind analysis/convergence.py.  Four pipelines, one core each, about 80 minutes.
#   u<n>_<A|B>   uniform chain at n cells, full observable library every n moves, started from the
#                final shape of production ladder A or B (independent, already equilibrated)
#   t<p>_<n>_*   chain tilted by (1-p)^perimeter, started from an exact sample
#   x<p>_<n>     exact independent percolation clusters of n cells at density p (src/leath.cpp)
#   g<n>_<k>     uniform polyominoes by population growth, no Markov chain (src/grow.cpp)
set -u
cd "$(dirname "$0")/.."
M=./manyomino_feat; D=conv; mkdir -p $D
PC=0.592746
uni() { $M -load prod/ckpt_$2_n$1.bin -pcp 0.4 -moves 1e18 -tsec $3 -every $1 -seed $4 -out $D/u$1_$2.ts -feat $D/u$1_$2.feat 2> $D/u$1_$2.log; }
tilt() { $M -load $D/x$1_$2_000$3.bin -p $1 -pcp 0.4 -moves 1e18 -tsec $4 -every $2 -seed $5 -out $D/t$1_$2_$3.ts -feat $D/t$1_$2_$3.feat 2> $D/t$1_$2_$3.log; }
exact() { ./leath -n $2 -p $1 -count $3 -tsec $4 -seed $5 -out $D/x$1_$2.feat -snap $D/x$1_$2 -nsnap 2 2> $D/x$1_$2.log; }
( uni 64000 A 4500 101 ) &
( uni 64000 B 4500 102 ) &
( exact $PC 1000 100000 2400 11; exact 0.55 500 20000 600 12; uni 16000 A 1500 103 ) &
( uni 16000 B 1500 104; uni 4000 A 400 105; uni 4000 B 400 106; uni 1000 A 200 107; uni 1000 B 200 108; uni 250 A 100 109; uni 250 B 100 110
  exact 0.5 250 100000 600 13
  tilt 0.5 250 0 120 201; tilt 0.5 250 1 120 202
  tilt 0.55 500 0 150 203; tilt 0.55 500 1 150 204
  tilt $PC 1000 0 300 205; tilt $PC 1000 1 300 206
  for k in 1 2 3 4 5 6; do ./grow -nmax 250 -pop 500000 -step 4 -sizes 250 -seed $k -out $D/g250_$k 2> $D/g250_$k.log; done ) &
wait
