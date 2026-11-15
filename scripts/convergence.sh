#!/bin/bash
# Runs behind analysis/convergence.py: four cores, about four hours.  (The committed results were
# produced by these commands launched in several pieces; this file collects them in one place.)
#
#   exact_<n>[_p<p>]  exact averages over every polyomino of n cells, optionally weighted by
#                     (1-p)^perimeter (src/enumerate.cpp);  e<n>[_p<p>]: the chain at that size
#   x<p>_<n>          exact independent percolation clusters of n cells at density p (src/leath.cpp)
#   t<p>_<n>_<k>      chain weighted by (1-p)^perimeter, started from one of those exact clusters
#   gk*, gs*          uniform polyominoes from grown populations, no Markov chain (src/grow.cpp)
#   u<n>_<A|B|C|D>    uniform chain at n cells, the full observable library every n moves, started
#                     from the final shape of a production ladder (already equilibrated)
#   s<n>_<bar|rect>   uniform chain started from a straight bar / a solid square
#   m250_<k>          eight uniform chains at 250 cells from one shape, different seeds
#   r0.5_250_<k>      eight more weighted chains at 250 cells, each from its own exact cluster
set -u
cd "$(dirname "$0")/.."
g++ -O2 -std=c++17 -o manyomino_feat src/manyomino.cpp
g++ -O2 -std=c++17 -o leath src/leath.cpp
g++ -O2 -std=c++17 -o grow src/grow.cpp
g++ -O2 -std=c++17 -o enumerate src/enumerate.cpp
M=./manyomino_feat; D=conv; mkdir -p $D
PC=0.592746
uni() { $M -load prod/ckpt_$5_n$1.bin -pcp 0.4 -moves 1e18 -tsec $3 -every $1 -seed $4 -out $D/u$1_$2.ts -feat $D/u$1_$2.feat 2> $D/u$1_$2.log; }   # n name seconds seed ladder
tilt() { $M -load $D/x$1_$2_000$(($3 % 2)).bin -p $1 -pcp 0.4 -moves 1e18 -tsec $4 -every $2 -seed $5 -out $D/t$1_$2_$3.ts -feat $D/t$1_$2_$3.feat 2> $D/t$1_$2_$3.log; }
exact() { ./leath -n $2 -p $1 -count $3 -tsec $4 -seed $5 -out $D/x$1_$2.feat -snap $D/x$1_$2 -nsnap 2 2> $D/x$1_$2.log; }
start() { $M -n $1 -shape $2 -pcp 0.4 -moves 1e18 -tsec $3 -every $1 -seed $4 -out $D/s$1_$2.ts -feat $D/s$1_$2.feat 2> $D/s$1_$2.log; }
grow() { ./grow -nmax $1 -pop $2 -kappa 1.8 -sizes $3 -seed $4 -out $D/$5 2> $D/$5.log; }

# ---- core 1 and 2: the long uniform chains
( uni 64000 A 4500 101 A; uni 256000 A 3300 201 A; uni 64000 C 2400 301 A; start 16000 bar 1300 621 ) &
( uni 64000 B 4500 102 B; uni 256000 B 3300 202 B; uni 64000 D 2400 302 B; start 16000 rect 1300 622 ) &

# ---- core 3: exact clusters, then the chains weighted to match them, then the control
( exact $PC 1000 100000 2400 11; exact 0.55 500 20000 600 12; exact 0.5 250 100000 600 13
  tilt 0.5 250 0 120 201; tilt 0.5 250 1 120 202; tilt 0.5 250 2 450 502; tilt 0.5 250 3 450 503; tilt 0.5 250 4 450 504
  tilt 0.55 500 0 150 203; tilt 0.55 500 1 150 204
  tilt $PC 1000 0 300 205; tilt $PC 1000 1 300 206
  # positive control: the chain at a slightly wrong density, compared with the exact clusters at 0.5
  $M -load $D/x0.5_250_0000.bin -p 0.51 -pcp 0.4 -moves 1e18 -tsec 120 -every 250 -seed 207 -out $D/t0.51_250_0.ts -feat $D/t0.51_250_0.feat 2> $D/t0.51_250_0.log
  uni 16000 A 1500 103 A; uni 16000 B 1500 104 B
  # more exact clusters at 250 cells, and a high-precision comparison at 60 cells
  for k in 2 3 4; do ./leath -n 250 -p 0.5 -count 10000000 -tsec 840 -seed $((40+k)) -out $D/x0.5_250_s$k.feat 2> $D/x0.5_250_s$k.log; done
  ./leath -n 60 -p 0.5 -count 100000000 -tsec 240 -seed 71 -out $D/x0.5_60.feat -snap $D/x0.5_60 -nsnap 2 2> $D/x0.5_60.log
  $M -load $D/x0.5_60_0000.bin -p 0.5 -pcp 0.4 -moves 1e18 -tsec 420 -every 240 -seed 721 -out /dev/null -feat $D/t0.5_60_0.feat 2> $D/t0.5_60_0.log
  ./leath -n 250 -p 0.5 -count 8 -seed 77 -out /dev/null -snap $D/y0.5_250 -nsnap 8
  for i in 0 1 2 3 4 5 6 7; do $M -load $D/y0.5_250_000$i.bin -p 0.5 -pcp 0.4 -moves 1e18 -tsec 240 -every 1000 -seed $((950+i)) -out /dev/null -feat $D/r0.5_250_$i.feat 2> $D/r0.5_250_$i.log; done
  for i in 1 2 3 4 5 6 7 8; do $M -load prod/ckpt_A_n250.bin -pcp 0.4 -moves 1e18 -tsec 240 -every 1000 -seed $((900+i)) -out /dev/null -feat $D/m250_$i.feat 2> $D/m250_$i.log; done ) &

# ---- core 4: enumeration, short uniform chains, grown populations
( for n in 14 16; do ./enumerate $n > $D/exact_$n.txt; $M -n $n -pcp 0.4 -moves 4.001e8 -burn 100000 -every 200 -seed $((300+n)) -out /dev/null -feat $D/e$n.feat 2> $D/e$n.log; done
  for p in 0.5 $PC; do ./enumerate 14 $p > $D/exact_14_p$p.txt; $M -n 14 -p $p -pcp 0.4 -moves 4.001e8 -burn 100000 -every 200 -seed 314 -out /dev/null -feat $D/e14_p$p.feat 2> $D/e14_p$p.log; done
  uni 4000 A 400 105 A; uni 4000 B 400 106 B; uni 1000 A 200 107 A; uni 1000 B 200 108 B; uni 250 A 100 109 A; uni 250 B 100 110 B
  # 500 and 2,000 cells have only ladder C; two chains from its final shape, first 5% dropped in the analysis
  uni 500 A 300 111 C; uni 500 B 300 112 C; uni 2000 A 300 113 C; uni 2000 B 300 114 C
  for k in 1 2 3 4 5 6 7 8; do grow 1000 200000 250,500,1000 $k gk1000_$k; done
  for k in 1 2 3 4 5 6; do grow 2000 100000 2000 $((20+k)) gk2000_$k; done
  # smaller populations, to show how the growth sampler degrades
  for k in $(seq 1 24); do grow 1000 12500 1000 $((70+k)) gs1000a_$k; done
  for k in $(seq 1 12); do grow 1000 50000 1000 $((100+k)) gs1000b_$k; done
  for k in $(seq 1 16); do grow 2000 25000 2000 $((40+k)) gs2000a_$k; done ) &
wait
python3 analysis/convergence.py
