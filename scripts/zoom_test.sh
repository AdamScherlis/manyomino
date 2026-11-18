#!/bin/bash
# Is a refined shape distinguishable from a sampled one?  32 equilibrium samples of 13,650 cells
# (four chains from a square, 3,000 moves per cell discarded, then a sample every 250) are each
# refined four times over by scripts/zoom_chain.sh, to about 1,030,000 cells, with T sweeps per
# level.  analysis/zoom_test.py sets their observables against the equilibrium runs at 1,024,000
# cells.  About 40 core-minutes for T = 64, in proportion for other T.
#   scripts/zoom_test.sh [T ...]        default: 64 16 32 256
cd $(dirname $0)/..
mkdir -p zoomtest; N0=13650
[ -x zoom ] || g++ -O3 -march=native -std=c++17 -o zoom src/zoom.cpp
[ -x manyomino ] || g++ -O3 -march=native -std=c++17 -o manyomino src/manyomino.cpp
for c in 0 1 2 3; do
  [ -f zoomtest/base_${c}_0007.bin ] && continue
  (./manyomino -n $N0 -shape rect -pcp 0.4 -seed $((9000 + c)) -moves $((3000 * N0)) -every $((3000 * N0)) -out zoomtest/burn_$c.txt -snap zoomtest/burn_$c.bin 2> /dev/null &&
   ./manyomino -n $N0 -load zoomtest/burn_$c.bin -pcp 0.4 -seed $((9100 + c)) -moves $((2000 * N0)) -every $((2000 * N0)) -out zoomtest/base_$c.txt \
     -snapEvery $((250 * N0)) -snap zoomtest/base_$c 2> /dev/null) &
done; wait
# the starting samples' own observables, and equilibrium at their size to set them against
for f in zoomtest/base_?_000?.bin; do
  [ -f ${f%.bin}.feat ] || ./manyomino -n $N0 -load $f -pcp 0 -moves 1 -every 1 -feat ${f%.bin}.feat -out ${f%.bin}.txt 2> /dev/null
done
mkdir -p zoomtest/ref
for c in 0 1 2 3; do
  [ -f zoomtest/ref/r$c.feat ] || ./manyomino -n $N0 -load zoomtest/base_${c}_0007.bin -pcp 0.4 -seed $((9300 + c)) -moves 1e15 -tsec 300 -every $N0 \
    -feat zoomtest/ref/r$c.feat -out zoomtest/ref/r$c.txt 2> /dev/null &
done; wait
for T in ${@:-64 16 32 256}; do for c in 0 1 2 3; do for i in 0 1 2 3 4 5 6 7; do
  echo "zoomtest/base_${c}_000$i.bin zoomtest/z${T}_${c}_$i $T 4 $((1000 * T + 100 * c + 10 * i))"
done; done; done | xargs -P 4 -L 1 scripts/zoom_chain.sh
