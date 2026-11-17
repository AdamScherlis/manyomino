#!/bin/bash
# The clump recording of scripts/clumps.sh at 1,024,000 cells: four chains continued from the
# ends of the production runs (prod/ckpt_{A,B,D,E}_n1024000.bin), sampled once per move per cell.
#   scripts/clumps_hero.sh T [A]    T seconds per chain; A as in clumps.sh (output in the same
#                                   directory, so analysis/clumps.py picks the size up)
cd $(dirname $0)/..
T=${1:-2700}; A=${2:-2}; OUT=clump; EX="-pcp 0.4"
[ "$A" != 2 ] && { OUT=clump_a$A; EX="$EX -relA $A"; }
mkdir -p $OUT
[ -x manyomino ] || g++ -O3 -march=native -std=c++17 -o manyomino src/manyomino.cpp
L() { setsid nohup "$@" > /dev/null 2>&1 & }
k=0
for tag in A B D E; do
  k=$((k + 1))
  BLOBS=1 L scripts/stage.sh $OUT $tag 1024000 $T $((8500 + 100 * k)) "-load prod/ckpt_${tag}_n1024000.bin" "$EX -every 1024000"
done
