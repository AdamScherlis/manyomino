#!/bin/bash
# How long do loop-connected clumps of a given size last?  Equilibrium chains at three sizes,
# started from the equilibrated checkpoints of the production ladder, recording at every move
# per cell how many clumps of at least each size the shape holds (-blobs).  analysis/clumps.py
# turns that into the size distribution and the persistence time of clumps against their size.
#   scripts/clumps.sh T [A]     T seconds per chain on four cores; A is the tail exponent of the
#                               relocate move's connectivity radius (-relA; 2, the rule of the
#                               published runs, writes to clump/, any other value to clump_aA/)
# Resumable: run it again after an interruption.
cd $(dirname $0)/..
T=${1:-10800}; A=${2:-2}; OUT=clump; EX="-pcp 0.4"
EX="$EX -relA $A"; [ "$A" != 2 ] && OUT=clump_a$A
mkdir -p $OUT
[ -x manyomino ] || g++ -O3 -march=native -std=c++17 -o manyomino src/manyomino.cpp
L() { setsid nohup "$@" > /dev/null 2>&1 & }
# with another exponent, two of the chains also record the full observable library
for job in "A 4000 8101 A" "A 16000 8201 A" "B 16000 8302 B $([ "$A" != 2 ] && echo 1)" "A 64000 8401 A $([ "$A" != 2 ] && echo 1)"; do
  set -- $job
  [ -f $OUT/start_$1_n$2.bin ] || cp prod/ckpt_$4_n$2.bin $OUT/start_$1_n$2.bin
  FEAT=${5:-} BLOBS=1 L scripts/stage.sh $OUT $1 $2 $T $3 "-load $OUT/start_$1_n$2.bin" "$EX -every $2"
done
