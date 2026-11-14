#!/bin/bash
# Size ladder: equilibrate at N0, then repeatedly inflate x2 (4x the cells) and re-equilibrate.
# Every stage is a production run in its own right.  Safe to re-run after an interruption.
#   ladder.sh TAG SEED N0 "T0 T1 T2 ..." [OUTDIR]      (seconds of wall-clock per stage)
TAG=$1; SEED=$2; N=$3; TIMES=($4); OUT=${5:-prod}
D=$(dirname $0); init="-shape rect"
for T in "${TIMES[@]}"; do
  snapEvery=$((1500*N)); [ $N -ge 60000 ] && snapEvery=$((150*N))
  $D/stage.sh $OUT $TAG $N $T $SEED "$init" "-pcp 0.4 -every $((3*N)) -snapEvery $snapEvery"
  init="-load $OUT/ckpt_${TAG}_n${N}.bin -inflate 2"
  N=$((4*N)); SEED=$((SEED+1000))
done
