#!/bin/bash
# Run one production stage for T seconds of wall-clock in total, resumably: each invocation
# writes a new part file and restarts from the rolling checkpoint of the previous part.
#   stage.sh OUTDIR TAG N T SEED INIT "EXTRA ARGS"
# INIT is either "-shape rect" or "-load file -inflate 2" and is used for the first part only.
OUT=$1; TAG=$2; N=$3; T=$4; SEED=$5; INIT=$6; EXTRA=$7
BIN=$(dirname $0)/../manyomino
base=$OUT/ts_${TAG}_n${N}; ck=$OUT/ckpt_${TAG}_n${N}.bin
while true; do
  done_s=0; parts=0
  for f in ${base}.part*.txt; do [ -f "$f" ] || continue; parts=$((parts+1)); s=$(tail -n 1 "$f" | awk '$1 !~ /#/ {print int($2)}'); done_s=$((done_s + ${s:-0})); done
  left=$((T - done_s)); [ $left -le 5 ] && break
  init="$INIT"; [ -f $ck ] && init="-load $ck"
  $BIN -n $N -seed $((SEED + parts)) -moves 1e15 -tsec $left $init $EXTRA -ckpt $ck -ckptSec 90 \
     -snap $OUT/snap_${TAG}_n${N}_p${parts} -out ${base}.part${parts}.txt 2>> $OUT/log_${TAG}_n${N}.txt
done
