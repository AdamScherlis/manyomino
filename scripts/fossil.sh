#!/bin/bash
# Runs behind analysis/fossil.py: from an equilibrated shape, save the shape every 5 moves per
# cell for 600 moves per cell, several seeds per size.  About an hour on one core.
cd "$(dirname "$0")/.."; mkdir -p surv
run() {   # n tag seed
  cp prod/ckpt_A_n$1.bin surv/s$1$2_start.bin; echo 5 > surv/s$1$2.step
  ./manyomino -load prod/ckpt_A_n$1.bin -pcp 0.4 -moves $(($1*600)) -every $(($1*600)) -snapEvery $(($1*5)) -snap surv/s$1$2 -seed $3 -out /dev/null 2> surv/s$1$2.log
}
run 16000 "" 77; run 64000 "" 77; run 256000 "" 77
i=0; for t in b c d e f g h; do i=$((i+1)); run 4000 $t $((80+i)); run 16000 $t $((80+i)); done
i=0; for t in b c d; do i=$((i+1)); run 64000 $t $((90+i)); done
