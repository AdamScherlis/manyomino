#!/bin/bash
# Refine a shape LEVELS times (src/zoom.cpp, T sweeps per level) and write the observables of the
# result (PREFIX.feat) and the shape itself (PREFIX.bin).  Nothing else is done to it: the number
# of cells is whatever the refinement arrives at.
#   scripts/zoom_chain.sh IN PREFIX T LEVELS SEED
cd $(dirname $0)/..
IN=$1; P=$2; T=$3; L=$4; SEED=$5
cur=$IN
for k in $(seq 1 $L); do
  ./zoom -load $cur -inflate 2 -sweeps $T -seed $((SEED + k)) -out ${P}_l$k.bin 2>> $P.log > /dev/null
  [ $k -gt 1 ] && rm -f ${P}_l$((k - 1)).bin
  cur=${P}_l$k.bin
done
mv $cur $P.bin
N=$(od -An -t d4 -N 4 $P.bin | tr -d ' ')
./manyomino -n $N -load $P.bin -pcp 0 -seed $SEED -moves 1 -every 1 -feat $P.feat -out $P.txt 2>> $P.log
