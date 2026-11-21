#!/bin/bash
# Coarsen a shape LEVELS times (src/zoom.cpp -coarsen, T sweeps per level) and write the
# observables of the result (PREFIX.feat) and the shape itself (PREFIX.bin).
#   scripts/coarsen_chain.sh IN PREFIX T LEVELS SEED        (KEEP=1 keeps every level as PREFIX_lK.bin)
cd $(dirname $0)/..
IN=$1; P=$2; T=$3; L=$4; SEED=$5
cur=$IN
for k in $(seq 1 $L); do
  ./zoom -load $cur -coarsen -sweeps $T -seed $((SEED + k)) -out ${P}_l$k.bin 2>> $P.log > /dev/null
  [ $k -gt 1 ] && [ -z "${KEEP:-}" ] && rm -f ${P}_l$((k - 1)).bin
  cur=${P}_l$k.bin
done
if [ -n "${KEEP:-}" ]; then cp $cur $P.bin; else mv $cur $P.bin; fi
N=$(od -An -t d4 -N 4 $P.bin | tr -d ' ')
./manyomino -n $N -load $P.bin -pcp 0 -seed $SEED -moves 1 -every 1 -feat $P.feat -out $P.txt 2>> $P.log
