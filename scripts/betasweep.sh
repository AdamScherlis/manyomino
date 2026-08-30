#!/bin/bash
# beta sweep at n=3000 from an equilibrated uniform-ensemble snapshot.
# usage: betasweep.sh <lane: neg|pos>
cd "$(dirname "$0")/.."
M=rust/target/release/manyomino
SNAP=gallery/raw/n3000_rect_150000000.txt
R2=5385   # <Rg^2>_0 at n=3000
if [ "$1" = neg ]; then BS="-2 -1.5 -1 -0.6 -0.3 -0.15 0"; else BS="0.15 0.3 0.6 1 1.5 2 3 5"; fi
for b in $BS; do
  beta=$(python3 -c "print($b/$R2)")
  tag=$(echo $b | tr '.-' 'pm')
  $M run --init-file $SNAP --steps 36000000 --record-every 2000 --seed 13$RANDOM --cp-inv 3 --pv-inv 2 --cp-cap 0 --beta $beta --out results/beta/n3000_b$tag.csv >/dev/null 2>results/beta/n3000_b$tag.log
  echo "done b=$b"
done
