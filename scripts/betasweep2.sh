#!/bin/bash
cd "$(dirname "$0")/.."
M=rust/target/release/manyomino
if [ "$1" = n1000 ]; then
  SNAP=gallery/raw/n1000_rect_200000000.txt; R2=1379; N=1000; STEPS=25000000
  BS="-3 -1 -0.3 0.3 1 3 10 30 100 300 1000"
else
  SNAP=gallery/raw/n3000_rect_150000000.txt; R2=5385; N=3000; STEPS=30000000
  BS="10 30 100 300 3000"
fi
for b in $BS; do
  beta=$(python3 -c "print($b/$R2)")
  tag=$(echo $b | tr '.-' 'pm')
  $M run --init-file $SNAP --steps $STEPS --record-every 2000 --seed 14$RANDOM --cp-inv 3 --pv-inv 2 --cp-cap 0 --beta $beta --out results/beta/n${N}_b$tag.csv >/dev/null 2>results/beta/n${N}_b$tag.log
  echo "done n=$N b=$b"
done
