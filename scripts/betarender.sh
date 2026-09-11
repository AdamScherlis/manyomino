#!/bin/bash
cd "$(dirname "$0")/.."
M=rust/target/release/manyomino
SNAP=gallery/raw/n100000_rect_40000000.txt
R2=476600
if [ "$1" = lane1 ]; then BS="30 30000"; else BS="-3 1000"; fi
for b in $BS; do
  beta=$(python3 -c "print($b/$R2)")
  tag=$(echo $b | tr '.-' 'pm')
  $M run --init-file $SNAP --steps 120000000 --record-every 100000 --seed 15$RANDOM --cp-inv 3 --pv-inv 2 --cp-cap 0 --beta $beta --out results/beta/n100000_b$tag.csv --dump-every 30000000 --dump-prefix gallery/raw/tilt_b$tag >/dev/null 2>results/beta/n100000_b$tag.log
  echo "done b=$b"
done
