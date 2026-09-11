#!/bin/bash
# Deep negative-b extension of the tilted sweep at n=1000, annealed:
# each run seeds from the previous (less negative) run's final state so
# the stretch transient is short.  Also two n=3000 universality checks.
cd "$(dirname "$0")/.."
M=rust/target/release/manyomino
R2=1379
PREV=gallery/raw/n1000_rect_200000000.txt
for b in -5 -8 -12 -20 -30; do
  beta=$(python3 -c "print($b/$R2)")
  tag=$(echo $b | tr '.-' 'pm')
  $M run --init-file $PREV --steps 25000000 --record-every 2000 --seed 19$RANDOM \
    --cp-inv 3 --pv-inv 2 --cp-cap 0 --beta $beta \
    --out results/beta/n1000_b$tag.csv --dump-final gallery/raw/nb_n1000_b$tag.txt \
    >/dev/null 2>results/beta/n1000_b$tag.log
  PREV=gallery/raw/nb_n1000_b$tag.txt
  echo "done n=1000 b=$b"
done
R2=5385
PREV=gallery/raw/n3000_rect_150000000.txt
for b in -5 -10; do
  beta=$(python3 -c "print($b/$R2)")
  tag=$(echo $b | tr '.-' 'pm')
  $M run --init-file $PREV --steps 36000000 --record-every 2000 --seed 19$RANDOM \
    --cp-inv 3 --pv-inv 2 --cp-cap 0 --beta $beta \
    --out results/beta/n3000_b$tag.csv --dump-final gallery/raw/nb_n3000_b$tag.txt \
    >/dev/null 2>results/beta/n3000_b$tag.log
  PREV=gallery/raw/nb_n3000_b$tag.txt
  echo "done n=3000 b=$b"
done
