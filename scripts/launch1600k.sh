#!/bin/bash
# Inflate the freshest converged 100k snapshots (one per lineage) and launch
# the n=1.6M pair.  Run from anywhere.
set -e
cd "$(dirname "$0")/.."
BAR_SRC=$(ls -t gallery/raw/n400000*bar*.txt | head -1)
RECT_SRC=$(ls -t gallery/raw/n400000*rect*.txt | head -1)
echo "inflating $BAR_SRC and $RECT_SRC"
python3 python/inflate.py "$BAR_SRC" gallery/raw/n1600000_seedA.txt
python3 python/inflate.py "$RECT_SRC" gallery/raw/n1600000_seedB.txt
cd rust
M=./target/release/manyomino
nohup $M run --init-file ../gallery/raw/n1600000_seedA.txt --steps 300000000 --record-every 100000 --seed 901 --cp-inv 3 --cp-cap 4000 --out ../results/prod/n1600000_bar.csv --dump-every 30000000 --dump-prefix ../gallery/raw/n1600000_bar --check-every 200000000 >/dev/null 2>../results/prod/n1600000_bar.log &
disown
nohup $M run --init-file ../gallery/raw/n1600000_seedB.txt --steps 300000000 --record-every 100000 --seed 902 --cp-inv 3 --cp-cap 4000 --out ../results/prod/n1600000_rect.csv --dump-every 30000000 --dump-prefix ../gallery/raw/n1600000_rect --check-every 200000000 >/dev/null 2>../results/prod/n1600000_rect.log &
disown
echo "1.6M pair launched"
