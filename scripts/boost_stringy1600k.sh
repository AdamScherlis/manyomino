#!/bin/bash
# The n=400k rung is published; retire its chains and give the cores to a
# second stringy-lineage n=1.6M chain (independent seed, poolable segment)
# to accelerate the stringy ESS toward the publication gate.
cd "$(dirname "$0")/.."
pkill -f "manyomino run.*n400000"
sleep 1
M=rust/target/release/manyomino
GEN=$(date +%s)
if ! pgrep -f "manyomino run.*n1600000_barB" >/dev/null; then
  latest=$(ls -t gallery/raw/n1600000_barB_*.txt 2>/dev/null | head -1)
  [ -z "$latest" ] && latest=$(ls -t gallery/raw/n1600000_bar_*.txt | head -1)
  nohup $M run --init-file "$latest" --steps 300000000 --record-every 100000 \
    --seed $((GEN + 7)) --cp-inv 3 --pv-inv 2 --cp-cap 0 \
    --out results/prod/n1600000_barB_g$GEN.csv \
    --dump-every 1000000 --dump-prefix gallery/raw/n1600000_barB \
    >/dev/null 2>results/prod/n1600000_barB_g$GEN.log &
  disown
fi
sleep 1
ls -t gallery/raw/n1600000_barB_*.txt 2>/dev/null | tail -n +4 | xargs -r rm -f
echo "1.6M chains: $(pgrep -fc 'manyomino ru[n].*n1600000')"
