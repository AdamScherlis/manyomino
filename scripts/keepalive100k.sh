#!/bin/bash
# Relaunch the four n=100k chains from the freshest lineage snapshots if
# they are not running (container restarts kill them).  Generation counter
# in the seed keeps RNG streams distinct.
cd "$(dirname "$0")/.."
RUNNING=$(pgrep -fc "manyomino run" || true)
if [ "$RUNNING" -ge 4 ]; then echo "chains alive ($RUNNING)"; exit 0; fi
pkill -f "manyomino run" 2>/dev/null; sleep 1
GEN=$(date +%s)
M=rust/target/release/manyomino
A=$(ls -t gallery/raw/n100000_bar_*.txt | head -1)
B=$(ls -t gallery/raw/n100000_rect_*.txt | head -1)
C=$(ls -t gallery/raw/n100000c_bar_*.txt | head -1)
D=$(ls -t gallery/raw/n100000c_rect_*.txt | head -1)
i=0
for SRC in "$A" "$B" "$C" "$D"; do
  i=$((i+1))
  case $i in
    1) TAG=bar; PFX=gallery/raw/n100000_bar;;
    2) TAG=rect; PFX=gallery/raw/n100000_rect;;
    3) TAG=bar_c; PFX=gallery/raw/n100000c_bar;;
    4) TAG=rect_c; PFX=gallery/raw/n100000c_rect;;
  esac
  OUT="results/prod/n100000_${TAG}_g${GEN}.csv"
  nohup $M run --init-file "$SRC" --steps 1000000000 --record-every 50000 --seed $((GEN+i)) --cp-inv 3 --cp-cap 0 --out "$OUT" --dump-every 20000000 --dump-prefix "$PFX" --check-every 250000000 >/dev/null 2>"${OUT%.csv}.log" &
  disown
done
sleep 1
echo "relaunched gen $GEN: $(pgrep -fc 'manyomino run') chains"
