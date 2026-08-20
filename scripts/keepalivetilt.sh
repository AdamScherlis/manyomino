#!/bin/bash
# Relaunch the four tilted n=100k render lanes and the n=10000 beta=0
# distribution run if dead, resuming each from its freshest dump with
# cumulative step accounting (--step-offset).  Safe to run repeatedly.
cd "$(dirname "$0")/.."
M=rust/target/release/manyomino
GEN=$(date +%s)
R2=476600
TOTAL=120000000

for b in 30 30000 -3 1000; do
  tag=$(echo $b | tr '.-' 'pm')
  if pgrep -f "manyomino run.*tilt_b$tag" >/dev/null; then continue; fi
  beta=$(python3 -c "print($b/$R2)")
  # freshest cumulative dump for this lane, else the equilibrium seed
  latest=$(ls -t gallery/raw/tilt_b${tag}_*.txt 2>/dev/null | head -1)
  if [ -n "$latest" ]; then
    base=$(basename "$latest" .txt); base=${base##*_}
    init="$latest"
  else
    base=0
    init=gallery/raw/n100000_rect_40000000.txt
  fi
  rem=$((TOTAL - base))
  if [ "$rem" -le 0 ]; then continue; fi
  nohup $M run --init-file "$init" --steps $rem --step-offset $base \
    --record-every 100000 --seed 17$((GEN % 10000)) --cp-inv 3 --pv-inv 2 \
    --cp-cap 0 --beta $beta --out results/beta/n100000_b${tag}_g$GEN.csv \
    --dump-every 5000000 --dump-prefix gallery/raw/tilt_b$tag \
    >/dev/null 2>results/beta/n100000_b${tag}_g$GEN.log &
  disown
  GEN=$((GEN + 1))
done

# n=10000 beta=0 distribution run (100M steps total)
if ! pgrep -f "manyomino run.*n10000_b0" >/dev/null; then
  latest=$(ls -t gallery/raw/n10000b0_*.txt 2>/dev/null | head -1)
  if [ -n "$latest" ]; then
    base=$(basename "$latest" .txt); base=${base##*_}
    init="$latest"
  else
    base=0
    init=gallery/raw/n10000_rect_600000000.txt
  fi
  rem=$((100000000 - base))
  if [ "$rem" -gt 0 ]; then
    nohup nice -n 5 $M run --init-file "$init" --steps $rem --step-offset $base \
      --record-every 5000 --seed 18$((GEN % 10000)) --cp-inv 3 --pv-inv 2 \
      --cp-cap 0 --beta 0 --out results/beta/n10000_b0_g$GEN.csv \
      --dump-every 10000000 --dump-prefix gallery/raw/n10000b0 \
      >/dev/null 2>results/beta/n10000_b0_g$GEN.log &
    disown
  fi
fi
sleep 1
echo "tilt lanes: $(pgrep -fc 'manyomino ru[n].*beta') total: $(pgrep -fc 'manyomino ru[n]')"
