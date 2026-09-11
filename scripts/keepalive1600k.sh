#!/bin/bash
# Resume the two n=1.6M chains from their freshest dumps if dead.
cd "$(dirname "$0")/.."
RUNNING=$(pgrep -fc "manyomino run.*n1600000" || true)
if [ "$RUNNING" -ge 2 ]; then echo "1.6M chains alive ($RUNNING)"; exit 0; fi
pkill -f "manyomino run.*n1600000" 2>/dev/null; sleep 1
GEN=$(date +%s)
M=rust/target/release/manyomino
for TAG in bar rect; do
  latest=$(ls -t gallery/raw/n1600000_${TAG}_*.txt 2>/dev/null | head -1)
  if [ -z "$latest" ]; then
    latest=$(ls gallery/raw/n1600000_seed[AB].txt 2>/dev/null | { [ "$TAG" = bar ] && head -1 || tail -1; })
  fi
  [ -z "$latest" ] && continue
  nohup $M run --init-file "$latest" --steps 300000000 --record-every 100000 \
    --seed $((GEN + $([ "$TAG" = bar ] && echo 1 || echo 2))) \
    --cp-inv 3 --pv-inv 2 --cp-cap 0 \
    --out results/prod/n1600000_${TAG}_g$GEN.csv \
    --dump-every 1000000 --dump-prefix gallery/raw/n1600000_${TAG} \
    >/dev/null 2>results/prod/n1600000_${TAG}_g$GEN.log &
  disown
done
sleep 1
echo "1.6M: $(pgrep -fc 'manyomino ru[n].*n1600000') running"
# prune old dumps (keep newest 3 per lineage; each is ~20MB)
for TAG in bar rect; do
  ls -t gallery/raw/n1600000_${TAG}_*.txt 2>/dev/null | tail -n +4 | xargs -r rm -f
done
