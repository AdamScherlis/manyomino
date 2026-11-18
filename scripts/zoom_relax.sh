#!/bin/bash
# How a doubled shape settles.  Equilibrium samples of 16,000 and 64,000 cells (snapshots of the
# production ladder in prod/) are doubled and relaxed by src/zoom.cpp for 1,024 sweeps, and the
# shape is saved after 1, 2, 4, ..., 1,024 sweeps.  analysis/zoom_relax.py compares those shapes
# with sampled ones.  About 20 core-minutes.
cd $(dirname $0)/..
mkdir -p zoomtest/relax
[ -x zoom ] || g++ -O3 -march=native -std=c++17 -o zoom src/zoom.cpp
k=0
for f in $(ls prod/snap_A_n16000_p*_000[13579]*.bin prod/snap_B_n16000_p*_000[13579]*.bin | head -12) \
         $(ls prod/snap_A_n64000_p* | awk 'NR%4==0' | head -6) $(ls prod/snap_B_n64000_p* | awk 'NR%4==0' | head -6); do
  k=$((k + 1)); t=$(basename $f .bin | sed 's/snap_//')
  echo "-load $f -inflate 2 -sweeps 1024 -seed $((100 + k)) -snaps zoomtest/relax/$t -snapAt '1 2 4 8 16 32 64 128 256 512 1024'"
done | xargs -P 4 -L 1 sh -c './zoom "$@" > /dev/null 2>&1' zoom
