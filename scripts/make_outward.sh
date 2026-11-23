#!/bin/sh
# Builds the zoom-out engine: page/outward.wasm for the browser (needs clang with the wasm32 target and
# wasm-ld), and ./outward_native for the tests of the engine itself.
set -e
cd "$(dirname "$0")/.."
clang --target=wasm32 -O2 -nostdlib -fno-builtin -Wl,--no-entry -Wl,--strip-all \
  -Wl,--export=init,--export=seed,--export=xs,--export=ys,--export=count,--export=capacity,--export=side \
  -Wl,--export=clear,--export=place,--export=hold,--export=pin,--export=run,--export=coarsen,--export=relax \
  -Wl,--export=pad,--export=grow,--export=link,--export=stat,--export=tube,--export=shadow,--export=partlabels,--export=parts,--export=orlabels,--export=add,--export=idle \
  -o page/outward.wasm src/outward.c
cc -O2 -DNATIVE -o outward_native src/outward.c
ls -l page/outward.wasm
