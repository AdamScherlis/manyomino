#!/bin/bash
# Ten times more running at the largest size.  Chains A and B continue from their checkpoints;
# D and E start from the final shapes of A and B with new seeds (their first 300 moves per cell
# are discarded as burn-in by analysis/summarize.py).  Resumable: run it again after an
# interruption.  About 7.5 hours on four cores.
cd $(dirname $0)/..
N=1024000; X=27000                      # seconds of new running per chain
[ -f prod/start_D_n$N.bin ] || cp prod/ckpt_A_n$N.bin prod/start_D_n$N.bin
[ -f prod/start_E_n$N.bin ] || cp prod/ckpt_B_n$N.bin prod/start_E_n$N.bin
EX="-pcp 0.4 -every $((3*N)) -snapEvery $((150*N))"
L() { FEAT=1 setsid nohup "$@" > /dev/null 2>&1 & }
L scripts/stage.sh prod A $N $((10840+X)) 6101 "" "$EX"          # 10,840 s were run before
L scripts/stage.sh prod B $N $((10840+X)) 6202 "" "$EX"
L scripts/stage.sh prod D $N $X 7101 "-load prod/start_D_n$N.bin" "$EX"
L scripts/stage.sh prod E $N $X 7202 "-load prod/start_E_n$N.bin" "$EX"
