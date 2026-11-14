#!/bin/bash
# (Re)launch every production job.  All of them resume from their checkpoints.
cd $(dirname $0)/..
pgrep -f "manyomino -n" > /dev/null && { echo "already running"; exit 0; }
L() { setsid nohup "$@" > /dev/null 2>&1 & }
L scripts/ladder.sh A 101 250 "40 80 300 800 1800 2700 5400"
L scripts/ladder.sh B 202 250 "40 80 300 800 1800 2700 5400"
L scripts/ladder.sh C 303 125 "40 60 160 500 1200 2400 4800"
L scripts/baseline.sh
# head-to-head from a filled square at n = 30,000: relocate-only vs the full chain
L scripts/stage.sh prod2 oldsq 30000 1500 9 "-shape rect" "-pcp 0 -relK 0 -every 200000 -full 0"
L scripts/stage.sh prod2 newsq 30000 400 9 "-shape rect" "-pcp 0.4 -every 100000 -full 0"
