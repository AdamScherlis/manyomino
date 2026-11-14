#!/bin/bash
# Uncontended speed of the production chain at each size, from equilibrated checkpoints.
# Run on an idle machine.  Writes analysis/bench.json  {n: moves per second}.
cd $(dirname $0)/..
echo "{" > analysis/bench.json; sep=""
for f in $(ls prod/ckpt_A_n*.bin | sort -t n -k 2 -n); do
  n=$(echo $f | sed 's/.*_n\([0-9]*\)\.bin/\1/')
  r=$(./manyomino -load $f -n $n -pcp 0.4 -seed 77 -moves 1e15 -tsec ${1:-20} -every $((3*n)) -full 0 -out /dev/null 2>&1 | sed -n 's/.*(\([0-9.]*\) Mmoves\/s).*/\1/p')
  echo "$sep\"$n\": $(python3 -c "print($r*1e6)")" >> analysis/bench.json; sep=","
done
echo "}" >> analysis/bench.json; cat analysis/bench.json
