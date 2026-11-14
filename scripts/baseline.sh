#!/bin/bash
# The relocate-only chain (exact global connectivity test), for the mixing-time comparison.
D=$(dirname $0); OUT=${1:-prod}
for spec in "50 100" "100 300" "200 900" "400 1800" "800 2400" "1600 3000"; do
  set -- $spec
  $D/stage.sh $OUT reloc $1 $2 5 "-shape rect" "-pcp 0 -relK 0 -every $(( $1 * $1 / 20 + 20 ))"
done
