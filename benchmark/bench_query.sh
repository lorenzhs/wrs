#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

# L and M are redundant (multi-sweep, multi-alias), run only multi-sweep
# Same with B and P (PSA, PSA2), run only PSA2
$DIR/run_alias_query.sh quni -S -D -M -B -n $((10**9)) -i 3 -s 11235813213455
$DIR/run_alias_query.sh qp1  -S -D -M -B -n $((10**9)) -i 3 -s 11235813213455 -p 1
$DIR/run_alias_query.sh qp2  -S -D -M -B -n $((10**9)) -i 3 -s 11235813213455 -p 2
$DIR/run_alias_query.sh qp05 -S -D -M -B -n $((10**9)) -i 3 -s 11235813213455 -p 0.5
