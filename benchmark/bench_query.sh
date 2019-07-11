#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

# DEGS sequential
# L and M are redundant (multi-sweep and multi-alias), run only multi-sweep
$DIR/run_alias_query.sh quni -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455
$DIR/run_alias_query.sh qp1  -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 1
$DIR/run_alias_query.sh qp2  -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 2
$DIR/run_alias_query.sh qp05 -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 0.5
