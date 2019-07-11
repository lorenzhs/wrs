#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

# DEGS sequential
# M and L are redundant (multi-alias and multi-sweep), run only multi-sweep
$DIR/run_query_strong.sh quni -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455
$DIR/run_query_strong.sh qp1  -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 1
#$DIR/run_query_strong.sh qp2  -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 2
#$DIR/run_query_strong.sh qp05 -D -E -G -S -M -n $((10**9)) -i 3 -s 11235813213455 -p 0.5
