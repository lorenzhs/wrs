#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

# M and L are redundant (multi-alias and multi-sweep), run only multi-sweep
# B and P are reduntant (psa and psa2), run only psa2
$DIR/run_query_weak.sh quni -S -E -B -D -G -M -n $((10**9)) -i 3 -s 11235813213455
$DIR/run_query_weak.sh qp1  -S -E -B -D -G -M -n $((10**9)) -i 3 -s 11235813213455 -p 1
#$DIR/run_query_weak.sh qp2  -S -E -B -D -G -M -n $((10**9)) -i 3 -s 11235813213455 -p 2
#$DIR/run_query_weak.sh qp05 -S -E -B -D -G -M -n $((10**9)) -i 3 -s 11235813213455 -p 0.5
