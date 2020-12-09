#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# no queries
$DIR/run_alias_64.sh 8 -Q -i 20 -r 25 -s 11235813213455
# Use fewer repetitions for larger run
$DIR/run_alias_64.sh 9 -Q -i 10 -r 15 -s 11235813213455
