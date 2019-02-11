#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# run queries as well
$DIR/run_alias.sh -n $((10**9)) -q $((10**7)) -i 10 -r 20 -s 11235813213455
