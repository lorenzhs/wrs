#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# no queries
$DIR/run_alias_ws.sh $((10**7)) -Q -i 10 -r 20 -s 11235813213455
