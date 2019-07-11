#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# no queries
$DIR/run_alias.sh 8 -Q -i 20 -r 25 -s 11235813213455
