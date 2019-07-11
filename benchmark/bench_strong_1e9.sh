#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# no queries
$DIR/run_alias.sh 9 -Q -i 10 -r 20 -s 11235813213455
