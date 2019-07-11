#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# run sequential algorithms with n=10^8
$DIR/run_seq.sh -n $((10**8)) -q $((10**8)) -i 10 -r 10 -s 11235813213455
