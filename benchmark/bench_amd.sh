#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

$DIR/bench_query_64.sh
$DIR/bench_weak_64.sh
$DIR/bench_strong_64.sh
