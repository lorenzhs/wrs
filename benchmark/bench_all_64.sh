#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Construction benchmarks
$DIR/bench_strong_64.sh
$DIR/bench_strong_1e9_64.sh
$DIR/bench_weak64.sh

# Query benchmarks
$DIR/bench_query_64.sh

# Sequential comparison
$DIR/bench_seq.sh
