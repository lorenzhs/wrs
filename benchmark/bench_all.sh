#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Construction benchmarks
$DIR/bench_strong.sh
$DIR/bench_strong_1e9.sh
$DIR/bench_weak.sh

# Query benchmarks
$DIR/bench_query.sh
$DIR/bench_query_strong.sh
$DIR/bench_query_weak.sh

# Sequential comparison
$DIR/bench_seq.sh
