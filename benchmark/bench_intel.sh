#!/bin/bash
set -e

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

$DIR/bench_query.sh
$DIR/bench_weak.sh
$DIR/bench_strong.sh
$DIR/bench_query_strong.sh
$DIR/bench_query_weak.sh
