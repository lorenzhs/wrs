#!/bin/bash
set -e

echo "Usage: ./run_seq.sh [args for alias]"

host=$(hostname)
filename=out/seqresults_${host}_$(date "+%F.%H-%M-%S").txt
linkname=out/seqresults_${host}.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on ${host} on $(date)" | tee -a $filename

# Sequential benchmark
# DEGS sequential, BLMOR parallel
./benchmark/alias -B -L -M -O -R -t 1 $@ | tee -a $filename

rm -f $linkname
ln -s $filename $linkname
