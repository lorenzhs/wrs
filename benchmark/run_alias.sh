#!/bin/bash
set -e

echo "Usage: ./run_alias.sh log10(n) [args for alias]"

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

logsize=$1
filename=out/results_${logsize}_$(date "+%F.%H-%M-%S").txt
linkname=out/results_${logsize}.txt

size=$((10**logsize))

# DEGS sequential, BLMOR parallel
./benchmark/alias -n ${size} ${@:2} -B -L -M -O -R -D -t 1 | tee -a $filename
for threads in {1,2,4,8,12,16,24,32,40,52,64,80,104,128,158}
do
    # O and R are redundant (only affects query)
    ./benchmark/alias -n ${size} ${@:2} -D -E -G -S -R -t $threads | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
