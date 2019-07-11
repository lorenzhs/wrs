#!/bin/bash
set -e

echo "Usage: ./run_alias_query.sh [basename] [args for alias]"

basename=$1

filename=/global_data/lorenz/wrs/${basename}results_$(date "+%F.%H-%M-%S").txt
linkname=/global_data/lorenz/wrs/${basename}results.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

for factor in {100,316,1000,3162,10000,31620,100000}
do
    reps=$((10**6/factor))
    ./benchmark/alias -t 158 -C -q $((factor*10**4)) ${@:2} -r ${reps} | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
