#!/bin/bash
set -e

echo "Usage: ./run_query_weak.sh [basename] [args for alias]"

basename=$1

filename=/global_data/lorenz/wrs/w${basename}results_$(date "+%F.%H-%M-%S").txt
linkname=/global_data/lorenz/wrs/w${basename}results.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

for threads in {1,2,4,8,12,16,24,32,40,52,64,80,104,128,158}
do
    #reps=$(echo "80*sqrt($threads)" | bc)
    qsize=$((10**6*threads))
    ./benchmark/alias -t ${threads} -r 1000 -C -q ${qsize} ${@:2} | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
