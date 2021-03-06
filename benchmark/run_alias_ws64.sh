#!/bin/bash
set -e

echo "Usage: ./run_alias_ws64.sh size-per-thread [other args for alias]"

size=$1

filename=/global_data/lorenz/wrs/wsresults64_$(date "+%F.%H-%M-%S").txt
linkname=/global_data/lorenz/wrs/wsresults64.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

# Sequential baseline for speedup calculation
# DEGS sequential, BLMORP parallel.
./benchmark/alias -B -L -M -O -R -P -D -n $size -t 1 ${@:2} | tee -a $filename
for threads in {1,2,4,6,8,12,16,20,24,32,40,48,56,62}
do
    # O and R are redundant (only affects query)
    ./benchmark/alias -D -E -G -S -R -n $((size*threads)) -t $threads ${@:2} | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
