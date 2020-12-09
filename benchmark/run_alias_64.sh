#!/bin/bash
set -e

echo "Usage: ./run_alias.sh log10(n) [other args for alias]"

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

logsize=$1
filename=/global_data/lorenz/wrs/results64_${logsize}_$(date "+%F.%H-%M-%S").txt
linkname=/global_data/lorenz/wrs/results64_${logsize}.txt

size=$((10**logsize))

# DEGS sequential, BLMORP parallel.
./benchmark/alias -n ${size} ${@:2} -B -L -M -O -R -P -D -t 1 | tee -a $filename
#for threads in {1,2,3,4,5,6,8,10,12,16,20,24,28,32,40,48,56,60,64,72,80,100,120,140,160}
for threads in {1,2,4,6,8,12,16,20,24,32,40,48,56,62}
do
    # O and R are redundant (only affects query)
    ./benchmark/alias -n ${size} ${@:2} -D -E -G -S -R -t $threads | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
