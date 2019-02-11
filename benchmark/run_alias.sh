#!/bin/bash
set -e

echo "Usage: ./run_alias.sh [args for alias]"

filename=results_$(date "+%F.%H-%M-%S").txt
linkname=results.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $HOST on $(date)" | tee -a $filename

./benchmark/alias $@ -P -M -t 1 | tee -a $filename
for threads in {1,2,3,4,5,6,8,10,12,16,20,24,28,32,40,48,56,60,64,72,80,100,120,140,160}
do
    ./benchmark/alias $@ -S -t $threads | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
