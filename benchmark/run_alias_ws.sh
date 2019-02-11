#!/bin/bash
set -e

echo "Usage: ./run_alias_ws.sh size-per-thread [other args for alias]"

size=$1

filename=wsresults_$(date "+%F.%H-%M-%S").txt
linkname=wsresults.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $HOST on $(date)" | tee -a $filename

# Sequential baseline for speedup calculation
./benchmark/alias -M -P -n $size -t $threads 1 | tee -a $filename
for threads in {1,2,3,4,5,6,8,10,12,16,20,24,28,32,40,48,56,60,64,72,80,100,120,140,160}
do
    ./benchmark/alias -S -n $((size*threads)) -t $threads ${@:2} | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
