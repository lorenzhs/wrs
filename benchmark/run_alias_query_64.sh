#!/bin/bash
set -e

echo "Usage: ./run_alias_query64.sh [basename] [args for alias]"

basename=$1

filename=out/${basename}results134_$(date "+%F.%H-%M-%S").txt
linkname=out/${basename}results134.txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $(hostname) on $(date)" | tee -a $filename

for factor in {100,316,1000,3162,10000,31623,100000,316228,1000000}
do
    reps=$((10**7/factor))
    ./benchmark/alias -t 62 -C -q $((factor*10**3)) ${@:2} -r ${reps} | tee -a $filename
done

rm -f $linkname
ln -s $filename $linkname
