#!/bin/bash
set -e

filename=/global_data/lorenz/wrs/gslresults_$(date "+%F.%H-%M-%S").txt

echo "Invocation: $0 $*" | tee $filename
echo "Running on $HOST on $(date)" | tee -a $filename

./benchmark/alias -t 1 -M -P -D -n $((10**7)) -q $((10**7)) -i 10 -r 20 -s 11235813213455 | tee -a $filename

linkname=/global_data/lorenz/wrs/gslresults.txt
rm -f $linkname
ln -s $filename $linkname
