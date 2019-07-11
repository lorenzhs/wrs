# Weighted Random Sampling

This is the companion repository to "Parallel Weighted Random Sampling" by Lorenz Hübschle-Schneider and Peter Sanders.

To reproduce our experiments, compile using cmake and execute [benchmark/bench_all.sh](benchmark/bench_all.sh) from your build directory.
Note that this script targets a machine with 80 cores (160 threads).  A version for 32 cores (64 threads) is provided in [benchmark/bench_all_64.sh](benchmark/bench_all_64.sh).
Adjust this to your machine by changing the values to the `-t` parameter in `benchmark/run_*.sh`.
You should also change the output paths in `benchmark_run_*.sh`

Licensed under GPLv3, see the LICENSE file.
