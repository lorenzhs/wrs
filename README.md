# Parallel Weighted Random Sampling

This code implements sequential and parallel methods for weighted random sampling, as described in our eponymous paper: *Hübschle-Schneider, L., & Sanders, P. (2019). Parallel Weighted Random Sampling. In 27th Annual European Symposium on Algorithms (ESA 2019).*

The current state of the repository represents the most up-to-date version and was used for the evaluation of the weighted sampling chapter in my dissertation (link to follow).  The version used for the experiments in our ESA paper is available as the [`esa` tag](https://github.com/lorenzhs/wrs/tree/esa).

The paper is freely accessible under a Creative Commons license (CC-BY) [on the publisher's website](https://drops.dagstuhl.de/opus/volltexte/2019/11180/).

If you use this library in the context of an academic publication, we ask that you cite our paper:
```bibtex
@InProceedings{HubSan2019weightedsampling,
  author =	{Lorenz H{\"u}bschle-Schneider and Peter Sanders},
  title =	{{Parallel Weighted Random Sampling}},
  booktitle =	{27th Annual European Symposium on Algorithms (ESA 2019)},
  pages =	{59:1--59:24},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-124-5},
  ISSN =	{1868-8969},
  year =	{2019},
  volume =	{144},
  editor =	{Michael A. Bender and Ola Svensson and Grzegorz Herman},
  publisher =	{Schloss Dagstuhl--Leibniz-Zentrum fuer Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{http://drops.dagstuhl.de/opus/volltexte/2019/11180},
  doi =		{10.4230/LIPIcs.ESA.2019.59},
}
```

### Building

Build with cmake (version 3.9.2 or later is required). Remember to fetch the submodules before compiling: `git submodule update --init`. A compiler compatible with C++17 is required.

**Optional Dependencies:** [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/mkl) for faster random variate generation, [libnuma](https://github.com/numactl/numactl) for Non-Uniform Memory Access (NUMA) awareness, [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/) for their alias table implementation (only for comparison in the experiments).

### Experiments

To reproduce our experiments, compile using cmake and execute the benchmark scripts. The scripts are customised to the machines we used for our experiments.  On our Intel machine with 80 cores and 160 threads, we ran [benchmark/bench_intel.sh](benchmark/bench_intel.sh) from the build directory.  On our AMD machine with 32 cores and 64 threads, we ran [benchmark/bench_all_64.sh](benchmark/bench_all_64.sh).  Adjust these to your machine by changing the thread counts over which the scripts iterate.  You might also want to change the output paths in `benchmark/run_*.sh`

**[License](/LICENSE):** GPLv3
