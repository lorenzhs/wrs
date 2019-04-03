/*******************************************************************************
 * wrs/parallel_do.cpp
 *
 * Pthreads parallelization helper
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <wrs/parallel_do.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

namespace wrs {

int g_num_numa_nodes;
int g_total_threads;
std::vector<tlx::ThreadPool*> global_pools;

void init_threads(int threads) {
    g_total_threads = threads;
#ifdef WRS_HAVE_NUMA
    g_num_numa_nodes = topoGetSystemNUMANodeCount();
    // if we're using less threads than we have numa nodes, cap number of nodes
    g_num_numa_nodes = std::min(g_total_threads, g_num_numa_nodes);
    int threads_per_node = (threads + g_num_numa_nodes - 1) / g_num_numa_nodes;

    int min = 0, max = threads_per_node;
    for (int i = 0; i < g_num_numa_nodes; i++) {
        int num_threads = max - min;
        sLOG1 << "Pool" << i << "has" << num_threads << "threads";
        global_pools.push_back(new tlx::ThreadPool(
                                   num_threads, [i]{ numa_run_on_node(i); }));
        min = max;
        max = std::min(threads, max + threads_per_node);
    }
#else
    g_num_numa_nodes = 1;
    global_pools.push_back(new tlx::ThreadPool(threads));
#endif
}

int get_num_threads() {
    return g_total_threads;
}


int get_num_nodes() {
    return g_num_numa_nodes;
}

tlx::ThreadPool* get_pool(size_t i) {
    return global_pools[i];
}

} // namespace wrs
