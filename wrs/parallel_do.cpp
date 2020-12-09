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

#include <algorithm>
#include <vector>

namespace wrs {

int g_num_numa_nodes;
int g_total_threads;
std::vector<tlx::ThreadPool*> global_pools;

void init_threads(int threads, bool numa) {
    g_total_threads = threads;
#ifdef WRS_HAVE_NUMA
    if (numa) {
        // pin main thread to numa node 0
        numa_run_on_node(0);

        g_num_numa_nodes = topoGetSystemNUMANodeCount();
        // if we're using less threads than we have numa nodes, cap number of nodes
        g_num_numa_nodes = std::min(g_total_threads, g_num_numa_nodes);
        int threads_per_node = (threads + g_num_numa_nodes - 1) / g_num_numa_nodes;

        int min = 0, max = threads_per_node;
        for (int i = 0; i < g_num_numa_nodes; i++) {
            int num_threads = max - min;
            sLOG1 << "Pool" << i << "has" << num_threads << "threads";
            global_pools.push_back(
                new tlx::ThreadPool(num_threads, [i] { numa_run_on_node(i); }));
            min = max;
            max = std::min(threads, max + threads_per_node);
        }
    } else {
#endif
        (void)numa; // suppress unused variable warnung
        g_num_numa_nodes = 1;
        global_pools.push_back(new tlx::ThreadPool(threads));
#ifdef WRS_HAVE_NUMA
    }
#endif
}

void release_threads() {
    // works for both numa and non-numa
    for (int i = 0; i < g_num_numa_nodes; i++) {
        global_pools[i]->loop_until_empty();
        global_pools[i]->terminate();
        delete global_pools[i];
        global_pools[i] = nullptr;
    }
    global_pools.clear();
    g_num_numa_nodes = 0;
    g_total_threads = 0;
}

int get_num_threads() {
    return g_total_threads;
}


int get_num_nodes() {
    return g_num_numa_nodes;
}

int get_threads_per_node() {
    return (get_num_threads() + get_num_nodes() - 1) / get_num_nodes();
}

tlx::ThreadPool* get_pool(size_t i) {
    return global_pools[i];
}

void wait_all_threads() {
    for (int i = 0; i < get_num_nodes(); i++) {
        if (get_pool(i)->size() > 0) {
            get_pool(i)->loop_until_empty();
        }
    }
}

} // namespace wrs
