/*******************************************************************************
 * wrs/numa_array.hpp
 *
 * Transparently distributed array using Silo
 *
 * Copyright (C) 2018 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <wrs/numa_array.hpp>
#include <wrs/memory.hpp>
#include <wrs/parallel_do.hpp>

#include <tlx/logger.hpp>

#ifdef WRS_HAVE_NUMA
#include <silo.h>
#include <topo.h>
#endif

namespace wrs {

#ifdef WRS_HAVE_NUMA
//! Private. Do not call directly.
static void* __numa_alloc(size_t bytes, bool align) {
    constexpr bool debug = true;

    size_t bytes_per_thread = (bytes + get_num_threads() - 1) / get_num_threads();
    int threads_per_node = (get_num_threads() + get_num_nodes() - 1) /
        get_num_nodes();

    sLOG << "Allocating" << bytes << "bytes on" << get_num_nodes()
         << "NUMA nodes";

    SSiloMemorySpec* specs = (SSiloMemorySpec*)malloc(
        sizeof(SSiloMemorySpec) * get_num_nodes());
    if (specs == nullptr) { // abort
        assert(false);
        return nullptr;
    }
    // Initialize the specifications.
    int min = 0, max = threads_per_node;
    for (int i = 0; i < get_num_nodes(); ++i) {
        // Align sizes to 2MB
        size_t size = bytes_per_thread * (max - min);
        if (align) size = align_size(size, 2048 * 1024);
        specs[i].size = size;
        specs[i].numaNode = i;
        sLOG << "Thereof" << specs[i].size << "bytes on node" << i;

        min = max;
        max = std::min(get_num_threads(), max + threads_per_node);
    }

    // Allocate the multi-node array. Uses transparent hugepages.
    void *buffer = siloMultinodeArrayAlloc(get_num_nodes(), specs);

    return buffer;
}
#endif

// Allocate an array distributed over the available NUMA nodes
void* numa_alloc(size_t bytes, bool align, bool local) {
#ifndef WRS_HAVE_NUMA
    (void) align;
    (void) local;
    return allocate(bytes);
#else
    constexpr bool debug = true;

    void* buffer;
    // Handle non-NUMA systems like local allocations
    if (local || get_num_nodes() == 1) {
        // Allocate on a single node
        size_t size = bytes;
        if (align) size = align_size(size, 2048*1024);
        buffer = siloSimpleBufferAllocLocal(size);
    } else {
        // Do the actual distributed NUMA allocation
        buffer = __numa_alloc(bytes, align);
    }

    // If for some reason the NUMA allocation failed, fall back to the simple
    // version
    if (buffer == nullptr) {
        LOG << "NUMA aware allocation failed; falling back to simple";
        buffer = allocate(bytes);
        assert(buffer != nullptr);
    }
    return buffer;
#endif
}

// Pointers allocated with silo need to be freed with `siloFree`.  This also
// frees the SSiloMemorySpec object (Silo tracks these internally)
void numa_free(void* ptr) {
    if (ptr != nullptr) {
#if WRS_HAVE_NUMA
        siloFree(ptr);
#else
        free(ptr);
#endif
    }
}

} // namespace wrs
