/*******************************************************************************
 * wrs/numa_array.hpp
 *
 * Transparently distributed array using Silo
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_NUMA_ARRY_HEADER
#define WRS_NUMA_ARRY_HEADER

#include <wrs/memory.hpp>
#include <wrs/parallel_do.hpp>

#include <tlx/logger.hpp>

#ifdef WRS_HAVE_NUMA
#include <silo.h>
#include <topo.h>
#endif

#include <memory>

namespace wrs {

#ifdef WRS_HAVE_NUMA
namespace __detail {
//! Private. Do not call directly.
void* __numa_alloc(size_t bytes, bool align) {
    constexpr bool debug = true;

    size_t bytes_per_thread = (bytes + g_total_threads - 1) / g_total_threads;
    int threads_per_node = (g_total_threads + g_num_numa_nodes - 1) /
        g_num_numa_nodes;

    sLOG << "Allocating" << bytes << "bytes on" << g_num_numa_nodes
         << "NUMA nodes";

    SSiloMemorySpec* specs = (SSiloMemorySpec*)malloc(
        sizeof(SSiloMemorySpec) * g_num_numa_nodes);
    if (specs == nullptr) { // abort
        assert(false);
        return nullptr;
    }
    // Initialize the specifications.
    int min = 0, max = threads_per_node;
    for (int i = 0; i < g_num_numa_nodes; ++i) {
        // Align sizes to 2MB
        size_t size = bytes_per_thread * (max - min);
        if (align) size = align_size(size, 2048 * 1024);
        specs[i].size = size;
        specs[i].numaNode = i;
        sLOG << "Thereof" << specs[i].size << "bytes on node" << i;

        min = max;
        max = std::min(g_total_threads, max + threads_per_node);
    }

    // Allocate the multi-node array. Uses transparent hugepages.
    void *buffer = siloMultinodeArrayAlloc(g_num_numa_nodes, specs);

    return buffer;
}
} // namespace __detail
#endif


// Allocate an array distributed over the available NUMA nodes
void* numa_alloc(size_t bytes, bool align = true, bool local = false) {
#ifndef WRS_HAVE_NUMA
    (void) align;
    (void) local;
    return allocate(bytes);
#else
    constexpr bool debug = true;

    void* buffer;
    // Handle non-NUMA systems like local allocations
    if (local || g_num_numa_nodes == 1) {
        // Allocate on a single node
        size_t size = bytes;
        if (align) size = align_size(size, 2048*1024);
        buffer = siloSimpleBufferAllocLocal(size);
    } else {
        // Do the actual distributed NUMA allocation
        buffer = __detail::__numa_alloc(bytes, align);
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

// A struct that fulfills the deleter requirements of std::unique_ptr
struct numa_deleter {
    template <typename T>
    void operator()(T* ptr) {
        numa_free((void*)ptr);
    }
};

// A type definition for easy use
template <typename T>
using numa_arr_ptr = std::unique_ptr<T[], numa_deleter>;

// Helper function to create a numa_arr_ptr akin to std::make_unique
template <typename T, typename ...Args>
numa_arr_ptr<T> make_numa_arr(size_t num_elems, Args...args) {
    T* ptr = static_cast<T*>(numa_alloc(num_elems * sizeof(T), std::forward<Args>(args)...));
    return numa_arr_ptr<T>(ptr);
}

} // namespace wrs

#endif // WRS_NUMA_ARRY_HEADER
