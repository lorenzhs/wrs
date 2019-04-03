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

#include <memory>
#include <utility>

namespace wrs {

// Allocate an array distributed over the available NUMA nodes
void* numa_alloc(size_t bytes, bool align = true, bool local = false);

// Pointers allocated with silo need to be freed with `siloFree`.  This also
// frees the SSiloMemorySpec object (Silo tracks these internally)
void numa_free(void* ptr);

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
