/*******************************************************************************
 * wrs/memory.hpp
 *
 * Memory helpers
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_MEMORY_HEADER
#define WRS_MEMORY_HEADER

#include <tlx/logger.hpp>

#include <sys/mman.h> // madvise
#include <cstdlib>
#include <memory>

namespace wrs {

constexpr size_t align_size(size_t size, size_t alignment) {
    return ((size + alignment - 1) / alignment) * alignment;
}

void* alloc_hugepage(size_t size) {
    constexpr bool debug = false;

    constexpr size_t alignment = 2 * 1024 * 1024;
    size_t bytes = align_size(size, alignment);
    sLOG << "Allocating" << bytes << "bytes instead of" << size;
    void* ptr = aligned_alloc(alignment, bytes);
    madvise(ptr, bytes, MADV_HUGEPAGE);
    return ptr;
}

// Allocate memory, using huge pages for allocations larger than 1MB
void* allocate(size_t size) {
    if (size >= 1024 * 1024) {
        return alloc_hugepage(size);
    } else {
        return malloc(size);
    }
}

struct deallocator {
    template <typename T>
    void operator()(T *ptr) {
        free((void*)ptr);
    }
};

// to avoid alloc-dealloc-mismatches
template <typename T>
using alloc_arr_ptr = std::unique_ptr<T[], deallocator>;

} // namespace wrs

#endif // WRS_MEMORY_HEADER
