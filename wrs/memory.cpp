/*******************************************************************************
 * wrs/memory.cpp
 *
 * Memory helpers
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <wrs/memory.hpp>

namespace wrs {

// Allocate memory using huge pages
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

} // namespace wrs
