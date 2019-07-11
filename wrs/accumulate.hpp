/*******************************************************************************
 * wrs/accumulate.hpp
 *
 * Accumulation
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_ACCUMULATE_HEADER
#define WRS_ACCUMULATE_HEADER

#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>

namespace wrs {

namespace _detail {

template <typename value_type>
value_type accumulate_clobber(value_type *first, value_type *last, value_type init)
{
    constexpr static size_t blocksize = 128;

    const size_t size = last - first;
    if (size < blocksize)
        return std::accumulate(first, last, init);

    const size_t num_blocks = (size + blocksize - 1) / blocksize;

    for (size_t block = 0; block < num_blocks; block++) {
        *(first + block) = std::accumulate(
            first + block * blocksize,
            std::min(first + (block + 1) * blocksize, last),
            value_type{});
    }
    return accumulate_clobber(first, first + num_blocks, init);
}

} // namespace _detail

template <typename Iterator,
          typename value_type = typename std::iterator_traits<Iterator>::value_type>
value_type accumulate(Iterator first, Iterator last, value_type init) {
    constexpr static size_t blocksize = 128;

    const size_t size = last - first;
    if (size < blocksize)
        return std::accumulate(first, last, init);

    const size_t num_blocks = (size + blocksize - 1) / blocksize;
    auto blocksum = std::make_unique<value_type[]>(num_blocks);

    for (size_t block = 0; block < num_blocks; block++) {
        blocksum[block] = std::accumulate(
            first + block * blocksize,
            std::min(first + (block + 1) * blocksize, last),
            value_type{});
    }
    // recurse without further allocations
    return _detail::accumulate_clobber(blocksum.get(), blocksum.get() + num_blocks, init);
}


} // namespace wrs

#endif // WRS_ACCUMULAvalue_typeE_HEADER
