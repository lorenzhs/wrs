/*******************************************************************************
 * wrs/prefix_sum.hpp
 *
 * prefix sums
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PREFIX_SUM_HEADER
#define WRS_PREFIX_SUM_HEADER

#include <wrs/parallel_do.hpp>

#include <type_traits>
#include <vector>

namespace wrs {

// Exclusive scan, out_begin may be equal to begin (in-place)
template <typename InputIterator, typename OutputIterator, typename Op,
          typename T = typename std::iterator_traits<InputIterator>::value_type>
OutputIterator exclusive_scan(InputIterator begin, InputIterator end,
                              OutputIterator out_begin, Op &&op, T initial,
                              bool include_last = false /* write last item? */) {
    const bool debug = false;

    sLOG << "exclusive_scan over" << end - begin << "elements"
         << (include_last ? "including" : "excluding") << "last write";
    if (begin == end)
        return out_begin;

    T temp = *begin;
    *out_begin = initial;
    InputIterator it = begin + 1;
    OutputIterator out_it = out_begin + 1;
    while (it != end) {
        T curr = *it++;
        *out_it++ = temp;
        LOG << "out[" << out_it - out_begin - 1 << "] = " << temp << " was " << curr;
        temp = op(temp, curr);
    }
    if (include_last) {
        *out_it++ = temp;
    }
    return out_it;
}

template <typename InputIterator, typename OutputIterator, typename Op>
void parallel_scan(InputIterator begin, InputIterator end,
                   OutputIterator out_begin, Op &&op) {
    size_t size = end - begin;
    using value_type = decltype(op(*begin, *begin));
    std::vector<value_type> temp(get_num_threads() + 1);
    auto phase1 = [&](size_t min, size_t max, int thread_id) {
        std::partial_sum(begin + min, begin + max, out_begin + min, op);
        temp[thread_id + 1] = *(out_begin + max - 1);
    };
    parallel_do_range(phase1, size);

    std::partial_sum(temp.begin(), temp.end(), temp.begin(), op);

    auto phase2 = [&](size_t min, size_t max, int thread_id) {
        OutputIterator it = out_begin + min, end = out_begin + max;
        value_type partial = temp[thread_id];
        while (it != end) {
            *it = op(partial, *it);
            it++;
        }
    };
    parallel_do_range(phase2, size);
}


} // namespace wrs

#endif // WRS_PREFIX_SUM_HEADER
