/*******************************************************************************
 * wrs/parallel_do.hpp
 *
 * Pthreads parallelization helper
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PARALLEL_DO_HEADER
#define WRS_PARALLEL_DO_HEADER

#include <tlx/thread_pool.hpp>

#ifdef WRS_HAVE_NUMA
#include <numa.h>
#include <silo.h>
#include <topo.h>
#endif

#include <numeric>
#include <thread>
#include <vector>

namespace wrs {

void init_threads(int threads);

int get_num_threads();

int get_num_nodes();

tlx::ThreadPool* get_pool(size_t i);

// Call function once on every NUMA node
template <typename F>
void do_per_numa_node(F&& callback) {
    for (int i = 0; i < get_num_nodes(); i++) {
        get_pool(i)->enqueue([callback, i]() { callback(i); });
    }
    for (int i = 0; i < get_num_nodes(); i++) {
        if (get_pool(i)->size() > 0) {
            get_pool(i)->loop_until_empty();
        }
    }
}

template <typename F>
void parallel_do_range(F&& callback, size_t max_, size_t min_ = 0) {
    int num_threads = get_num_threads();
    int threads_per_node = (num_threads + get_num_nodes() - 1) / get_num_nodes();
    size_t elems_per_thread = (max_ - min_ + num_threads - 1) / num_threads;
    size_t min = min_, max = min_ + elems_per_thread;

    for (int i = 0; i < num_threads; i++) {
        get_pool(i / threads_per_node)->enqueue(
            [callback, min, max, i]() { callback(min, max, i); });
        min = max;
        max = std::min(max_, max + elems_per_thread);
        // don't spawn empty threads
        if (min == max_) break;
    }
    for (int i = 0; i < get_num_nodes(); i++) {
        if (get_pool(i)->size() > 0) {
            get_pool(i)->loop_until_empty();
        }
    }
}

template <typename F, typename G>
void parallel_do(F&& init, G&& callback, size_t max_, size_t min_ = 0) {
    auto worker = [&init, &callback](size_t min, size_t max, int thread_id)
        noexcept(noexcept(callback) && noexcept(init))
    {
        auto state = init(min, max, thread_id);

        for (size_t i = min; i < max; ++i) {
            callback(state, i);
        }
    };
    parallel_do_range(worker, max_, min_);
}

template <typename F>
void parallel_do(F&& callback, size_t num_elems) {
    parallel_do(
        [](auto, auto, auto){ return nullptr; },
        [&callback](auto, size_t i) { callback(i); },
        num_elems);
}

template <typename InputIterator, typename OutputIterator>
void parallel_copy(InputIterator begin, InputIterator end,
                   OutputIterator out_begin) {
    if (end - begin <= (1L << 22)) {
        std::copy(begin, end, out_begin);
        return;
    }
    parallel_do_range(
        [&begin, &out_begin](size_t min, size_t max, auto) {
            std::copy(begin + min, begin + max, out_begin);
        }, end - begin);
}

template <typename Iterator, typename value_type =
          typename std::iterator_traits<Iterator>::value_type>
value_type parallel_accumulate(Iterator begin, Iterator end,
                               value_type initial) {
    size_t size = end - begin;
    std::vector<value_type> results(get_num_threads());
    auto worker =
        [&](size_t min, size_t max, int thread_id) noexcept {
            results[thread_id] =
                std::accumulate(begin + min, begin + max, value_type{});
        };
    parallel_do_range(worker, size);
    return std::accumulate(results.begin(), results.end(), initial);
}


template <typename InputIterator, typename OutputIterator, typename Op>
void parallel_scan(InputIterator begin, InputIterator end,
                   OutputIterator out_begin, Op &&op)
{
    size_t size = end - begin;
    using value_type = decltype(op(*begin, *begin));
    std::vector<value_type> temp(get_num_threads() + 1);
    auto phase1 = [&](size_t min, size_t max, int thread_id) {
        std::partial_sum(begin + min, begin + max, out_begin + min, op);
        temp[thread_id + 1] = *(out_begin + max - 1);
    };
    parallel_do_range(phase1, size);

    std::partial_sum(temp.begin(), temp.end(), temp.begin(), op);

    auto phase2 = [&](auto partial, size_t index) {
        *(out_begin + index) = op(partial, *(out_begin + index));
    };
    parallel_do(
        [&](size_t, size_t, int thread_id) { return temp[thread_id]; },
        phase2, size);
}

} // namespace wrs

#endif // WRS_PARALLEL_DO_HEADER
