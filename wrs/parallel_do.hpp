/*******************************************************************************
 * wrs/parallel_do.hpp
 *
 * Pthreads parallelization helper
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PARALLEL_DO_HEADER
#define WRS_PARALLEL_DO_HEADER

#include <wrs/accumulate.hpp>
#include <wrs/timer.hpp>
#include <wrs/util.hpp>

#include <tlx/logger.hpp>
#include <tlx/math/aggregate.hpp>
#include <tlx/thread_pool.hpp>

#ifdef WRS_HAVE_NUMA
#include <numa.h>
#include <silo.h>
#include <topo.h>
#endif

#include <algorithm>
#include <numeric>
#include <thread>
#include <vector>

namespace wrs {

void init_threads(int threads, bool numa = true);

void release_threads();

int get_num_threads();

int get_num_nodes();

int get_threads_per_node();

tlx::ThreadPool* get_pool(size_t i);

void wait_all_threads();

// Call function once on every NUMA node
template <typename F>
void do_per_numa_node(F&& callback) {
    for (int i = 0; i < get_num_nodes(); i++) {
        get_pool(i)->enqueue([callback, i]() { callback(i); });
    }
    wait_all_threads();
}

template <typename F>
void do_for_thread(F&& callback, int thread) {
    int node = thread / get_threads_per_node();
    get_pool(node)->enqueue([callback]() { callback(); });
}

template <typename F, typename size_type>
void parallel_do_block(F&& callback, size_type max_, size_type min_ = 0,
                       size_type block_size = 65536) {
    const bool debug = false;
    static_assert(std::is_integral_v<size_type>, "size_type must be integral");
    int num_threads = get_num_threads(), num_nodes = get_num_nodes(),
        threads_per_node = get_threads_per_node();
    size_type total_blocks = (max_ - min_ + block_size - 1) / block_size;

    const size_type min_blocks = num_threads;
    if (total_blocks < min_blocks) {
        sLOG << "Only got" << total_blocks << "blocks for" << num_threads
             << "threads! min=" << min_ << "max=" << max_;
        block_size = (max_ - min_ + min_blocks - 1) / min_blocks;
        total_blocks = (max_ - min_ + block_size - 1) / block_size;
        sLOG << "Changed block size to" << block_size << "->" << total_blocks
             << "blocks";
    }
    size_type blocks_per_node = (total_blocks + num_nodes - 1) / num_nodes;

    // If you build a machine with more than 16 NUMA nodes, you deserve this
    assert(num_nodes < 16);
    size_type last_block[16];
    std::atomic<size_type> next_block[16];

    for (int i = 0; i < num_nodes; i++) {
        last_block[i] = std::min(total_blocks, (i + 1) * blocks_per_node);
        next_block[i] = (i == 0) ? 0 : last_block[i - 1];
        sLOG << "Node" << i << "handling blocks" << next_block[i] << "to"
             << last_block[i];
    }

    auto worker = [&, callback](int numa_node, int thread) {
        tlx::Aggregate<double> block_stats;
        size_type my_next_block, my_blocks = 0;
        wrs::timer t_block, t_total;
        while ((my_next_block = next_block[numa_node].fetch_add(1)) <
               last_block[numa_node]) {
            sLOG << "Thread" << thread << "node" << numa_node
                 << "handling block" << my_next_block << "of" << total_blocks
                 << "last for this node:" << last_block[numa_node];

            ++my_blocks;
            size_type min = min_ + my_next_block * block_size,
                      max = std::min<size_type>(min + block_size, max_);
            t_block.reset();
            callback(min, max, thread);
            block_stats.add(t_block.get());
        }
        double total = t_total.get();
        // clang-format off
        LOG0 << "RESULT type=dynthread"
             << " total=" << total
             << " blocks=" << my_blocks
             << " blockavg=" << block_stats.avg()
             << " blockdev=" << block_stats.stdev()
             << " thread=" << thread
             << " node=" << numa_node
             << " threads=" << num_threads;
        // clang-format on
    };

    for (int i = 0; i < num_threads && static_cast<size_type>(i) < total_blocks;
         i++) {
        int numa_node = i / threads_per_node;
        get_pool(numa_node)->enqueue(
            [worker, i, numa_node]() { worker(numa_node, i); });
    }
    for (int i = 0; i < get_num_nodes(); i++) {
        if (get_pool(i)->size() > 0) {
            get_pool(i)->loop_until_empty();
        }
    }
}


template <typename F, typename size_type>
void parallel_do_range(F&& callback, size_type max_, size_type min_ = 0) {
    static_assert(std::is_integral_v<size_type>, "size_type must be integral");
    int num_threads = get_num_threads();
    int threads_per_node = get_threads_per_node();
    size_type elems_per_thread = (max_ - min_ + num_threads - 1) / num_threads;
    size_type min = min_, max = min_ + elems_per_thread;

    for (int i = 0; i < num_threads; i++) {
        get_pool(i / threads_per_node)->enqueue([callback, min, max, i]() {
            callback(min, max, i);
        });
        min = max;
        max = std::min(max_, max + elems_per_thread);
        // don't spawn empty threads
        if (min == max_)
            break;
    }
    for (int i = 0; i < get_num_nodes(); i++) {
        if (get_pool(i)->size() > 0) {
            get_pool(i)->loop_until_empty();
        }
    }
}

template <typename F, typename size_type>
void parallel_do(F&& callback, size_type max_, size_type min_ = 0) {
    static_assert(std::is_integral_v<size_type>, "size_type must be integral");
    parallel_do_range(
        [&callback](size_type min, size_type max, int /* thread */) {
            for (size_type i = min; i < max; i++)
                callback(i);
        },
        max_, min_);
}

template <typename InputIterator, typename OutputIterator>
void parallel_copy(InputIterator begin, InputIterator end, OutputIterator out_begin) {
    if (end - begin <= (1L << 22)) {
        std::copy(begin, end, out_begin);
        return;
    }
    parallel_do_range(
        [&begin, &out_begin](auto min, auto max, auto) {
            std::copy(begin + min, begin + max, out_begin);
        },
        end - begin);
}

template <typename Iterator,
          typename value_type = typename std::iterator_traits<Iterator>::value_type>
value_type parallel_accumulate(Iterator begin, Iterator end, value_type initial) {
    size_t size = end - begin;
    std::atomic<value_type> result = initial;
    auto worker = [&](size_t min, size_t max, int /* thread_id */) noexcept {
        value_type block_result =
            wrs::accumulate(begin + min, begin + max, value_type{});
        wrs::atomic_fetch_add(&result, block_result);
    };
    parallel_do_block(worker, size, /* min */ static_cast<size_t>(0),
                      /* blocksize */ static_cast<size_t>(1) << 18);
    return result;
}

} // namespace wrs

#endif // WRS_PARALLEL_DO_HEADER
