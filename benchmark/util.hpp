/*******************************************************************************
 * benchmark/util.hpp
 *
 * Benchmark utilities
 *
 * Copyright (C) 2018-2020 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_UTIL_HEADER
#define BENCHMARK_UTIL_HEADER

#include <tlx/logger.hpp>
#include <tlx/math/aggregate.hpp>

#include <cassert>
#include <vector>

template <typename Benchmark>
std::vector<tlx::Aggregate<double>>
run_benchmark(Benchmark &&benchmark, int repetitions, int warmup_reps = 1,
              bool is_warmup = false, int iteration = -1) {
    std::vector<tlx::Aggregate<double>> timers;

    for (int i = (-1 + is_warmup) * warmup_reps; i < repetitions; i++) {
        const bool warmup = i < 0 || is_warmup;
        std::vector<double> results = benchmark(!warmup, iteration, i);

        if (timers.empty()) {
            timers.resize(results.size());
        }
        assert(timers.size() == results.size());
        if (!warmup) {
            // skip warmup repetitions
            for (size_t j = 0; j < results.size(); ++j) {
                timers[j].add(results[j]);
            }
        }
    }
    return timers;
}


template <typename Benchmark, typename Callback, typename Init>
std::vector<tlx::Aggregate<double>>
run_benchmark(Benchmark &&benchmark, Callback &&callback, Init &&init,
              int iterations, int repetitions, int warmup_its = 1,
              int warmup_reps = 1) {
    std::vector<tlx::Aggregate<double>> stats;

    for (int i = -warmup_its; i < iterations; i++) {
        const bool warmup = i < 0;
        init();
        int reps = repetitions;
        // Limit number of warmup reps
        if (warmup) {
            if (reps > 30)
                reps /= 2;
            if (reps >= 100)
                reps /= 5;
        }
        auto it_stats = run_benchmark(benchmark, reps, warmup_reps, warmup, i);
        if (warmup)
            continue; // skip stats collection and callback

        if (stats.empty()) {
            stats = it_stats;
        } else {
            for (size_t j = 0; j < stats.size(); ++j) {
                stats[j] += it_stats[j];
            }
        }
        callback(i, it_stats);
    }
    return stats;
}


namespace tlx {

template <typename T>
std::ostream &operator<<(std::ostream &os, const Aggregate<T> &x) {
    if (x.count() > 1) {
        os << "avg=" << x.avg() << " stdev=" << x.stdev() << " range=["
           << x.min() << ".." << x.max() << "]";
    } else {
        os << x.avg();
    }
    return os;
}

} // namespace tlx

#endif // BENCHMARK_UTIL_HEADER
