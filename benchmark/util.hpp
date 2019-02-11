/*******************************************************************************
 * benchmark/util.hpp
 *
 * Benchmark utilities
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_UTIL_HEADER
#define BENCHMARK_UTIL_HEADER

#include <wrs/aggregate.hpp>

#include <vector>

template <typename Benchmark>
std::vector<wrs::Aggregate<double>>
run_benchmark(Benchmark &&benchmark, int repetitions, int warmup_reps = 1,
              bool is_warmup = false) {
    std::vector<wrs::Aggregate<double>> timers;

    for (int i = (-1 + is_warmup) * warmup_reps; i < repetitions; i++) {
        const bool warmup = i < 0 || is_warmup;
        std::vector<double> results = benchmark(!warmup);

        if (timers.empty()) {
            timers.resize(results.size());
        }
        assert(timers.size() == results.size());
        if (!warmup) {
            // skip warmup repetitions
            for (size_t j = 0; j < results.size(); ++j) {
                timers[j].Add(results[j]);
            }
        }
    }
    return timers;
}


template <typename Benchmark, typename Callback, typename Init>
std::vector<wrs::Aggregate<double>>
run_benchmark(Benchmark &&benchmark, Callback &&callback, Init &&init,
              int iterations, int repetitions,
              int warmup_its = 1, int warmup_reps = 1)
{
    std::vector<wrs::Aggregate<double>> stats;

    for (int i = -warmup_its; i < iterations; i++) {
        const bool warmup = i < 0;
        init();
        auto it_stats = run_benchmark(benchmark, repetitions,
                                      warmup_reps, warmup);
        if (warmup) continue; // skip stats collection and callback

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


#endif // BENCHMARK_UTIL_HEADER
