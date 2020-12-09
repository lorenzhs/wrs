/*******************************************************************************
 * tests/alias_common.hpp
 *
 * Alias table tests: utilities
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

#pragma once
#ifndef TESTS_ALIAS_COMMON_HEADER
#define TESTS_ALIAS_COMMON_HEADER

// Make sure NDEBUG isn't set
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <wrs/alias.hpp>
#include <wrs/dedup.hpp>
#include <wrs/gsl.hpp>
#include <wrs/multi_alias.hpp>
#include <wrs/outsens.hpp>
#include <wrs/par_scan_alias.hpp>
#include <wrs/par_scan_alias2.hpp>
#include <wrs/simple_scan_alias.hpp>
#include <wrs/verify.hpp>

#include <gtest/gtest.h>

#include <vector>


template <typename alias_t, typename count_t = int>
void do_sample(alias_t &table, std::vector<count_t> &out, size_t samples,
               size_t seed) {
    using result_type = typename alias_t::result_type;
    const typename alias_t::result_type empty = -1;

    size_t size = out.size();
    if constexpr (alias_t::yields_single_sample) {
        wrs::generators::select_t gen(seed);
        for (size_t i = 0; i < samples; i++) {
            result_type sample;
            if constexpr (alias_t::pass_rand == 1)
                sample = table.sample(gen.next());
            else if constexpr (alias_t::pass_rand == 2)
                sample = table.sample(gen.next(), gen.next());
            else
                sample = table.sample();
            EXPECT_NE(sample, empty) << alias_t::name << " invalid sample -1";
            EXPECT_LT(sample, size) << alias_t::name << " sample range fail";
            if (sample == empty || static_cast<size_t>(sample) >= size)
                continue;
            out[sample]++;
        }
    } else {
        auto callback = [&](auto sample, auto mult) {
            ASSERT_NE(sample, empty) << alias_t::name << " invalid sample -1";
            ASSERT_LT(sample, size) << alias_t::name << " sample range fail";
            if (sample == empty || static_cast<size_t>(sample) >= size)
                return;
            out[sample] += mult;
        };
        StochasticLib1 rng(seed);
        table.sample(rng, callback, samples);
    }
}

template <typename alias_t>
void construct(alias_t &table, const std::vector<double> &weights) {
    if constexpr (alias_t::init_with_seed) {
        table.init(weights.size(), 42);
    } else {
        table.init(weights.size());
    }
    table.construct(weights.begin(), weights.end());
}

template <template <typename> typename runner>
struct run_for_all {
    void operator()(bool skip_dedup = false) {
        runner<wrs::alias<>>()();
        if (!skip_dedup)
            runner<wrs::dedup<wrs::alias<>>>()();
        runner<wrs::multi_alias<wrs::alias>>()();
#ifdef WRS_HAVE_GSL
        runner<wrs::gsl_alias>()();
#endif
        runner<wrs::simple_scan_alias<>>()();
        runner<wrs::par_scan_alias<>>()();
        runner<wrs::par_scan_alias2<>>()();
        runner<wrs::outsens<int64_t, true, 256>>()();
        runner<wrs::outsens<uint64_t, true, 64>>()();
        runner<wrs::outsens<int32_t, true, 8>>()();
        runner<wrs::outsens<int32_t, false, 8>>()();
        runner<wrs::par_outsens<int32_t, true, 128>>()();
        runner<wrs::par_outsens<int32_t, false, 32>>()();
    }
};


// uniformly random in [0, 1)
struct uniform_gen {
    static constexpr const char *name = "uniform";
    uniform_gen(size_t seed = 42) : seed_(seed) {}
    double operator()(std::vector<double> &out, size_t size) const;
    size_t seed_;
};


// [1.1, 1-0.1/(n-1), ..., 1-0.1/(n-1)]
struct one_heavy_gen {
    static constexpr const char *name = "one_heavy";
    double operator()(std::vector<double> &out, size_t size) const;
};

// [0.9, 1+0.1/(n-1), ..., 1+0.1/(n-1)]
struct one_light_gen {
    static constexpr const char *name = "one_light";
    double operator()(std::vector<double> &out, size_t size) const;
};

// [n, 1, ..., 1]
struct one_heavy2_gen {
    static constexpr const char *name = "one_heavy2";
    double operator()(std::vector<double> &out, size_t size) const;
};


// [1^-exp, 2^-exp, 3^-exp, ..., n^-exp]
struct powerlaw_gen {
    static constexpr const char *name = "powerlaw";
    powerlaw_gen(double exp) : exp_(-exp) {}
    double operator()(std::vector<double> &out, size_t size) const;
    const double exp_;
};


// [1^-exp, 2^-exp, 3^-exp, ..., n^-exp], but shuffled
struct powerlaw_shuffle_gen {
    static constexpr const char *name = "powerlaw_shuffle";
    powerlaw_shuffle_gen(double exp) : exp_(exp) {}
    double operator()(std::vector<double> &out, size_t size) const;
    const double exp_;
};


#endif // TESTS_ALIAS_COMMON_HEADER
