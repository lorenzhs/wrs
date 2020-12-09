/*******************************************************************************
 * benchmark/alias.cpp
 *
 * Parallel Weighted Random Sampling benchmarks
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <benchmark/util.hpp>

#include <wrs/alias.hpp>
#include <wrs/dedup.hpp>
#include <wrs/generators/select.hpp>
#include <wrs/multi_alias.hpp>
#include <wrs/outsens.hpp>
#include <wrs/par_scan_alias.hpp>
#include <wrs/par_scan_alias2.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/simple_scan_alias.hpp>
#include <wrs/timer.hpp>
#include <wrs/verify.hpp>

// GSL is optional, only include it if it was found
#ifdef WRS_HAVE_GSL
#include <wrs/gsl.hpp>
#endif

#include <tlx/cmdline_parser.hpp>
#include <tlx/logger.hpp>
#include <tlx/math/aggregate.hpp>
#include <tlx/thread_pool.hpp>

#include <algorithm>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>

#ifdef WRS_HAVE_NUMA
static constexpr bool numa = true;
#else
static constexpr bool numa = false;
#endif

static constexpr bool debug = true;

template <typename alias_table, typename input_generator>
tlx::Aggregate<double> run_alias(input_generator &&input_gen, size_t input_size,
                                 int iterations, int repetitions, size_t seed,
                                 const std::string &name) {
    LOG << "";
    LOG << "Constructing alias table using " << name << " method";

    wrs::numa_arr_ptr<double> weights;
    alias_table table;

    auto init = [&]() {
        weights = input_gen(seed, input_size);
        if constexpr (alias_table::init_with_seed) {
            table.init(input_size, seed);
        } else {
            table.init(input_size);
        }
    };
    auto runner = [&](bool no_warmup, int, int) {
        const bool debug = no_warmup;
        wrs::timer timer;
        table.construct(weights.get(), weights.get() + input_size);
        double duration = timer.get();
        wrs::verify(weights.get(), weights.get() + input_size, table);
        auto timers = table.get_timers();

        // clang-format off
        LOG << "RESULT type=construction"
            << " total=" << duration
            << " timers=" << timers
            << " size=" << input_size
            << " method=" << name
            << " threads=" << wrs::get_num_threads() // only relevant for parallel...
            << " numa=" << numa
            << " seed=" << seed;
        // clang-format on

        timers.insert(timers.begin(), duration);
        return timers;
    };
    auto logger = [&](int it, auto &stats) {
        LOG << "It " << it << " construction: " << std::fixed
            << std::setprecision(3) << stats[0] << " seed=" << seed;
        ++seed;
    };
    auto stats = run_benchmark(runner, logger, init, iterations, repetitions,
                               /* warmup its */ 1, /* warmup reps */ 1);
    LOG << "Construction: " << std::fixed << std::setprecision(3) << stats[0] << " ("
        << iterations << " iterations x " << repetitions << " repetitions)";
    if (stats.size() == 9) { // psa2
        LOG << "\tPreprocessing:  " << stats[1];
        LOG << "\tGreedy fill:    " << stats[2];
        LOG << "\tClassification: " << stats[3];
        LOG << "\tClass. prefsum: " << stats[4];
        LOG << "\tIndex comp:     " << stats[5];
        LOG << "\tWeight prefsum: " << stats[6];
        LOG << "\tSplitting:      " << stats[7];
        LOG << "\tBase cases:     " << stats[8];
    } else if (stats.size() == 8) { // parallel scan
        LOG << "\tPreprocessing:  " << stats[1];
        LOG << "\tClassification: " << stats[2];
        LOG << "\tClass. prefsum: " << stats[3];
        LOG << "\tIndex comp:     " << stats[4];
        LOG << "\tWeight prefsum: " << stats[5];
        LOG << "\tSplitting:      " << stats[6];
        LOG << "\tBase cases:     " << stats[7];
    } else if (stats.size() == 3) { // multi
        LOG << "\tSubtable constr: " << stats[1];
        LOG << "\tTop-level table: " << stats[2];
    } else if (stats.size() == 4) { // sequential
        LOG << "\tPreprocessing:  " << stats[1];
        LOG << "\tClassification: " << stats[2];
        LOG << "\tAssignment:     " << stats[3];
    } else if (stats.size() == 2) { // GSL
        LOG << "\tConstruction:   " << stats[1];
    }
    return stats[0];
}

template <typename alias_table, typename input_generator>
tlx::Aggregate<double> run_queries(input_generator &&input_gen, size_t input_size,
                                   int iterations, int repetitions, size_t queries,
                                   size_t seed, const std::string &name) {
    LOG << "";
    LOG << "Querying alias table using " << name << " method";

    // construct table
    wrs::numa_arr_ptr<double> weights;
    alias_table table;

    auto init = [&]() {
        wrs::timer timer, total_timer;
        weights = input_gen(seed, input_size);
        double t_gen = timer.get_and_reset();

        if constexpr (alias_table::init_with_seed) {
            table.init(input_size, seed);
        } else {
            table.init(input_size);
        }
        double t_init = timer.get_and_reset();

        table.construct(weights.get(), weights.get() + input_size);
        double t_cons = timer.get_and_reset();

        wrs::verify(weights.get(), weights.get() + input_size, table);
        double t_verify = timer.get_and_reset();

        // clang-format off
        LOG1 << "RESULT type=qcons"
             << " total=" << total_timer.get()
             << " gen=" << t_gen
             << " init=" << t_init
             << " cons=" << t_cons
             << " verify=" << t_verify
             << " size=" << input_size
             << " method=" << name
             << " threads=" << wrs::get_num_threads()
             << " numa=" << numa
             << " seed=" << seed;
        // clang-format on
    };

    auto runner = [&](const bool debug, int, int) {
        wrs::timer timer;

        // Run queries in parallel
        size_t result_size = queries;
        if constexpr (alias_table::yields_single_sample) {
            wrs::parallel_do_block(
                [&](size_t min, size_t max, auto) {
                    std::vector<typename alias_table::result_type> dummy(1);
                    wrs::generators::select_t rng(seed + min + 1);
                    for (size_t i = min; i < max; i++) {
                        if constexpr (alias_table::pass_rand == 1) {
                            dummy[0] = table.sample(rng.next());
                        } else if constexpr (alias_table::pass_rand == 2) {
                            dummy[0] = table.sample(rng.next(), rng.next());
                        } else {
                            dummy[0] = table.sample();
                        }
                    }
                },
                queries);
        } else {
            thread_local std::pair<typename alias_table::result_type, size_t> dummy;
            auto callback = [](auto sample, auto mult) {
                sLOG0 << "Element" << sample << "multiplicity" << mult;
                dummy = std::make_pair(sample, mult);
            };
            StochasticLib1 rng(seed);
            result_size = table.sample(rng, callback, queries);
        }
        double duration = timer.get();

        // clang-format off
        LOG << "RESULT type=query"
            << " total=" << duration
            << " size=" << input_size
            << " queries=" << queries
            << " osize=" << result_size
            << " method=" << name
            << " threads=" << wrs::get_num_threads()
            << " numa=" << numa
            << " seed=" << seed;
        // clang-format on

        return std::vector<double>{duration};
    };
    auto logger = [&](int it, auto &stats) {
        LOG << "It " << it << " queries: " << std::fixed << std::setprecision(3)
            << stats[0] << " seed=" << seed;
        ++seed;
    };
    auto stats = run_benchmark(runner, logger, init, iterations, repetitions,
                               /* warmup its */ 0, /* warmup reps */ 1);
    LOG << "Queries: " << std::fixed << std::setprecision(3) << stats[0];
    return stats[0];
}


int main(int argc, char *argv[]) {
    tlx::CmdlineParser clp;

    size_t input_size = 100, queries = 1000, seed = 0, num_threads = 0;
    int iterations = 2, repetitions = 5;
    bool normalize = false, no_seq = false, no_dedup = false, no_multi = false,
         no_multisweep = false, no_simplescan = false, no_psa = false,
         no_psa2 = false, no_outsens = false, no_outsensnd = false,
         no_gsl = false, no_queries = false, no_construct = false;
    double powerlaw_exp = 0;
    clp.add_size_t('n', "input-size", input_size, "input size");
    clp.add_int('i', "iterations", iterations, "iterations");
    clp.add_int('r', "repetitions", repetitions, "repetitions per iteration");
    clp.add_size_t('q', "queries", queries, "number of queries");
    clp.add_size_t('s', "seed", seed, "random generator seed (0 for random)");
    clp.add_bool('z', "normalize", normalize, "normalize inputs (for debuggability");
    clp.add_double('p', "powerlaw", powerlaw_exp,
                   "powerlaw distribution exponent (instead of uniform)");
    clp.add_bool('Q', "noqueries", no_queries, "don't run query benchmark");
    clp.add_bool('C', "noconstruct", no_construct,
                 "don't run construction benchmark");
    clp.add_bool('S', "noseq", no_seq, "don't run sequential version");
    clp.add_bool('D', "nodedup", no_dedup,
                 "don't run deduplicating sequential version");
    clp.add_bool('M', "nomulti", no_multi, "don't run multi+alias version");
    clp.add_bool('L', "nomultisweep", no_multisweep,
                 "don't run multi+sweep version");
    clp.add_bool('E', "nosimplescan", no_simplescan,
                 "don't run simple scanning version (sweep)");
    clp.add_bool('B', "nopsa", no_psa, "don't run parallel scanning version");
    clp.add_bool('P', "nopsa2", no_psa2, "don't run fast parallel scanning version");
    clp.add_bool('O', "nooutsens", no_outsens, "don't run output sensitive version");
    clp.add_bool('R', "nooutsensnd", no_outsensnd,
                 "don't run non-deduplicating output-sensitive version");
    clp.add_bool('G', "nogsl", no_gsl, "don't run GSL version");
    clp.add_size_t('t', "threads", num_threads, "number of threads");

    if (!clp.process(argc, argv))
        return -1;

    if (seed == 0) {
        seed = std::random_device{}();
    }
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    clp.print_result();

    wrs::init_threads(num_threads);

    // Generate input weights
    LOG << "Selected " << wrs::generators::select_t::name << " generator";
    auto input_gen = [&normalize, &powerlaw_exp](auto seed, size_t input_size) {
        wrs::timer timer;
        auto weights = wrs::make_numa_arr<double>(input_size);
        wrs::parallel_do_block(
            [data = weights.get(), &seed, &powerlaw_exp](size_t min, size_t max, auto) {
                if (powerlaw_exp > 0) {
                    double exp = -powerlaw_exp;
                    for (size_t i = min; i < max; i++) {
                        data[i] = std::pow(static_cast<double>(i + 1), exp);
                    }
                } else {
                    wrs::generators::select_t generator(seed + min + 1);
                    // TODO size must be even
                    generator.generate_block(data + min, max - min,
                                             /* left_open */ true);
                }
            },
            input_size, static_cast<size_t>(0), 1UL << 18);
        if (powerlaw_exp > 0) {
            // shuffle
            wrs::generators::select_t generator(seed);
            for (size_t i = input_size - 1; i != 0; i--) {
                size_t j = generator.next_int<size_t>(0, i);
                std::swap(weights[i], weights[j]);
            }
        }

        LOG << "generated " << input_size << " input values in "
            << timer.get_and_reset() << "ms";

        // normalize
        if (normalize) {
            double sum =
                std::accumulate(weights.get(), weights.get() + input_size, 0.0);
            sum /= input_size;
            wrs::parallel_do([&weights, &sum](size_t i) { weights[i] /= sum; },
                             input_size);

            LOG << "normalized input in " << timer.get_and_reset() << "ms";
        }

        if (input_size <= 100)
            LOG << std::vector<double>(weights.get(), weights.get() + input_size);

        return weights;
    };

    // clang-format off
    if (!no_construct && !no_seq)
        run_alias<wrs::alias<>>(input_gen, input_size, iterations,
                                repetitions, seed, "sequential");
    if (!no_construct && !no_multi)
        run_alias<wrs::multi_alias<wrs::alias>>(input_gen, input_size, iterations,
                                      repetitions, seed, "multi");
    if (!no_construct && !no_multisweep)
        run_alias<wrs::multi_alias<wrs::simple_scan_alias>>(input_gen, input_size, iterations,
                                      repetitions, seed, "multisweep");
    if (!no_construct && !no_simplescan)
        run_alias<wrs::simple_scan_alias<>>(input_gen, input_size, iterations,
                                            repetitions, seed, "simplescan");
    if (!no_construct && !no_psa)
        run_alias<wrs::par_scan_alias<>>(input_gen, input_size, iterations,
                                         repetitions, seed, "psa");
    if (!no_construct && !no_psa2)
        run_alias<wrs::par_scan_alias2<>>(input_gen, input_size, iterations,
                                          repetitions, seed, "psa2");
    if (!no_construct && !no_outsens)
        run_alias<wrs::par_outsens<>>(input_gen, input_size, iterations,
                                      repetitions, seed, "outsens");
    if (!no_construct && !no_outsensnd)
        run_alias<wrs::par_outsens<int32_t, false, 32>>(
            input_gen, input_size, iterations, repetitions, seed, "outsensnd");
#ifdef WRS_HAVE_GSL
    // For the GSL experiments, for fairness reasons, you need to configure the
    // other algorithms to use size_t as their size_type, too.  Do this by
    // passing size_t as template parameter wherever there are default template
    // arguments in the list above, and appending it to multi_alias'es template
    // argument list. The result should be something like wrs::alias<size_t> or
    // wrs::multi_alias<wrs::scan_alias, size_t>
    if (!no_construct && !no_gsl)
        run_alias<wrs::gsl_alias>(
            input_gen, input_size, iterations, repetitions, seed, "gsl");
#endif

    if (!no_queries && !no_seq)
        run_queries<wrs::alias<>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "sequential");
    if (!no_queries && !no_dedup) {
        run_queries<wrs::dedup<wrs::alias<>>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "seqdedup");

        run_queries<wrs::store_vec<wrs::alias<>>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "seqstore");
    }
    if (!no_queries && !no_multi)
        run_queries<wrs::multi_alias<wrs::alias>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "multi");
    if (!no_queries && !no_multisweep)
        run_queries<wrs::multi_alias<wrs::simple_scan_alias>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "multisweep");
    if (!no_queries && !no_simplescan)
        run_queries<wrs::simple_scan_alias<>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "simplescan");
    if (!no_queries && !no_psa)
        run_queries<wrs::par_scan_alias<>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "psa");
    if (!no_queries && !no_psa2)
        run_queries<wrs::par_scan_alias2<>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "psa2");
    if (!no_queries && !no_outsens)
        run_queries<wrs::par_outsens<>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "outsens");
    if (!no_queries && !no_outsensnd)
        run_queries<wrs::par_outsens<int32_t, false, 128>>(
            input_gen, input_size, iterations, repetitions, queries, seed, "outsensnd");
#ifdef WRS_HAVE_GSL
    if (!no_queries && !no_gsl)
        run_queries<wrs::gsl_alias>(
            input_gen, input_size, iterations, repetitions, queries, seed, "gsl");
#endif
    // clang-format on

    wrs::release_threads();
    return 0;
}
