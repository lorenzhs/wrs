/*******************************************************************************
 * benchmark/alias.cpp
 *
 * Weighted Random Permutation via Alias methods
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <benchmark/util.hpp>
#include <wrs/alias.hpp>
#include <wrs/aggregate.hpp>
#include <wrs/generators/select.hpp>
#include <wrs/multi_alias.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/par_alias.hpp>
#include <wrs/timer.hpp>

// GSL is optional, only include it if it was found
#ifdef WRS_HAVE_GSL
#include <wrs/gsl.hpp>
#endif

#include <tlx/cmdline_parser.hpp>
#include <tlx/logger.hpp>
#include <tlx/thread_pool.hpp>

#include <iomanip>
#include <vector>

#ifdef WRS_HAVE_NUMA
static constexpr bool numa = true;
#else
static constexpr bool numa = false;
#endif

static constexpr bool debug = true;

template <typename alias_table, typename input_generator>
wrs::Aggregate<double> run_alias(
    input_generator && input_gen, size_t input_size,
    int iterations, int repetitions,
    size_t seed, const std::string &name)
{
    LOG << "";
    LOG << "Constructing alias table using " << name << " method";

    wrs::numa_arr_ptr<double> weights;
    alias_table table;

    auto init = [&]() {
        weights = input_gen(seed, input_size);
        table.init(input_size);
    };
    auto runner = [&](bool no_warmup) {
        const bool debug = no_warmup;
        wrs::timer timer;
        table.construct(weights.get(), weights.get() + input_size);
        double duration = timer.get();
        table.verify(weights.get(), weights.get() + input_size);
        auto timers = table.get_timers();

        LOG << "RESULT type=construction"
            << " total=" << duration
            << " timers=" << timers
            << " size=" << input_size
            << " method=" << name
            << " threads=" << wrs::get_num_threads() // only relevant for parallel...
            << " numa=" << numa
            << " seed=" << seed;

        timers.insert(timers.begin(), duration);
        return timers;
    };
    auto logger = [&](int it, auto& stats) {
        LOG << "It " << it << " construction: " << std::fixed
        << std::setprecision(3) << stats[0] << " seed=" << seed;
        ++seed;
    };
    auto stats = run_benchmark(runner, logger, init, iterations, repetitions,
                               /* warmup its */ 1, /* warmup reps */ 1);
    LOG << "Construction: " << std::fixed << std::setprecision(3) << stats[0]
        << " (" << iterations << " iterations x " << repetitions << " repetitions)";
    if (stats.size() == 9) { // parallel
        LOG << "\tPreprocessing:  " << stats[1];
        LOG << "\tClassification: " << stats[2];
        LOG << "\tB prefix sum:   " << stats[3];
        LOG << "\tAssign small:   " << stats[4];
        LOG << "\tZero F and T:   " << stats[5];
        LOG << "\tF prefix sum:   " << stats[6];
        LOG << "\tL prefix sum:   " << stats[7];
        LOG << "\tAssign large:   " << stats[8];
    } else if (stats.size() == 3) { // multi
        LOG << "\tSubtable constr: " << stats[1];
        LOG << "\tTop-level table: " << stats[2];
    } else if (stats.size() == 4) { // sequential
        LOG << "\tPreprocessing:  " << stats[1];
        LOG << "\tClassification: " << stats[2];
        LOG << "\tAssignment:     " << stats[3];
    } else { // GSL
        LOG << "\tConstruction:   " << stats[1];
    }
    return stats[0];
}

template <typename alias_table, typename input_generator>
wrs::Aggregate<double> run_queries(
    input_generator && input_gen, size_t input_size,
    int iterations, int repetitions,
    size_t queries, size_t seed, bool duplicate, const std::string &name)
{
    LOG << "";
    LOG << "Querying alias table using " << name << " method";

    // construct table
    wrs::numa_arr_ptr<double> weights;
    alias_table table;

    auto init = [&]() {
        weights = input_gen(seed, input_size);
        table.init(input_size);
        table.construct(weights.get(), weights.get() + input_size);
    };

    auto runner = [&](const bool debug) {
        wrs::timer timer;

        // Copy tables onto each NUMA node
        std::vector<std::unique_ptr<alias_table>> tables(wrs::g_num_numa_nodes);
        // Only duplicate if >1 NUMA nodes
        if (duplicate && wrs::g_num_numa_nodes > 1) {
            wrs::do_per_numa_node(
                [&tables, &table](int id) {
                    tables[id] = std::make_unique<alias_table>(table);
                });
        }

        double copy_duration = timer.get_and_reset();

        // Run queries in parallel
        if (duplicate && wrs::g_num_numa_nodes > 1)
            wrs::parallel_do_range(
                [&tables, &seed](size_t min, size_t max, int thread_id) {
                    const int numa_node = thread_id /
                        (wrs::g_total_threads + wrs::g_num_numa_nodes - 1)
                        / wrs::g_num_numa_nodes;

                    wrs::generators::select_t rng(seed + min + 1);
                    std::vector<typename alias_table::result_type> samples(1);
                    for (size_t i = min; i < max; i++) {
                        if constexpr (alias_table::pass_rand) {
                            samples[0] = tables[numa_node]->sample(rng.next());
                        } else {
                            samples[0] = tables[numa_node]->sample();
                        }
                    }
                }, queries);
        else
            wrs::parallel_do_range(
                [&table, &seed](size_t min, size_t max, auto) {
                    wrs::generators::select_t rng(seed + min + 1);
                    std::vector<typename alias_table::result_type> samples(1);
                    for (size_t i = min; i < max; i++) {
                        if constexpr (alias_table::pass_rand) {
                            samples[0] = table.sample(rng.next());
                        } else {
                            samples[0] = table.sample();
                        }
                    }
                }, queries);


        double duration = timer.get();

        LOG << "RESULT type=query"
            << " total=" << duration
            << " tcopy=" << copy_duration
            << " size=" << input_size
            << " queries=" << queries
            << " method=" << name
            << " threads=" << wrs::get_num_threads()
            << " numa=" << numa
            << " seed=" << seed;

        return std::vector<double>{duration};
    };
    auto logger = [&](int it, auto& stats) {
        LOG << "It " << it << " queries: " << std::fixed
        << std::setprecision(3) << stats[0] << " seed=" << seed;
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
    bool normalize = false, no_seq = false, no_par = false, no_multi = false,
        no_gsl = false,
        no_queries = false, no_construct = false, no_duplicate = false;
    clp.add_size_t('n', "input-size", input_size, "input size");
    clp.add_int('i', "iterations", iterations, "iterations");
    clp.add_int('r', "repetitions", repetitions, "repetitions per iteration");
    clp.add_size_t('q', "queries", queries, "number of queries");
    clp.add_size_t('s', "seed", seed, "random generator seed (0 for random)");
    clp.add_bool('z', "normalize", normalize, "normalize inputs (for debuggability");
    clp.add_bool('D', "nodup", no_duplicate,
                 "don't duplicate tables onto every NUMA node for queries");
    clp.add_bool('Q', "noqueries", no_queries, "don't run query benchmark");
    clp.add_bool('C', "noconstruct", no_construct, "don't run construction benchmark");
    clp.add_bool('S', "noseq", no_seq, "don't run sequential version");
    clp.add_bool('P', "nopar", no_par, "don't run parallel version");
    clp.add_bool('M', "nomulti", no_multi, "don't run multi version");
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
    auto input_gen = [&normalize](auto seed, size_t input_size) {
        wrs::timer timer;
        auto weights = wrs::make_numa_arr<double>(input_size);
        wrs::parallel_do_range(
            [data = weights.get(), &seed](size_t min, size_t max, auto) {
                wrs::generators::select_t generator(seed + min + 1);
                // TODO size must be even
                generator.generate_block(data + min, max - min);
            }, input_size);

        LOG << "generated " << input_size << " input values in "
        << timer.get_and_reset() << "ms";

        // normalize
        if (normalize) {
            double sum = std::accumulate(weights.get(),
                                         weights.get() + input_size, 0.0);
            sum /= input_size;
            wrs::parallel_do(
                [&weights, &sum](size_t i) { weights[i] /= sum; }, input_size);

            LOG << "normalized input in " << timer.get_and_reset() << "ms";
        }

        if (input_size <= 100)
            LOG << std::vector<double>(weights.get(), weights.get() + input_size);

        return weights;
    };

    if (!no_construct && !no_par)
        run_alias<wrs::par_alias<size_t>>(input_gen, input_size, iterations,
                                    repetitions, seed, "parallel");
    if (!no_construct && !no_seq)
        run_alias<wrs::alias<size_t>>(input_gen, input_size, iterations,
                                repetitions, seed, "sequential");
    if (!no_construct && !no_multi)
        run_alias<wrs::multi_alias<size_t>>(input_gen, input_size, iterations,
                                      repetitions, seed, "multi");
#ifdef WRS_HAVE_GSL
    if (!no_construct && !no_gsl)
        run_alias<wrs::gsl_alias>(input_gen, input_size, iterations,
                                  repetitions, seed, "gsl");
#endif

    if (!no_queries && !no_par)
        run_queries<wrs::par_alias<size_t>>(
            input_gen, input_size, iterations, repetitions, queries, seed,
            !no_duplicate, "parallel");
    if (!no_queries && !no_seq)
        run_queries<wrs::alias<size_t>>(
            input_gen, input_size, iterations, repetitions, queries, seed,
            !no_duplicate, "sequential");
    if (!no_queries && !no_multi)
        run_queries<wrs::multi_alias<size_t>>(
            input_gen, input_size, iterations, repetitions, queries, seed,
            !no_duplicate, "multi");
#ifdef WRS_HAVE_GSL
    if (!no_queries && !no_gsl)
        run_queries<wrs::gsl_alias>(
            input_gen, input_size, iterations, repetitions, queries, seed,
            !no_duplicate, "gsl");
#endif

    return 0;
}
