/*******************************************************************************
 * benchmark/rmat.cpp
 *
 * R-MAT graph generator example
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#ifdef WRS_HAVE_NUMA
#undef WRS_HAVE_NUMA
#endif

#include <benchmark/graph_generator.hpp>
#include <benchmark/rmat.hpp>
#include <benchmark/rmat_fd.hpp>
#include <benchmark/rmat_graph500.hpp>
#include <benchmark/rmat_networkit.hpp>
#include <benchmark/util.hpp>
#include <wrs/generators/select.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/timer.hpp>

#include <tlx/cmdline_parser.hpp>
#include <tlx/logger.hpp>

#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <unistd.h>
#include <vector>

static constexpr bool debug = true;

struct arguments {
    size_t num_edges;
    int log_nodes;
    int iterations;
    int repetitions;
    int warmup_its;
    int warmup_reps;
    double par_a, par_b, par_c;
    size_t block_size;
};


template <typename RNG, typename rmat_t, typename initfn>
void run(initfn && init_gen, const arguments &args, size_t seed) {
    LOG << "";
    LOG << "R-MAT generation using " << rmat_t::name << " method";

    RNG gen(seed);
    rmat_t r(gen, args.log_nodes, args.par_a, args.par_b, args.par_c);

    auto init = [&]() {
        init_gen(r)();
    };

    auto runner = [&](bool no_warmup) {
        const bool debug = no_warmup;

        wrs::timer timer;
        graph_generator<rmat_t, RNG> graphgen(r, args.block_size);
        graphgen.get_edges(args.num_edges, seed, debug);

        double duration = timer.get();
        auto depth_stats = r.get_depth_stats(),
            sample_stats = r.get_sample_stats();

        LOG << "RESULT type=sampling"
            << " total=" << duration
            << " edges=" << args.num_edges
            << " logn=" << args.log_nodes
            << " scramble=" << rmat_t::scramble_ids
            << " method=" << rmat_t::name
            << " paths=" << r.table_size()
            << " seed=" << seed
            << " threads=" << wrs::get_num_threads()
            << " depth=" << depth_stats.Avg()
            << " bps=" << sample_stats.Avg();
        return std::vector<double>{ duration, depth_stats.Avg(),
                sample_stats.Avg() };
    };
    auto logger = [&](int it, auto& stats) {
        LOG << "It " << it << " sampling: " << std::fixed
            << std::setprecision(3) << stats[0] << " depth: " << stats[1]
            << " samples: " << stats[2]
            << " seed=" << seed;
        ++seed;
    };
    auto stats = run_benchmark(runner, logger, init,
                               args.iterations, args.repetitions,
                               args.warmup_its, args.warmup_reps);
    LOG << "RMAT: " << std::fixed << std::setprecision(3) << stats[0]
        << " depth: " << stats[1] << " samples: " << stats[2]
        << " (" << args.iterations << " iterations x "
        << args.repetitions << " repetitions)";
}

int main(int argc, char *argv[]) {
    tlx::CmdlineParser clp;

    size_t edges = 1000, capacity = 262144, seed = 0, depth = 9,
        num_threads = 0, block_size = 0;
    int log_n = 30, iterations = 2, repetitions = 5;
    double par_a = 0.57, par_b = 0.19, par_c = 0.19, cutoff = 0;
    bool no_vardepth = false, no_fixdepth = false, no_g500 = false,
        no_networkit = false, warmup = false, scramble = false;
    clp.add_bool('V', "novar", no_vardepth, "don't run variable depth generator");
    clp.add_bool('F', "nofix", no_fixdepth, "don't run fixed depth generator");
    clp.add_bool('G', "nog500", no_g500, "don't run graph500 generator");
    clp.add_bool('N', "nonetworkit", no_networkit, "don't run networkit generator");
    clp.add_bool('w', "warmup", warmup, "run a warmup iteration");
    clp.add_bool('p', "scramble", scramble, "permute node IDs (scramble)");
    clp.add_int('n', "logn", log_n, "log2 of number of nodes");
    clp.add_size_t('m', "edges", edges, "number of edges to generate");
    clp.add_double('a', "par_a", par_a, "R-MAT parameter a");
    clp.add_double('b', "par_b", par_b, "R-MAT parameter b");
    clp.add_double('c', "par_c", par_c, "R-MAT parameter c");
    clp.add_double('l', "cutoff", cutoff, "path probability cutoff");
    clp.add_size_t('x', "capacity", capacity, "maximum number of paths");
    clp.add_size_t('d', "depth", depth, "depth of generated paths (variant 2)");
    clp.add_size_t('s', "seed", seed, "random generator seed (0 for random)");
    clp.add_int('i', "iterations", iterations, "number of iterations");
    clp.add_int('r', "repetitions", repetitions, "repetitions per iteration");
    clp.add_size_t('g', "blocksize", block_size, "block size for parallelisation");
    clp.add_size_t('t', "threads", num_threads, "number of threads");

    if (!clp.process(argc, argv))
        return -1;

    if (par_a < 0 || par_b < 0 || par_c < 0 || par_a + par_b + par_c > 1) {
        LOG1 << "Condition a + b + c < 1, a > 0, b > 0, c > 0 violated";
        return -1;
    }
    if (log_n > 62) {
        LOG1 << "More than 2^62 nodes are not supported";
        return -1;
    }
    if (seed == 0) {
        seed = std::random_device{}();
    }
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    if (block_size == 0) {
        block_size = size_t{1} << 16;
    }
    clp.print_result();

    wrs::init_threads(num_threads);

    using RNG = wrs::generators::select_t;
    LOG << "Selected " << RNG::name << " generator";

    const arguments args { edges, log_n, iterations, repetitions,
                           static_cast<int>(warmup), static_cast<int>(warmup),
                           par_a, par_b, par_c, block_size};

    if (!no_vardepth) {
        auto init = [&](auto &r) {
            return [&]() {
                r.init(cutoff, capacity);
            };
        };
        if (scramble)
            run<RNG, rmat<true>>(init, args, seed);
        else
            run<RNG, rmat<false>>(init, args, seed);
    }

    if (!no_fixdepth) {
        auto init = [&](auto &r) {
            return [&]() {
                r.init(depth);
            };
        };
        if (scramble)
            run<RNG, rmat_fd<true>>(init, args, seed);
        else
            run<RNG, rmat_fd<false>>(init, args, seed);
    }

    if (!no_g500) {
        auto init = [&](auto &) { return [&]() {};  };
        if (scramble)
            run<RNG, rmat_graph500<true>>(init, args, seed);
        else
            run<RNG, rmat_graph500<false>>(init, args, seed);
    }

    if (!no_networkit) {
        auto init = [&](auto &) { return [&]() {};  };
        if (scramble)
            run<RNG, rmat_networkit<true>>(init, args, seed);
        else
            run<RNG, rmat_networkit<false>>(init, args, seed);
    }

    return 0;
}
