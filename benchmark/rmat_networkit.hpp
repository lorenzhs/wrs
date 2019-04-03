/*******************************************************************************
 * benchmark/rmat_networkit.hpp
 *
 * NetworKit R-MAT generator
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * NetworKit code:
 *
 * Copyright (C) 2013-2019 NetworKit contributors, MIT License
 * Original Authors: Henning Meyerhenke, Christian Staudt
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_RMAT_NETWORKIT_HEADER
#define BENCHMARK_RMAT_NETWORKIT_HEADER

#include <benchmark/rmat_graph500.hpp>
#include <wrs/aggregate.hpp>

#include <utility>
#include <stdint.h>

template <bool Scramble_IDs>
class rmat_networkit {
public:
    static constexpr const char* name = "networkit";
    static constexpr bool scramble_ids = Scramble_IDs;

    using node = int64_t;
    // definitions from networkit Globals.hpp
    using count = uint64_t;
    using index = uint64_t;

    template <typename RNG>
    rmat_networkit(RNG &rng, int log_n, double a_, double b_, double c_):
        scale(log_n),
        a(a_),
        ab(a_ + b_),
        abc(a_ + b_ + c_),
        scramble_state(graph500::init_scramble_state(rng, log_n))
    {
    }

    template <typename RNG>
    std::pair<node, node> get_edge(size_t /* ignored */, RNG& rng) const {
        auto [src, dst] = drawEdge(rng);
#ifdef RMAT_CLIPFLIP
        if (src > dst) std::swap(src, dst);
#endif // RMAT_CLIPFLIP
        if constexpr (scramble_ids) {
            graph500::scramble_two(src, dst, scramble_state);
        }
        return std::make_pair(src, dst);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback && callback, size_t min, size_t max, RNG& rng) const {
        node src, dst;
        for (size_t i = min; i < max; i++) {
            std::tie(src, dst) = drawEdge(rng);
#ifdef RMAT_CLIPFLIP
            if (src > dst) std::swap(src, dst);
#endif // RMAT_CLIPFLIP
            if constexpr (scramble_ids) {
                graph500::scramble_two(src, dst, scramble_state);
            }
            callback(src, dst);
        }
    }

    size_t table_size() const {
        return 0;
    }

    wrs::Aggregate<double> get_depth_stats() const {
        return wrs::Aggregate<double>().Add(scale);
    }

    wrs::Aggregate<double> get_sample_stats() const {
        return wrs::Aggregate<double>().Add(1);
    }

protected:
    // adapted from NetworKit
    constexpr count quadrant(double r) const {
        if (r <= a) {
            return 0;
        }
        else if (r <= ab) {
            return 1;
        } else if (r <= abc) {
            return 2;
        } else
            return 3;
    }

    // adapted from NetworKit
    template <typename RNG>
    std::pair<node, node> drawEdge(RNG &rng) const {
        node u = 0;
        node v = 0;
        for (index i = 0; i < scale; ++i) {
            count q = quadrant(rng.next());
            u = u << 1;
            v = v << 1;
            u = u | (q >> 1);
            v = v | (q & 1);
        }
        return std::make_pair(u, v);
    }

    count scale;
    double a, ab, abc;

    graph500::scramble_state scramble_state;

};

#endif // BENCHMARK_RMAT_NETWORKIT_HEADER
