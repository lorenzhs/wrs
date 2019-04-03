/*******************************************************************************
 * benchmark/rmat_base.hpp
 *
 * Base functions for R-MAT graph generator using weighted random sampling
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_RMAT_BASE_HEADER
#define BENCHMARK_RMAT_BASE_HEADER

#include <benchmark/rmat_graph500.hpp>

#include <wrs/alias_key.hpp>
#include <wrs/aggregate.hpp>

#include <tlx/logger.hpp>
#include <tlx/math.hpp>

#include <array>
#include <ios>
#include <utility>

template <bool enabled>
struct rmat_stats_t;

template<>
struct rmat_stats_t<false> {
    void record_sample_drawn() { }
    void record_edge_done() { }
    void record_sample_bits(unsigned) { }
    wrs::Aggregate<double> samples_per_edge;
    wrs::Aggregate<double> bits_per_sample;
};

template<>
struct rmat_stats_t<true> {
    void record_sample_drawn() { lvl++; }
    void record_edge_done() {
        samples_per_edge.Add(lvl);
        lvl = 0;
    }
    void record_sample_bits(unsigned bits) {
        bits_per_sample.Add(bits);
    }

    unsigned lvl = 0;
    wrs::Aggregate<double> samples_per_edge;
    wrs::Aggregate<double> bits_per_sample;
};

template <typename Derived>
class rmat_base {
public:
    using node = int64_t;
    using prefix = uint32_t;
    using path = uint64_t;
    using edge = std::pair<node, node>;
    using entry = std::pair<prefix, double>;

    static constexpr unsigned prefix_half_bits = 4 * sizeof(prefix);
    static constexpr unsigned path_half_bits = 4 * sizeof(path);
    static constexpr prefix dst_mask = (prefix{1} << prefix_half_bits) - 1;

    static constexpr bool debug = false;
    static constexpr bool collect_stats = false; // true;
    static constexpr bool scramble_ids = Derived::scramble_ids;

    template <typename RNG>
    rmat_base(RNG &rng, int log_n_, double a_, double b_, double c_) :
        quad({a_, b_, c_, 1.0 - (a_ + b_ + c_)}),
        log_n(log_n_),
        node_mask((path{1} << log_n_) - 1),
        scramble_state(graph500::init_scramble_state(rng, log_n))
    {
        sLOG << "RMAT: need" << 2 * log_n << "path bits";
        sLOG << "RMAT: node extraction mask" << std::hex << node_mask;
    }

    template <typename RNG>
    std::pair<node, node> get_edge(size_t /* ei = ignored */, RNG &rng) const {
        return get_edge(rng);
    }

    template <typename RNG>
    std::pair<node, node> get_edge(RNG &rng) const {
        path result = 0;
        int bits_per_half = 0;

        do {
            prefix next = table.sample(rng.next());
            stats.record_sample_drawn();

            path src_part, dst_part;
            unsigned bits_in_prefix = static_cast<const Derived*>(this)->
                split_prefix(next, src_part, dst_part);

            sLOG << "Prefix" << std::hex << next << "with"
                 << std::dec << bits_in_prefix  << "bits per half; src half:"
                 << std::hex << src_part << "dst:" << dst_part
                 << "old tentative result:" << result;

            bits_per_half += bits_in_prefix;
            result <<= bits_in_prefix;
            result |= (src_part << path_half_bits);
            result |= dst_part;
        } while (bits_per_half < log_n);
        // remove unneeded bits
        unsigned shift = bits_per_half - log_n;
        sLOG << "got enough bits in" << std::hex << result << "-- removing"
             << std::dec << shift << "of" << bits_per_half;
        result >>= shift;
        node src = static_cast<node>(result >> path_half_bits);
        node dst = static_cast<node>(result & node_mask);
        sLOG << "Extracted nodes" << std::hex << src << dst << "from" << result;
        stats.record_edge_done();

        if constexpr (scramble_ids) {
            // scrambles in-place
            graph500::scramble_two(src, dst, scramble_state);
        }

        return std::make_pair(src, dst);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback && callback, size_t min, size_t max, RNG &rng) const {
        get_edges(callback, max - min, rng);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback && callback, size_t num_edges, RNG &rng) const {
        path result = 0;
        int bits_per_half = 0;

        for (size_t i = 0; i < num_edges; /* nothing! */) {
            prefix next = table.sample(rng.next());
            stats.record_sample_drawn();

            path src_part, dst_part;
            unsigned bits_in_prefix = static_cast<const Derived*>(this)->
                split_prefix(next, src_part, dst_part);

            sLOG << "Prefix" << std::hex << next << "with"
                 << std::dec << bits_in_prefix << "bits per half; src half:"
                 << std::hex << src_part << "dst:" << dst_part
                 << "old tentative result:" << result;

            bits_per_half += bits_in_prefix;
            result <<= bits_in_prefix;
            result |= (src_part << path_half_bits);
            result |= dst_part;

            if (bits_per_half >= log_n) {
                sLOG << "got enough bits in" << std::hex << result << "-- have"
                     << std::dec << bits_per_half << "need" << log_n;
                // we have enough bits for an edge, extract it
                int excess = bits_per_half - log_n;
                path tmp = (result >> excess),
                    removal_mask = (path{1} << excess) - 1;
                node src = (tmp >> path_half_bits),
                    dst = tmp & node_mask;

                sLOG << "Extracted nodes" << std::hex << src << dst << "from" << result;
                // implement clip-and-flip
#ifdef RMAT_CLIPFLIP
                if (src > dst) std::swap(src, dst);
#endif // RMAT_CLIPFLIP
                // emit the edge
                if constexpr (scramble_ids) {
                    // scrambles in-place
                    graph500::scramble_two(src, dst, scramble_state);
                }
                stats.record_edge_done();
                callback(src, dst);
                // this is placed here because compilers' optimisers are daft
                removal_mask |= (removal_mask << path_half_bits);
                i++;
                // Reuse the remaining `leftover` bits
                bits_per_half -= log_n;
                result &= removal_mask;
                sLOG << bits_per_half << "bits remaining, tentative result:"
                     << std::hex << result;
            }
        }
    }

    size_t table_size() const {
        return table.size();
    }

    wrs::Aggregate<double> get_depth_stats() const {
        return stats.samples_per_edge;
    }

    wrs::Aggregate<double> get_sample_stats() const {
        return stats.bits_per_sample;
    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    constexpr unsigned split_prefix(const prefix& in, path &out_src, path &out_dst) const {
        assert(in >= (prefix{1} << prefix_half_bits));
        unsigned marker_pos = msb(in);
        unsigned bits_in_prefix = marker_pos - prefix_half_bits;

        assert(bits_in_prefix != static_cast<unsigned>(-1));
        stats.record_sample_bits(bits_in_prefix);

        // need to unset marker bit only for src
        out_src = ((in & ~(1 << marker_pos)) >> prefix_half_bits);
        out_dst = (in & dst_mask);

        assert(out_dst == 0 || msb(out_dst) < bits_in_prefix);

        sLOG << "splitting prefix" << std::hex << in << "marker bit at pos"
             << std::dec << marker_pos << "->" << bits_in_prefix << "bits, src part:"
             << std::hex << out_src << "dst" << out_dst;

        return bits_in_prefix;
    }


protected:
    template <typename integral>
    constexpr unsigned msb(integral i) const {
        return 8 * sizeof(integral) - tlx::clz(i) - 1;
    }

    wrs::alias_key<prefix> table;
    std::array<double, 4> quad;
    const int log_n;
    const path node_mask;

    graph500::scramble_state scramble_state;

    mutable rmat_stats_t<collect_stats> stats;
};

#endif // BENCHMARK_RMAT_BASE_HEADER
