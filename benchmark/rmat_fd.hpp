/*******************************************************************************
 * benchmark/rmat.hpp
 *
 * R-MAT graph generator using weighted random sampling
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_RMAT_FD_HEADER
#define BENCHMARK_RMAT_FD_HEADER

#include <benchmark/rmat_base.hpp>

#include <tlx/logger.hpp>

#include <queue>
#include <utility>
#include <vector>


template <bool Scramble_IDs = true>
struct rmat_fd : public rmat_base<rmat_fd<Scramble_IDs>> {
    static constexpr const char* name = "fixdepth";
    using supertype = rmat_base<rmat_fd<Scramble_IDs>>;
    friend supertype;
    using entry = typename supertype::entry;
    using path = typename supertype::path;
    using prefix = typename supertype::prefix;
    using supertype::prefix_half_bits;
    using supertype::log_n;
    using supertype::dst_mask;
    using supertype::table;
    using supertype::quad;
    using supertype::stats;
    static constexpr bool scramble_ids = Scramble_IDs;

    static constexpr bool debug = false;

    template <typename RNG>
    rmat_fd(RNG &rng_, int log_n_, double a_, double b_, double c_) :
        supertype(rng_, log_n_, a_, b_, c_) {}

    void init(int max_depth) {
        wrs::timer t;

        table_depth = max_depth;
        marker_pos = prefix_half_bits + table_depth;
        marker_removal_mask = ~(prefix{1} << marker_pos);

        entries.clear();

        enum_items(0, 1.0, 0);

        sLOG1 << "generated" << entries.size() << "path prefixes in"
              << t.get_and_reset() << "ms";

        table.init(entries.size());
        sLOG1 << "init table in" << t.get_and_reset() << "ms";

        table.construct(entries.begin(), entries.end(), /* is_dist */ true);
        sLOG1 << "construct table in" << t.get_and_reset() << "ms";

    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    constexpr unsigned split_prefix(const prefix &in, path &out_src, path &out_dst) const {
        assert(in != 0 && in < (path{1} << (table_depth + prefix_half_bits)));
        // We know how many bits are in each prefix, no need to calculate it
        out_dst = (in & dst_mask);
        assert(out_dst < (path{1} << table_depth));
        out_src = (in >> prefix_half_bits);
        stats.record_sample_bits(table_depth);
        return table_depth;
    }

protected:
    void enum_items(prefix pref, double prob, int depth) {
        for (size_t d = 0; d < quad.size(); ++d) {
            prefix new_prefix = pref;
            unsigned dst_pos = depth;
            unsigned src_pos = depth + prefix_half_bits;

            sLOG << "Considering bits" << src_pos << "and" << dst_pos
                 << "of" << pref
                 << "src half:" << (new_prefix >> prefix_half_bits)
                 << "dst half:" << (new_prefix & dst_mask)
                 << "direction" << d;

            // Then add this position and the new marker bit
            new_prefix |= ((d >> 1) << src_pos);
            new_prefix |= ((d & 1) << dst_pos);

            sLOG << "\tresult:" << new_prefix << "src half:"
                 << (new_prefix >> prefix_half_bits)
                 << "dst half:" << (new_prefix & dst_mask);

            const double new_prob = prob * quad[d];
            if (depth + 1 < table_depth) {
                enum_items(new_prefix, new_prob, depth + 1);
            } else {
                entries.emplace_back(new_prefix, new_prob);
            }
        }
    }

    std::vector<entry> entries;
    unsigned marker_pos;
    prefix marker_removal_mask;
    int table_depth;
};

#endif // BENCHMARK_RMAT_FD_HEADER
