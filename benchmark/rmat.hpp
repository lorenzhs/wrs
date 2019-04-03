/*******************************************************************************
 * benchmark/rmat.hpp
 *
 * R-MAT graph generator using weighted random sampling
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_RMAT_HEADER
#define BENCHMARK_RMAT_HEADER

#include <benchmark/rmat_base.hpp>

#include <tlx/logger/all.hpp>

#include <limits>
#include <queue>
#include <utility>
#include <vector>


template <bool Scramble_IDs = true>
class rmat : public rmat_base<rmat<Scramble_IDs>> {
public:
    static constexpr const char* name = "vardepth";
    using supertype = rmat_base<rmat<Scramble_IDs>>;
    friend supertype;
    using entry = typename supertype::entry;
    using prefix = typename supertype::prefix;
    using supertype::prefix_half_bits;
    using supertype::dst_mask;
    using supertype::table;
    using supertype::quad;
    static constexpr bool scramble_ids = Scramble_IDs;

    template <typename T>
    struct entry_less {
        constexpr bool operator()(const T &lhs, const T &rhs) const {
            return lhs.second < rhs.second;
        }
    };

    using queue_t = std::priority_queue<entry, std::vector<entry>,
                                        entry_less<entry>>;

    static constexpr bool debug = false;

    template <typename RNG>
    rmat(RNG &rng_, int log_n_, double a_, double b_, double c_) :
        supertype(rng_, log_n_, a_, b_, c_) {}

    void init(double cutoff, size_t max_size) {
        wrs::timer t;

        entries.clear();
        assert(queue.empty());

        queue.emplace((1 << prefix_half_bits), 1.0);
        while (!queue.empty() && entries.size() + queue.size() + 3 < max_size) {
            entry curr = queue.top();
            queue.pop();

            for (size_t d = 0; d < quad.size(); ++d) {
                prefix new_prefix = curr.first << 1;
                new_prefix |= ((d >> 1) << prefix_half_bits);
                new_prefix |= (d & 1);
                double prob = curr.second * quad[d];

                sLOG << "Got" << std::hex << curr.first << "with src half:"
                     << (curr.first >> prefix_half_bits) << "dst half:"
                     << (curr.first & dst_mask) << "for direction" << d;
                sLOG << "\tresult:" << std::hex << new_prefix
                     << "src half:" << (new_prefix >> prefix_half_bits)
                     << "dst half:" << (new_prefix & dst_mask)
                     << "probability" << prob;

                if (prob <= cutoff ||
                    new_prefix >= (std::numeric_limits<prefix>::max() >> 1)) {
                    entries.emplace_back(new_prefix, prob);
                } else {
                    queue.emplace(new_prefix, prob);
                }
            }
        }
        while (!queue.empty()) {
            entries.push_back(queue.top());
            queue.pop();
        }


        // Sort entries by probability to increase locality of table
        std::sort(entries.begin(), entries.end(), entry_less<entry>());

        sLOG1 << "generated" << entries.size() << "path prefixes in"
              << t.get_and_reset() << "ms";

        table.init(entries.size());
        sLOG1 << "init table in" << t.get_and_reset() << "ms";

        table.construct(entries.begin(), entries.end(), /* is_dist */ true);
        sLOG1 << "construct table in" << t.get_and_reset() << "ms";
    }
protected:
    std::vector<entry> entries;
    queue_t queue;
};

#endif // BENCHMARK_RMAT_HEADER
