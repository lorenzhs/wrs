/*******************************************************************************
 * wrs/multi_alias.hpp
 *
 * Shared memory implementation of the simple distributed algorithm
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_MULTI_ALIAS_HEADER
#define WRS_MULTI_ALIAS_HEADER

#include <wrs/alias.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/timer.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <cassert>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace wrs {

template <template <typename> typename b_alias_type, typename size_type = uint32_t>
struct multi_alias {
    static constexpr const char* name = "multi";
    static constexpr bool yields_single_sample = true;
    static constexpr bool init_with_seed = false;
    static constexpr int pass_rand = 2;
    using alias_type = b_alias_type<size_type>;
    using result_type = size_type;

    static_assert((sizeof(double) / sizeof(size_type)) * sizeof(size_type) == sizeof(double),
                  "size_type has weird size");

    static constexpr bool debug = false;
    static constexpr bool time = false;

    static constexpr double eps = 1e-7;

    multi_alias() {}

    template <typename Iterator>
    multi_alias(Iterator w_begin, Iterator w_end) {
        init(w_end - w_begin);
        construct(w_begin, w_end);
    }

    multi_alias & operator = (const multi_alias &) = default;
    multi_alias(const multi_alias &) = default;
    //! delete move-constructor
    multi_alias(multi_alias &&) = delete;
    //! delete move-assignment
    multi_alias & operator = (multi_alias &&) = delete;

    void init(size_t size) {
        num_subtables_ = get_num_threads();
        size_= size;
        subtables_.resize(num_subtables_);
        subtable_offsets_.resize(num_subtables_);
        subtable_weights_.resize(num_subtables_);

        // Allocate top-level table and subtables
        top_level_.init(num_subtables_);
        parallel_do_range(
            [&](size_t min, size_t max, int thread_id) {
                subtables_[thread_id].init(max - min);
            }, size_);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        timer t;
        if (end - begin != static_cast<ssize_t>(size_)) {
            sLOG1 << "Error: tried to construct alias table of incorrect size!"
                  << "Expected" << size_ << "got" << end - begin;
            return;
        }
        if (size_ > static_cast<size_t>(std::numeric_limits<size_type>::max())) {
            sLOG1 << "Error: size_type cannot hold" << size_
                  << "- please specify an appropriate type (e.g. size_t)";
            return;
        }
        timers_.clear();

        // construct sub-tables
        auto build_subtable = [&](size_t min, size_t max, int thread_id) {
            constexpr bool debug = false;
            sLOG << "Thread" << thread_id << "constructing subtable of size"
            << max - min;
            timer t;
            subtable_offsets_[thread_id] = min;
            subtables_[thread_id].construct(begin + min, begin + max);
            sLOG << "Thread" << thread_id << "took" << t.get() << "ms";
        };
        parallel_do_range(build_subtable, size_);

        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: sub-tables built in " << timers_.back() << "ms";

        // construct top-level table
        for (size_t i = 0; i < num_subtables_; i++) {
            subtable_weights_[i] = subtables_[i].total_weight();
        }
        top_level_.construct(subtable_weights_.begin(),
                             subtable_weights_.end());
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 2: top-level table construction took "
                   << timers_.back() << "ms";
    }

    size_type sample(double uniform1, double uniform2) const {
        size_type subtable = top_level_.sample(uniform1);
        size_t offset = subtable_offsets_[subtable];
        assert(subtables_[subtable].size() != static_cast<size_t>(-1));
        return offset + subtables_[subtable].sample(uniform2);
    }

    // Sum up all weights in the subtables
    void verify_helper(std::vector<double> &weights, size_t offset) const {
        for (size_t i = 0; i < num_subtables_; ++i) {
            subtables_[i].verify_helper(weights, offset + subtable_offsets_[i]);
        }
    }

    size_t size() const {
        return size_;
    }

    std::vector<double> get_timers() const {
        return timers_;
    }

    double total_weight() const {
        return top_level_.total_weight();
    }


    template <typename Callback>
    void find(size_type item, Callback && callback) const {
        for (size_t table = 0; table < num_subtables_; ++table) {
            size_type offset = subtable_offsets_[table];
            if (item < offset) continue;

            auto cb =
                [&](int, size_type idx, double w, auto &tableitem) {
                    callback(table, idx + offset, w, tableitem);
                };
            subtables_[table].find(item - offset, cb);

            if (table + 1 < num_subtables_ && item < subtable_offsets_[table + 1])
                break;
        }
    }

protected:
    size_t num_subtables_, size_;
    alias_type top_level_;
    std::vector<alias_type> subtables_;
    std::vector<size_t> subtable_offsets_;
    std::vector<double> subtable_weights_;
    std::vector<double> timers_;
};

} // namespace wrs

#endif // WRS_MULTI_ALIAS_HEADER
