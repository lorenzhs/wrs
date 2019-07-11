/*******************************************************************************
 * wrs/psa/psa_base.hpp
 *
 * Common pieces for parallel 1-alias table construction algorithms
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PSA_PSA_BASE_HEADER
#define WRS_PSA_PSA_BASE_HEADER

#include <wrs/psa/subproblem.hpp>

#include <wrs/numa_array.hpp>
#include <wrs/util.hpp>
#include <wrs/verify.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace wrs {
namespace psa {

template <template <typename> typename Derived, typename size_type>
class psa_base {
protected:
    static constexpr bool debug = false;
    static constexpr bool time = false;

    // relative epsilon, will be multiplied with number of items
    static constexpr double relative_epsilon = 1e-10;

public:
    using derived = Derived<size_type>;
    static constexpr const char* name = derived::name;

    static constexpr bool yields_single_sample = true;
    static constexpr bool init_with_seed = false;
    static constexpr int pass_rand = 1;
    using result_type = size_type;
    using Subproblem = subproblem<size_type>;
    using Splitter = splitter<size_type>;
    static constexpr size_type empty = static_cast<size_type>(-1);


    struct tableitem {
        double p; // probability
        size_type a; // alias
        tableitem() : p(0), a(static_cast<size_type>(-1)) {}
        tableitem(double p_, size_type a_) : p(p_), a(a_) {}
    };
    friend std::ostream &operator << (std::ostream &os, const tableitem &item) {
        return os << '(' << item.p << ',' << item.a << ')';
    }


    psa_base() : table_(nullptr), indices_(nullptr), classifier_(nullptr)
               , size_(-1), num_light_(-1), W_(0), W_n_(0), eps_(0) {}

    psa_base & operator = (const psa_base & other) {
        LOG0 << "psa copy assignment/constructor called";
        size_ = other.size_;
        num_light_ = other.num_light_;
        W_ = other.W_;
        W_n_ = other.W_n_;
        eps_ = other.eps_;
        table_ = make_numa_arr<tableitem>(size_);
        memcpy(table_.get(), other.table_.get(), size_ * sizeof(tableitem));
        indices_ = make_numa_arr<size_type>(size_);
        memcpy(indices_.get(), other.indices_.get(), size_ * sizeof(size_type));
        prefsum_ = make_numa_arr<double>(size_);
        memcpy(prefsum_.get(), other.prefsum_.get(), size_ * sizeof(double));
        classifier_ = make_numa_arr<size_type>(size_);
        subproblems_ = std::make_unique<Subproblem[]>(wrs::get_num_threads());
        splitters_ = std::make_unique<Splitter[]>(wrs::get_num_threads() + 1);
        timers_ = other.timers_;
        return *this;
    }
    psa_base(const psa_base &other) {
        *this = other;
    }
    //! delete move-constructor
    psa_base(psa_base &&) = delete;
    //! delete move-assignment
    psa_base & operator = (psa_base &&) = delete;

    void init(size_t /* not size_type */ size) {
        if (size > static_cast<size_t>(std::numeric_limits<size_type>::max())) {
            sLOG1 << "Error: size_type cannot hold" << size_
                  << "- please specify an appropriate type (e.g. size_t)";
            return;
        }
        size_ = size;

        table_ = make_numa_arr<tableitem>(size_);
        indices_ = make_numa_arr<size_type>(size_);
        prefsum_ = make_numa_arr<double>(size_);
        classifier_ = make_numa_arr<size_type>(size_);
        subproblems_ = std::make_unique<Subproblem[]>(wrs::get_num_threads());
        splitters_ = std::make_unique<Splitter[]>(wrs::get_num_threads() + 1);
    }


    // Query given a uniformly distributed [0,1) random value
    size_type sample(double uniform) const {
        double rand = uniform * size_;
        size_t index = rand;
        tableitem& item = table_[index];

        // sanity check: if alias == -1, then prob == 1
        assert(item.a != empty || std::abs(item.p - W_n_) < 1e-7);

        rand = (rand - index) * W_n_;
        assert(rand < item.p || item.a != empty);
        return (rand < item.p) ? index : item.a;
    }


    // Sum up all weights in the table, adding `offset` to every item ID.
    // This is used internally by verify() and by multi_alias' verify()
    void verify_helper(std::vector<double> &weights, size_t offset) const {
        wrs::sum_table<size_type>(
            table_.get(), size_, W_, W_n_, weights, offset);
    }

    std::vector<double> get_timers() const {
        return timers_;
    }

    size_t size() const {
        return size_;
    }

    double total_weight() const {
        return W_;
    }

    template <typename Callback>
    void find(size_type item, Callback && callback) const {
        if (size_ == static_cast<size_type>(-1)) return;
        callback(0, item, table_[item].p, table_[item]);
        for (size_type i = 0; i < size_; i++) {
            if (table_[i].a == item) {
                callback(0, i, W_n_ - table_[i].p, table_[i]);
            }
        }
    }


protected:
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool is_light(const double& d) const {
        return d <= W_n_;
    }

    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool is_light(tableitem& i) const {
        return is_light(i.p);
    }

    std::pair<size_type, double>
    get_light(size_type local_index, const Subproblem &p) const {
        assert(local_index < p.num_light);
        size_type idx = indices_[p.l_begin + local_index];
        return std::make_pair(idx, table_[idx].p);
    }

    std::pair<size_type, double>
    get_heavy(size_type local_index, const Subproblem &p) const {
        // prefetch indices array because boundaries are unlikely
        size_type idx = p.h_begin + local_index;
        __builtin_prefetch(&indices_[idx], 0, 1);

        if (TLX_UNLIKELY(p.has_left_boundary && local_index == 0)) {
            sLOG << "get_heavy left boundary:" << indices_[p.h_begin - 1];
            return std::make_pair(
                indices_[p.h_begin], p.left_boundary_weight);
        } else if (TLX_UNLIKELY(p.has_right_boundary &&
                                local_index == p.num_heavy)) {
            LOG << "get_heavy right boundary: local idx " << local_index;
            return std::make_pair(
                indices_[p.h_end - 1], p.right_boundary_weight);
        }
        assert(local_index < (p.h_end - p.h_begin));
        idx = indices_[idx];
        sLOG << "get_heavy mid: local" << local_index << "idx" << idx;
        return std::make_pair(idx, table_[idx].p);
    }

    void set_heavy(double weight, size_type local_index, Subproblem &p) {
        // no prefetching needed here, because item is already in cache
        if (TLX_UNLIKELY(p.has_left_boundary && local_index == 0)) {
            sLOG << "set_heavy setting left boundary weight to" << weight;
            p.left_boundary_weight = weight;
            // left boundary item is always owned, so don't return here
        } else if (TLX_UNLIKELY(p.has_right_boundary &&
                                local_index == p.num_heavy)) {
            sLOG << "set_heavy setting right boundary weight to" << weight;
            p.right_boundary_weight = weight;
            return;
        }
        assert(local_index < p.num_heavy ||
               (p.has_left_boundary && local_index == 0));
        size_type idx = indices_[p.h_begin + local_index];
        sLOG << "set_heavy setting bucket" << idx << "weight to" << weight;
        assert(weight > 0);
        table_[idx].p = weight;
    }

    std::pair<size_type, double>
    get_item(size_type local_index, const Subproblem &p) const {
        LOG << "get_item(" << local_index << ", " << p << ")";
        if (local_index < p.num_light)
            return get_light(local_index, p);
        else
            return get_heavy(local_index - p.num_light, p);
    }


    Splitter split(const Subproblem &problem, double fraction) const {
        LOG << "\nsplit(" << problem << ", " << fraction << ")";
        assert(0 < fraction && fraction < 1);
        /* yikes */
        double first_light_w = problem.l_begin == 0 ? 0 :
            prefsum_[problem.l_begin - 1];
        double last_light_w = problem.l_end == 0 ? 0 :
            prefsum_[problem.l_end - 1];
        double first_heavy_w = problem.h_begin == num_light_ ? 0 :
            prefsum_[problem.h_begin - 1];
        double last_heavy_w = problem.h_end == 0 ? 0 :
            prefsum_[problem.h_end - 1];

        size_type num_left = static_cast<size_type>(problem.size * fraction),
            num_right = problem.size - num_left;
        double target_weight = num_left * W_n_;
        sLOG << "split() targeting fraction" << fraction << "of" << problem.size
             << "=" << problem.size * fraction << "rounded to" << num_left
             << "target weight" << target_weight << "with W/n =" << W_n_;

        size_type left_heavy_idx = 0, right_heavy_idx = problem.num_heavy - 1,
                  heavy_split_idx = 0;
        double split_weight;
        while (left_heavy_idx < right_heavy_idx) {
            heavy_split_idx = (right_heavy_idx + left_heavy_idx) / 2;
            size_type light_split_idx = std::max<size_type>(0, num_left - heavy_split_idx - 1);
            light_split_idx = std::min(light_split_idx, problem.num_light - 1);

            double left_light_weight = light_split_idx < 0 ? 0.0 :
                prefsum_[problem.l_begin + light_split_idx] - first_light_w;
            double left_heavy_weight =
                prefsum_[problem.h_begin + heavy_split_idx] - first_heavy_w;

            split_weight = left_light_weight + left_heavy_weight;

            sLOG << "binsearch: L =" << left_heavy_idx << "mid = " << heavy_split_idx
                 << "R =" << right_heavy_idx << "split weight =" << split_weight
                 << "=" << left_light_weight << "+" << left_heavy_weight
                 << "vs target of" << target_weight;

            if (split_weight > target_weight) {
                right_heavy_idx = heavy_split_idx;
            } else {
                heavy_split_idx++;
                left_heavy_idx = heavy_split_idx;
            }
        }

        size_type heavy_left = heavy_split_idx,
            heavy_right = problem.num_heavy - heavy_left;
        sLOG << "found split target, have" << heavy_left << "heavy items"
             << "on left, 1 split, and" << heavy_right << "on right side";

        size_type light_left = num_left - heavy_left,
            light_right = problem.num_light - light_left;

        assert(light_left <= num_light_);
        assert(light_right <= num_light_);

        // Compute weights of subproblems.  While this looks redundant to the
        // computations of the binary search above, it is not (the split index
        // can be incremented after the last recomputation in the binary search)
        double light_split_w = light_left == 0 ? first_light_w :
                prefsum_[problem.l_begin + light_left - 1],
            // without (w1) and with (w2) split item, to properly exclude it
            heavy_split_w1 = heavy_left == 0 ? first_heavy_w :
                prefsum_[problem.h_begin + heavy_left - 1],
            heavy_split_w2 = heavy_left == 0 ? 0 :
                prefsum_[problem.h_begin + heavy_left];
        double left_light_weight = light_split_w - first_light_w;
        double left_heavy_weight = heavy_split_w1 - first_heavy_w;
        double right_light_weight = last_light_w - light_split_w;
        double right_heavy_weight = last_heavy_w - heavy_split_w2;
        sLOG << "split: left" << num_left << "items are" << light_left
             << "light with weight" << left_light_weight << "and" << heavy_left
             << "heavy with weight" << left_heavy_weight;
        sLOG << "split: right" << num_right << "items are" << light_right
             << "light with weight" << right_light_weight << "and"
             << heavy_right << "heavy with weight" << right_heavy_weight
             << "(+ split item)";

        // Figure out splitting of boundary item
        auto [heavy_split_item, heavy_split_weight] = get_heavy(heavy_split_idx, problem);
        LOG << "get_heavy(" << heavy_split_idx << ") returned ("
            << heavy_split_item << ", " << heavy_split_weight << ")";
        double left_partial_weight = num_left * W_n_ - (left_light_weight + left_heavy_weight);
        double right_partial_weight = heavy_split_weight - left_partial_weight;
        bool right_owns_split = true;
        if (right_partial_weight > num_right * W_n_) {
            sLOG << "Capping right partial weight" << right_partial_weight
                 << "at num_right_ * W_n_ =" << num_right * W_n_;
            right_partial_weight = num_right * W_n_;
            right_owns_split = false;
        }

        sLOG << "left regular weight:" << left_light_weight + left_heavy_weight
             << "of" << num_left * W_n_ << "getting" << left_partial_weight
             << "right:" << right_light_weight + right_heavy_weight << "of"
             << num_right * W_n_ << "getting" << right_partial_weight;
        assert(left_partial_weight >= 0);
        assert(right_partial_weight >= 0);

        Splitter result =
        {
            problem.l_begin + light_left, // l_end
            problem.h_begin + heavy_left, // h_split
            left_partial_weight, right_partial_weight,
            (left_partial_weight > 1e-10), // has_split
            right_owns_split
        };
        LOG << "splitter: " << result;

        return result;
    }

    void verify_subproblem_feasibility(const Subproblem &problem) const {
        (void)problem;
#ifndef NDEBUG
        assert(problem.sanity_check());
        assert(problem.size == problem.num_light + problem.num_heavy);
        double sum = 0;
        for (size_type i = 0; i < problem.num_light; i++) {
            auto [idx, weight] = get_light(i, problem);
            sLOG << "light item" << i << "idx" << idx << "weight" << weight;
            assert(is_light(weight));
            sum += weight;
        }
        size_type max_heavy = problem.num_heavy;
        max_heavy += problem.has_right_boundary;

        for (size_type i = 0; i < max_heavy; i++) {
            auto [idx, weight] = get_heavy(i, problem);
            sLOG << "heavy item" << i << "idx" << idx << "weight" << weight;
            assert(!is_light(weight) || (i == 0 && problem.has_left_boundary)
                   || (i == problem.num_heavy && problem.has_right_boundary));
            sum += weight;
        }
        double exp_sum = problem.size * W_n_;
        double diff = exp_sum - sum;
        assert(std::abs(diff) < 1e-10 * std::max(1.0, sum));
#endif
    }

    void base_case(Subproblem &problem) {
        sLOG << "Base case has" << problem.num_light << "light items,"
             << problem.num_heavy << "heavy ones"
             << (problem.has_right_boundary ? "plus one shared right item" :
                 "no right boundary");

        // fill heavy items
        size_type light_idx = 0, heavy_idx = 0;
        const size_type max_heavy_idx = problem.num_heavy + problem.has_right_boundary;
        while (heavy_idx < max_heavy_idx && light_idx < problem.size) {
            auto [heavy_id, heavy_weight] = get_heavy(heavy_idx, problem);
            sLOG << "now considering heavy item" << heavy_id << "local idx"
                 << heavy_idx << "remaining weight" << heavy_weight;
            //tableitem &heavy = table[heavy_id];
            while (light_idx < problem.size &&
                   (heavy_weight > W_n_ ||
                    (problem.has_right_boundary && heavy_idx == problem.num_heavy)))
            {
                auto [light_id, light_weight] = get_item(light_idx, problem);
                if (TLX_UNLIKELY(light_id == heavy_id)) {
                    LOG << "light_id == heavy_id == " << light_id
                        << " -- skipping to next light item";
                    light_idx++;
                    continue;
                }
                tableitem &light = table_[light_id];
                assert(light.a == empty);
                light.a = heavy_id;
                assert(light.p <= W_n_ + eps_);
                heavy_weight += (light.p - W_n_);

                // we can prefetch the next item right away, surprisingly, this
                // improves base case performance by a few percent
                __builtin_prefetch(&light + 1, 1);

                set_heavy(heavy_weight, heavy_idx, problem);
                sLOG << "filling bucket" << light_id << "local index"
                     << light_idx << "with" << W_n_ - light.p << "of heavy item" << heavy_id
                     << "local index" << heavy_idx << "remainder" << heavy_weight
                     << table_[heavy_id]
                     << "bucket" << light;
                assert(heavy_weight >= -eps_);
                light_idx++;
            }
            sLOG << "inner loop ended with light_idx" << light_idx
                 << "heavy_idx" << heavy_idx << "heavy id" << heavy_id
                 << "weight" << heavy_weight << (heavy_weight > W_n_);
            heavy_idx++;
        }
    }


    Subproblem subproblem_from_splitters(const Splitter &left,
                                         const Splitter &right) const
    {
        sLOG << "combining" << left << "and" << right;

        // there can be valid subproblems without a light item if the part of the
        // split left item is light
        assert(left.l_end <= right.l_end);
        assert(left.h_split <= right.h_split);
        size_type l_begin = left.l_end, l_end = right.l_end,
            h_begin = left.h_split,
            h_end = right.h_split + right.has_split_item;
        size_type num_light = l_end - l_begin,
            num_heavy = right.h_split - left.h_split,
            size = num_light + num_heavy;
        double left_boundary_weight = left.split_weight_right,
            right_boundary_weight = right.split_weight_left;
        bool has_left_boundary = left.has_split_item,
            has_right_boundary = right.has_split_item;
        // condition on left.has_split_item to avoid false positives at the start
        if (left.h_split == right.h_split && left.has_split_item) {
            // we only have a part of one heavy item
            sLOG << "only have a partial heavy item! Setting as right boundary";
            left_boundary_weight = 0;
            right_boundary_weight = left.split_weight_right - right.split_weight_right;
            has_left_boundary = false;
            has_right_boundary = true;
            sLOG << "settled on taking on" << right_boundary_weight << "of"
                 << left.h_split;
        }
        assert(left_boundary_weight >= 0);
        assert(right_boundary_weight >= 0);
        return {l_begin, l_end, h_begin, h_end, size,
                num_light, num_heavy,
                left_boundary_weight, right_boundary_weight,
                has_left_boundary, has_right_boundary};
    }

    void fixup_boundary(const Subproblem &problem) {
        if (problem.has_right_boundary) {
            size_type id = indices_[problem.h_end - 1];
            double weight = problem.right_boundary_weight;
            sLOG << "Fixup right boundary: subtracting" << weight
                 << "from table entry" << id << table_[id];
            wrs::atomic_fetch_sub(&(table_[id].p), weight);
        }
    }


    wrs::numa_arr_ptr<tableitem> table_;
    wrs::numa_arr_ptr<size_type> indices_;
    wrs::numa_arr_ptr<double> prefsum_;
    wrs::numa_arr_ptr<size_type> classifier_;
    std::unique_ptr<Subproblem[]> subproblems_;
    std::unique_ptr<Splitter[]> splitters_;
    size_type size_;
    size_type num_light_;
    double W_, W_n_, eps_;
    std::vector<double> timers_;


};


} // namespace psa
} // namespace wrs

#endif // WRS_PSA_PSA_BASE_HEADER
