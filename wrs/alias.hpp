/*******************************************************************************
 * wrs/alias.hpp
 *
 * Alias method implementations
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_ALIAS_HEADER
#define WRS_ALIAS_HEADER

#include <wrs/accumulate.hpp>
#include <wrs/memory.hpp>
#include <wrs/timer.hpp>
#include <wrs/util.hpp>
#include <wrs/verify.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace wrs {

template <typename size_type = uint32_t>
struct alias {
    static constexpr const char *name = "alias";
    static constexpr bool yields_single_sample = true;
    static constexpr bool init_with_seed = false;
    static constexpr int pass_rand = 1;
    using result_type = size_type;
    using pair = std::pair<size_type, double>;

    static_assert((sizeof(double) / sizeof(size_type)) * sizeof(size_type) ==
                      sizeof(double),
                  "size_type has weird size");
    static constexpr size_t ti_ptroffset = sizeof(double) / sizeof(size_type);

    // Alias table member item, used to be a std::tuple<double, size_type, size_type>
    struct tableitem {
        double p; // probability
        size_type a; // alia

        tableitem() : p(0), a(-1) {}
        tableitem(double p_, size_type a_) : p(p_), a(a_) {}
    };

    friend std::ostream &operator<<(std::ostream &os, const tableitem &item) {
        return os << '(' << item.p << ',' << item.a << ')';
    }

    static constexpr bool debug = false;
    static constexpr bool time = false;

    alias() : table_(nullptr), work_(nullptr), size_(-1), W_(0) {}

    template <typename Iterator>
    alias(Iterator w_begin, Iterator w_end) {
        init(w_end - w_begin);
        construct(w_begin, w_end);
    }

    alias &operator=(const alias &other) {
        LOG0 << "alias copy assignment/constructor called";
        size_ = other.size_;
        W_ = other.W_;
        W_n_ = other.W_n_;
        if (size_ != static_cast<size_t>(-1)) {
            table_ = make_alloc_arr<tableitem>(size_);
            memcpy(table_.get(), other.table_.get(), size_ * sizeof(tableitem));
        } else {
            table_ = nullptr;
        }
        timers_ = other.timers_;
        return *this;
    }
    alias(const alias &other) {
        *this = other;
    }
    //! delete move-constructor
    alias(alias &&) = delete;
    //! delete move-assignment
    alias &operator=(alias &&) = delete;

    void init(size_t size) {
        table_ = make_alloc_arr<tableitem>(size);
        work_ = make_alloc_arr<pair>(size);
        size_ = size;
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

        W_ = wrs::accumulate(begin, end, 0.0);
        W_n_ = W_ / size_;
        LOG << "W = " << W_ << ", W/n = " << W_n_;
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 0: preprocessing took " << timers_.back() << "ms";

        // Classify into small and large items
        ssize_t i_small = 0, i_large = size_ - 1;
        for (Iterator it = begin; it != end; ++it) {
            if (is_small(*it)) {
                LOG << *it << " is small";
                work_[i_small++] = {it - begin, *it};
            } else {
                LOG << *it << " is large";
                work_[i_large--] = {it - begin, *it};
            }
        }
        assert(i_small > i_large); // must meet
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: classification took " << timers_.back() << "ms";

        // Assign items
        auto s_begin = work_.get(), s_end = work_.get() + i_small, l_begin = s_end,
             l_end = work_.get() + size_, small = s_begin, large = l_begin;
        LOG_ARR(work_.get(), "work");
        while (small != s_end && large != l_end) {
            LOG << "table_[" << small->first << "] = (" << small->second << ", "
                << small->first << ", " << large->first << ")";
            table_[small->first] = tableitem(small->second, large->first);
            large->second = (large->second + small->second) - W_n_;

            if (is_small(large)) {
                sLOG << "large item became small, remaining weight:" << *large
                     << "moving to end of small items";
                if (large > s_end) {
                    std::swap(large, s_end);
                } else {
                    LOG << "no need to swap, is first";
                }
                s_end++;
                // l_begin++; // not really needed, just for bookkeping
                large++;
            } else {
                sLOG << "large item" << *large << "remained large, advancing small";
            }
            small++;
        }
        while (large != l_end) {
            size_t i_large = large->first;
            LOG << "large item left at the end: " << *large;
            LOG << "table_[" << i_large << "] = (" << W_n_ << ", " << i_large
                << ", -1)";

            table_[i_large] = tableitem(W_n_, i_large);
            large++;
        }
        while (small != s_end) {
            size_t i_small = small->first;
            sLOG << "encountered some numerical instability!"
                 << "Item" << i_small << "weight" << *small;
            LOG << "table_[" << i_small << "] = (" << W_n_ << ", " << i_small
                << ", -1)";

            table_[i_small] = tableitem(W_n_, i_small);
            small++;
        }
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 2: assignment took " << timers_.back() << "ms";
    }

    // Query given a uniformly distributed [0,1) random value
    size_type sample(double rand) const {
        // scale up to table size
        rand *= size_;

        // Copy the item's index and its alias into temp_
        size_type candidates[2];
        candidates[0] = static_cast<size_type>(rand); // bucket ID (index)
        candidates[1] = table_[candidates[0]].a; // fetch bucket alias

        // reference to the item, for convenience
        tableitem &item = table_[candidates[0]];
        // sanity check: if the item doesn't have an alias, then it must fill
        // the bucket on its own
        assert(item.a != static_cast<size_type>(-1) ||
               std::abs(item.p - W_n_) < 1e-7);

        // scale remaining randomness to range of bucket weights
        rand = (rand - candidates[0]) * W_n_;
        return candidates[rand >= item.p];
    }

    // Sum up all weights in the table, adding `offset` to every item ID.
    // This is used internally by verify() and by multi_alias' verify()
    void verify_helper(std::vector<double> &weights, size_t offset) const {
        wrs::sum_table<size_type>(table_.get(), size_, W_, W_n_, weights, offset);
    }

    template <typename Iterator>
    void verify(Iterator begin, Iterator end) const {
        (void)begin;
        (void)end;
#ifndef NDEBUG
        assert(end - begin == static_cast<ssize_t>(size_));
        std::vector<double> weights(size_);
        verify_helper(weights, 0);
        if (size_ < 100) {
            LOG_ARR(table_.get(), "table");
            LOG << "reconstructed: " << weights;
        }
        wrs::check_weights(begin, end, weights, *this);
#endif
    }

    std::vector<double> get_timers() const {
        return timers_;
    }

    double total_weight() const {
        return W_;
    }

    size_t size() const {
        return size_;
    }

    template <typename Callback>
    void find(size_type item, Callback &&callback) const {
        if (size_ == static_cast<size_t>(-1))
            return;
        if (item < size_)
            callback(0, item, table_[item].p, table_[item]);
        for (size_t i = 0; i < size_; i++) {
            if (table_[i].a == item) {
                callback(0, i, W_n_ - table_[i].p, table_[i]);
            }
        }
    }

protected:
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool is_small(const double &d) const {
        return d <= W_n_;
    }

    template <typename Iterator>
    TLX_ATTRIBUTE_ALWAYS_INLINE bool is_small(Iterator it) const {
        return it->second <= W_n_;
    }

    alloc_arr_ptr<tableitem> table_;
    alloc_arr_ptr<pair> work_;
    std::vector<double> timers_;
    size_t size_;
    double W_, W_n_;
};

} // namespace wrs

#endif // WRS_ALIAS_HEADER
