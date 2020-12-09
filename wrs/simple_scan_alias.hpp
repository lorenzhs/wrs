/*******************************************************************************
 * wrs/simple_scan_alias.hpp
 *
 * Simpler alias table construction without auxiliary arrays
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_SIMPLE_SCAN_ALIAS_HEADER
#define WRS_SIMPLE_SCAN_ALIAS_HEADER

#include <wrs/accumulate.hpp>
#include <wrs/memory.hpp>
#include <wrs/timer.hpp>
#include <wrs/util.hpp>
#include <wrs/verify.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace wrs {

template <typename size_type = uint32_t>
class simple_scan_alias {
public:
    static constexpr const char *name = "simplescan";
    static constexpr bool yields_single_sample = true;
    static constexpr bool init_with_seed = false;
    static constexpr int pass_rand = 1;
    using result_type = size_type;
    static constexpr size_type empty = static_cast<size_type>(-1);

    // Alias table member item
    struct tableitem {
        double p; // probability
        size_type a; // alias

        tableitem() : p(0), a(-1) {}
        tableitem(double p_, size_type a_) : p(p_), a(a_) {}
    };

    friend std::ostream &operator<<(std::ostream &os, const tableitem &item) {
        return os << '(' << item.p << ',' << item.a << ')';
    }

    static constexpr bool debug = false;
    static constexpr bool time = false;

    simple_scan_alias() : table_(nullptr), size_(-1), W_(0) {}

    template <typename Iterator>
    simple_scan_alias(Iterator w_begin, Iterator w_end) {
        init(w_end - w_begin);
        construct(w_begin, w_end);
    }

    simple_scan_alias &operator=(const simple_scan_alias &other) {
        LOG0 << "simple_scan_alias copy assignment/constructor called";
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
    simple_scan_alias(const simple_scan_alias &other) {
        *this = other;
    }
    //! delete move-constructor
    simple_scan_alias(simple_scan_alias &&) = delete;
    //! delete move-assignment
    simple_scan_alias &operator=(simple_scan_alias &&) = delete;

    void init(size_t size) {
        table_ = make_alloc_arr<tableitem>(size);
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

        Iterator i = begin, j = begin;
        while (i != end && !is_light(*i))
            i++; // light item
        while (j != end && is_light(*j))
            j++; // heavy item
        double w = (j != end) ? *j : 0; // residual weight of current heavy item
        if (j == end) {
            // Only light items, meaning every item has to be of weight W/n
            // (most likely, there is only one item), and this isn't even really
            // weighted sampling
            sLOG << "No heavy item found of a total of" << size_ << "items";
            while (i != end) {
                size_type iidx = i - begin;
                assert(std::abs(W_n_ - *i) < 1e-8);
                table_[iidx].p = *i;
                table_[iidx].a = iidx;
                i++;
            }
        }
        // Can't fill if only light or heavy items
        while (TLX_LIKELY(j != end)) {
            if (is_light(w)) {
                // pack a heavy bucket
                Iterator jprime = j + 1;
                size_type jidx = j - begin;
                // find next heavy item
                while (TLX_LIKELY(jprime != end) && is_light(*jprime))
                    jprime++;
                if (TLX_UNLIKELY(jprime == end)) {
                    // no more heavy items
                    sLOG << "residual heavy item" << jidx << "is light"
                         << "but search for next heavy item exceeded end"
                         << "- finishing up light items;"
                         << "i =" << i - begin << "j =" << jidx << "w =" << w;
                    assert(std::abs(W_n_ - w) < 1e-8);
                    table_[jidx].p = W_n_;
                    table_[jidx].a =
                        jidx; // w is W/n, this is just for numerical weirdness
                    for (Iterator k = i; k != end; ++k) {
                        assert(is_light(*k));
                        assert(std::abs(W_n_ - *k) < 1e-8);
                        size_type kidx = k - begin;
                        table_[kidx].p = *k;
                        table_[kidx].a = empty;
                    }
                    break;
                }
                table_[jidx].p = w;
                table_[jidx].a = jprime - begin;
                w = (w + *jprime) - W_n_;
                j = jprime;
            } else {
                // pack a light bucket with piece of heavy item
                if (TLX_UNLIKELY(i == end)) {
                    // the "light" piece of the heavy item that remains is W/n,
                    // so pack its bucket.
                    // TODO can we move this out of the loop?
                    assert(std::abs(W_n_ - w) < 1e-10 * W_);
                    size_type jidx = j - begin;
                    table_[jidx].p = W_n_;
                    table_[jidx].a = jidx;
                    break;
                }
                size_type iidx = i - begin;
                sLOG << "packing light bucket" << iidx << "with heavy item"
                     << j - begin << "residual weight" << w << "->"
                     << (w + *i) - W_n_;
                table_[iidx].p = *i;
                table_[iidx].a = j - begin;
                w = (w + *i) - W_n_;
                do {
                    i++;
                } while (i != end && !is_light(*i));
                if (i == end) {
                    sLOG << "residual heavy item" << j - begin
                         << "was still heavy but search for next light item"
                         << "after" << iidx << "exceeded end, i =" << i - begin
                         << "j =" << j - begin << "w =" << w;
                }
            }
        }

        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: table construction took " << timers_.back() << "ms";

        if constexpr (debug) {
            LOG_ARR(table_.get(), "table");
        }
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
        assert(item.a != empty || std::abs(item.p - W_n_) < 1e-7);

        // scale remaining randomness to range of bucket weights
        rand = (rand - candidates[0]) * W_n_;
        return candidates[rand >= item.p];
    }

    // Sum up all weights in the table, adding `offset` to every item ID.
    // This is used internally by verify() and by multi_alias' verify()
    void verify_helper(std::vector<double> &weights, size_t offset) const {
        wrs::sum_table<size_type>(table_.get(), size_, W_, W_n_, weights, offset);
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
    bool is_light(const double &d) const {
        return d <= W_n_;
    }

    alloc_arr_ptr<tableitem> table_;
    std::vector<double> timers_;
    size_t size_;
    double W_, W_n_;
};

} // namespace wrs

#endif // WRS_SIMPLE_SCAN_ALIAS_HEADER
