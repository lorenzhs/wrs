/*******************************************************************************
 * wrs/par_alias.hpp
 *
 * Parallel alias method
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PAR_ALIAS_HEADER
#define WRS_PAR_ALIAS_HEADER

#include <wrs/numa_array.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/timer.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>
#include <tlx/logger/all.hpp>

#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <utility>

// Log messages that should be sLOG1 if assertions are enabled
#ifdef NDEBUG
#define sDLOG sLOG
#else
#define sDLOG sLOG1
#endif

namespace wrs {

template <typename size_type = uint32_t>
struct par_alias {
    static constexpr const char* name = "parallel";
    static constexpr bool pass_rand = true;
    using result_type = size_type;
    // (prob_1, prob_2, item, alias_1, alias_2, iteration)
    // we get the space for the iteration for free because of alignment rules
    using tableitem = std::tuple<double, double, size_type, size_type, size_type, uint32_t>;
    using pair = std::pair<size_type, double>;

    static_assert((sizeof(double) / sizeof(size_type)) * sizeof(size_type) == sizeof(double),
                  "size_type has weird size");
    static constexpr size_t ti_ptroffset = 2 * sizeof(double) / sizeof(size_type);

    static constexpr bool debug = false;
    static constexpr bool time = false;

    // relative epsilon, will be multiplied with total weight sum
    static constexpr double rel_eps = 1e-10;

    par_alias() {}

    template <typename Iterator>
    par_alias(Iterator w_begin, Iterator w_end) {
        init(w_end - w_begin);
        construct(w_begin, w_end);
    }

    par_alias & operator = (const par_alias &other) {
        LOG0 << "par_alias copy assignment/constructor called";
        size_ = other.size_;
        // Allocate table locally
        table_ = make_numa_arr<tableitem>(
            size_, /* align */ true, /* local */ true);
        memcpy(static_cast<void*>(table_.get()), other.table_.get(),
               size_ * sizeof(tableitem));
        timers_ = other.timers_;
        return *this;
    }
    par_alias(const par_alias &other) {
        *this = other;
    }
    //! delete move-constructor
    par_alias(par_alias &&) = delete;
    //! delete move-assignment operator
    par_alias & operator = (par_alias &&) = delete;

    void init(size_t size) {
        table_ = make_numa_arr<tableitem>(size);
        classifier = make_numa_arr<size_type>(size);
        F = make_numa_arr<double>(size);
        L = make_numa_arr<pair>(size);
        it_ = 0;
        size_ = size;

        std::fill(table_.get(), table_.get() + size,
                  std::make_tuple(0, 0, -1, -1, -1, 0));
    }

    template <typename Iterator>
    void construct(Iterator w_begin, Iterator w_end) {
        timer t;
        if (w_end - w_begin != static_cast<ssize_t>(size_)) {
            sLOG1 << "Error: tried to construct alias table of incorrect size!"
                  << "Expected" << size_ << "got" << w_end - w_begin;
            return;
        }
        if (size_ > static_cast<size_t>(std::numeric_limits<size_type>::max())) {
            sLOG1 << "Error: size_type cannot hold" << size_
                  << "- please specify an appropriate type (e.g. size_t)";
            return;
        }
        timers_.clear();

        // Start a new iteration
        ++it_;

        // Step 0: calculate total weight and weight per bucket
        W_ = parallel_accumulate(w_begin, w_end, 0.0);
        W_n_ = W_ / size_;

        eps = rel_eps * W_;
        LOG << "Using epsilon: " << eps;

        // try to improve stability for normalized inputs
        if (std::abs(W_n_ - 1.0) < 1e-8) { // need fixed epsilon here
            LOG << "Setting W/n to 1 (was off by " << W_n_ - 1.0 << ")";
            W_n_ = 1.0;
        }
        LOG << "W = " << W_ << ", W/n = " << W_n_;
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 0: preprocessing took " << timers_.back() << "ms";

        // Step 1: classify
        parallel_do([&](size_t index) {
                classifier[index] = *(w_begin + index) <= W_n_;
            }, size_);
        log_arr(classifier.get(), "classifier");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: classification took " << timers_.back() << "ms";

        // Step 2: prefix sum
        parallel_scan(classifier.get(), classifier.get() + size_,
                      classifier.get(), std::plus<>());
        // This implies the number of large items
        L_size = size_ - classifier[size_ - 1];
        log_arr(classifier.get(), "post scan");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 2: prefsum took " << timers_.back() << "ms";

        // Step 3: assign small items and extract F, L at the same time
        parallel_do_range(
            [&](size_t min, size_t max, int) {
                size_type last = (min > 0) ? classifier[min - 1] : 0;
                for (size_t i = min; i < max; i++) {
                    auto pos = classifier[i];
                    assert(pos >= last);
                    if (pos != last) {
                        const double &weight = *(w_begin + i);
                        // offset -1 because the prefix sum is obviously 1-based
                        table_[pos - 1] = { weight, weight, i, -1, -1, it_ };
                        F[pos - 1] = W_n_ - weight;
                        sLOG << "Processing small item" << pos << "aka item" << i
                             << "of weight" << weight
                             << "writing to bucket" << pos-1;
                        last = pos;
                    } else {
                        // it's a large item
                        sLOG << "Processing large item" << i - pos << "aka item" << i
                             << "of weight" << *(w_begin + i);
                        L[i - pos] = std::make_pair(i, *(w_begin + i));
                    }
                }
            }, size_);

        log_arr(L.get(), "L", 50, L_size);
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 3: small item assignment took " << timers_.back() << "ms";

/******************************************************************************/

        // Step 4: prefix sum over remaining space (F) and large items' weights (L)
        auto pair_plus = [](const pair &a, const pair &b) {
            auto res = std::make_pair(b.first, a.second + b.second);
            return res;
        };
        // 4a: assign F array for and initialize empty buckets
        parallel_do_range(
            [&](size_t min, size_t max, int) {
                for (size_t i = min; i < max; i++) {
                    F[i] = W_n_;
                }
            }, size_, classifier[size_ - 1]);
        log_arr(F.get(), "F");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 4a: F & T filling took " << timers_.back() << "ms";

        // 4b: prefix sum over F
        parallel_scan(F.get(), F.get() + size_, F.get(), std::plus<>());
        log_arr(F.get(), "F (psum)");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 4b: F prefsum took " << timers_.back() << "ms";

        // 4d: prefix sum over L
        parallel_scan(L.get(), L.get() + L_size, L.get(), pair_plus);
        log_arr(L.get(), "L (psum)", 50, L_size);
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 4d: L prefsum took " << timers_.back() << "ms";

        // TODO increase numerical stability of prefix sum computation
        double diff = F[size_ - 1] - L[L_size - 1].second;
        if (std::abs(diff) >= eps) {
            sDLOG << "F and L don't add up to the same thing. F:"
                  << F[size_ - 1] << "L:" << L[L_size - 1].second
                  << "diff =" << diff;
            assert(false);
        }

/******************************************************************************/

        // Step 5: assign large items
        auto assign_large = [&](size_t min, size_t max, int thread_id) {
            sLOG << "Thread" << thread_id << "processing F range" << min << max;

            const auto F_begin = F.get() + min, F_end = F.get() + max;
            double needle_begin = *F_begin;
            auto needle_cmp = [](const pair &p, const double &needle) {
                return p.second < needle;
            };
            auto L_begin = std::lower_bound(
                L.get(), L.get() + L_size, needle_begin, needle_cmp),
            L_end = L.get() + L_size; // for last PE
            if (max < size_) { // if not last PE, find needle
                double needle_end = *(F_end - 1);
                // The +1 means that the next thread's L_begin is likely our
                // L_end - 1. That's okay, this large item needs to be
                // handled by both because bucket ownership is exclusive.
                L_end = std::lower_bound(
                    L_begin, L.get() + L_size, needle_end, needle_cmp) + 1;

                sLOG << "Thread" << thread_id << "found L range"
                     << L_begin - L.get() << "to" << L_end - L.get() << "of"
                     << L_size << "needles" << needle_begin << needle_end;
            } else {
                sLOG << "Thread" << thread_id << "found L range"
                << L_begin -L.get() << "to" << L_size
                << "(last thread) for needle" << needle_begin;
            }
            sLOG << "Thread" << thread_id
            << "F from" << *F_begin << "to" << *(F_end - 1)
            << "L from" << *L_begin << "to" << *(L_end - 1);
            assert(*(F_end - 1) <= (L_end - 1)->second + eps);

            // Check how much of *L_begin is handled by the preceding thread
            double needed = *(w_begin + L_begin->first),
            taken_care_of = 0.0;
            if (L_begin > L.get() &&
                L_begin->second > *(F_begin - 1) + eps)
            {
                auto prev_F = F_begin - 1;
                // Check if our L_begin was overly optimistic
                if ((L_begin - 1)->second > *prev_F) {
                    sLOG << "Thread" << thread_id << "takes a step back";
                    --L_begin;
                }

                sLOG << "Thread" << thread_id << "found a bit of an overlap"
                     << "considering new L_begin =" << *L_begin
                     << "and preceding PE's last bucket," << *prev_F;

                needed = L_begin->second - *prev_F;
                double weight = *(w_begin + L_begin->first);
                taken_care_of = weight - needed;

                assert(needed > 0);
                sLOG << "Thread" << thread_id
                     << "prev L weight:" << L_begin->second
                     << "prev F cap:" << *prev_F
                     << "still needed:" << needed
                     << "taken care of:" << taken_care_of;
            }

            auto F_it = F_begin;
            auto L_it = L_begin;
            while (F_it != F_end && L_it != L_end) {
                size_t F_idx = F_it - F.get();
                tableitem &item = table_[F_idx];
                assert(needed > 0);
                assert(taken_care_of >= 0);
                assert(std::get<0>(item) <= std::get<1>(item) + eps);
                sLOG << "Thread" << thread_id << "item" << L_it - L.get()
                     << "/" << L_it->first << "needs" << needed
                     << "space, using bucket" << F_idx << "which has"
                     << (std::get<5>(item) == it_ ?
                         W_n_ - std::get<1>(item) : W_n_)
                     << "left";

                auto [filled, remaining] =
                    fill_item(F_idx, L_it->first, needed);
                taken_care_of += filled;
                needed -= filled;

                if (remaining > eps || needed < eps) {
                    // there is remaining capacity in F
                    sLOG << "Thread" << thread_id << "bucket"
                         << F_idx << "has" << remaining
                         << "cap remaining after exceeding item"
                         << L_it - L.get();
                    // Everything should be taken care of, nothing should be
                    // needed any more
                    assert(std::abs(taken_care_of - *(w_begin + L_it->first)) < eps);
                    assert(needed < eps);
                    ++L_it;
                    if (L_it - L.get() != static_cast<ssize_t>(L_size)) {
                        // don't access beyond the end if this was the end
                        needed = *(w_begin + L_it->first);
                        taken_care_of = 0.0;
                    }

                    if (needed < eps) {
                        sLOG << "Thread" << thread_id << "bucket" << F_idx
                             << "is exhausted at the same time! Advancing.";
                        ++F_it;
                    }
                } else {
                    sLOG << "Thread" << thread_id << "item"
                         << L_it - L.get() << "filled bucket"
                         << F_it - F.get() << "proceeding with next bucket";
                    ++F_it;
                }
            }

            if constexpr (debug) {
                if (max < size_) { // we're not the last PE
                    double remaining = W_n_ - std::max(
                        std::get<0>(table_[F_end - F.get() - 1]),
                        std::get<1>(table_[F_end - F.get() - 1]));
                    sLOG << "Thread" << thread_id << "last bucket has"
                    << remaining << "free capacity, last item"
                    << L_it - L.get() << *L_it << "needs another"
                    << needed << "in next thread's allocation, already took care of"
                    << taken_care_of << "L at end:" << (L_it == L_end)
                    << "F at end:" << (F_it == F_end);
                }
            }
        };
        parallel_do_range(assign_large, size_);
        log_arr(table_.get(), "table");
        sLOG << "Table end:" << table_[size_ - 2] << table_[size_ - 1];
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 5a: large item assignment took " << timers_.back() << "ms";
    }

    size_type sample(double uniform) const {
        double rand = uniform * size_;
        size_t index = rand;
        tableitem& item = table_[index];

        // sanity check
        if constexpr(debug) {
            auto [p1, p2, i1, i2, i3, it] = item;
            assert(it == it_);
            if (i2 == static_cast<size_type>(-1)) {
                assert(std::abs(p1 - W_n_) < eps);
            }
            if (i3 == static_cast<size_type>(-1)) {
                assert(p1 < p2 + eps && std::abs(p2 - W_n_) < eps);
            }
        }

        rand = (rand - index) * W_n_;
        double &p1 = std::get<0>(item), &p2 = std::get<1>(item);
        size_t offset = (rand >= p1) + (rand >= p2);

        // more sanity checks
        assert(rand < p1 || std::get<3>(item) != static_cast<size_type>(-1));
        assert(rand < p2 || std::get<4>(item) != static_cast<size_type>(-1));

        // Ugly hack, but faster than a conditional
        return *(reinterpret_cast<size_type*>(&item) + ti_ptroffset + offset);
    }


    template <typename Iterator>
    void verify(Iterator w_begin, Iterator w_end) {
        sanity_check();
        (void) w_begin;
        (void) w_end;
#ifndef NDEBUG
        assert(w_end - w_begin == static_cast<ssize_t>(size_));
        std::vector<double> weights(size_);
        for (size_t i = 0; i < size_; i++) {
            auto [p1, p2, i1, i2, i3, it] = table_[i];
            assert(p1 > 0);
            assert(i1 != static_cast<size_type>(-1));
            weights[i1] += p1;
            if (p1 < p2) {
                assert(i2 != static_cast<size_type>(-1));
                weights[i2] += p2 - p1;
            }
            if (p2 + eps < W_n_) {
                if (i3 == static_cast<size_type>(-1)) {
                    sLOG1 << "verify: numerical instability in bucket" << i
                          << "have W/n - p2 = " << W_n_ - p2 << "- skipping";
                    continue;
                }
                weights[i3] += W_n_ - p2;
            }
        }
        for (size_t i = 0; i < size_; ++i) {
            double should = *(w_begin + i);
            double have = weights[i];
            assert(std::abs(should - have) < eps);
        }
        LOG1 << "Verification succeeded!";
#endif
    }

    size_t size() const {
        return size_;
    }

    std::vector<double> get_timers() const {
        return timers_; // copy
    }

    double total_weight() const {
        return W_;
    }

protected:
    // returns std::make_pair(how much was filled, how much still remains)
    std::pair<double, double> fill_item(size_t tableidx, size_t item_id,
                                        double needed) {
        tableitem &item = table_[tableidx];
        auto [p1, p2, i1, i2, i3, it] = item;
        const bool is_current_it = (it == it_);

        assert(!is_current_it || p1 < p2 + eps);
        // if alias1 is -1, p2 must be p1
        assert(!is_current_it || i2 != static_cast<size_type>(-1) || p2 < p1 + eps);
        // if p2 = 0, alias2 must be -1 or the item unitialized (all 0)
        assert(!is_current_it || p2 > eps || i3 == static_cast<size_type>(-1) ||
               (p1 < eps && p2 < eps && i1 == 0 && i2 == 0 && i3 == 0));
        // If the item is still dirty from the previous iteration, it's empty
        double remaining = is_current_it ? W_n_ - p2 : W_n_;

        if constexpr (debug) {
            if (is_current_it &&
                (i1 == item_id || i2 == item_id || i3 == item_id)) {
                sLOG1 << "fill_item: aborting, item" << tableidx
                      << "already has some item" << item_id << item;
                if (i1 == item_id) return std::make_pair(p1, remaining);
                if (i2 == item_id) return std::make_pair(p2, remaining);
                else return std::make_pair(remaining, 0.0);
            }
        } else {
            assert(!is_current_it ||
                   (i1 != item_id && i2 != item_id && i3 != item_id));
        }

        sLOG << "F @" << tableidx << "has" << remaining << "remaining, need"
             << needed << "-- filling with" << item_id;

        double to_fill = std::min(remaining, needed);
        if (!is_current_it) {
            // The item is dirty from the previous iteration, so it's empty and
            // we can fill the primary slot
            sLOG << "F @" << tableidx << "using primary slot - empty item" << to_fill;
            item = { to_fill, to_fill, item_id, -1, -1, it_ };
        } else if (p2 > p1 + eps) {
            sLOG << "F @" << tableidx << "using alias2 slot";
            // fill the remaining space as alias2
            assert(remaining <= needed);
            // item = std::make_tuple(p1, p2, i1, i2, item_id);
            std::get<4>(item) = item_id;
            return std::make_pair(remaining, 0.0);
        } else {
            assert(p1 > 0);
            sLOG << "F @" << tableidx << "using alias1 slot" << to_fill;
            // we're alias1, fill as much as we need
            // item = std::make_tuple(p1, p1 + to_fill, i1, item_id, -1, it_);
            std::get<1>(item) = std::get<0>(item) + to_fill;
            std::get<3>(item) = item_id;
        }
        return std::make_pair(to_fill, remaining - to_fill);
    }

    void sanity_check() const {
        for (size_t i = 0; i < size_; i++) {
            auto [p1, p2, i1, i2, i3, it] = table_[i];
            assert(it == it_);
            assert(p1 <= p2 && p2 <= W_n_ + eps);
            if (i == size_ - 1 && i2 == static_cast<size_type>(-1)) {
                // We can also have some instability at the boundary between
                // threads, but there's not much we can do about that
                if (std::abs(p2 - W_n_) > eps) {
                    LOG << "Numerical instability adds up to "
                        << (W_n_ - p2) * 100 / W_n_
                        << "% of last bucket wasted";
                }
            } else if (i1 == static_cast<size_type>(-1)) {
                assert(p1 < eps && p2 < eps);
            } else if (i2 == static_cast<size_type>(-1)) {
                assert(std::abs(p2 - W_n_) < eps && std::abs(p1 - W_n_) < eps);
            } else if (i3 == static_cast<size_type>(-1)) {
                assert(p1 < p2 + eps && std::abs(p2 - W_n_) < eps);
            }
        }
        LOG << "Sanity check passed";
    }


    template <typename Array>
    void log_arr(const Array &arr, const std::string &desc,
                 size_t max_size = 50, size_t size = 0)
    {
        using value_type = std::remove_reference_t<decltype(arr[0])>;
        if (size == 0) size = size_;
        LOG << desc << ": " << std::vector<value_type>(
            arr, arr + std::min(max_size, size));
    }

    numa_arr_ptr<tableitem> table_;
    numa_arr_ptr<size_type> classifier;
    numa_arr_ptr<double> F;
    numa_arr_ptr<pair> L;
    std::vector<double> timers_; // for measurements
    size_t size_, L_size;
    double W_, W_n_, eps;
    uint32_t it_;
};

} // namespace wrs

#endif // WRS_PAR_ALIAS_HEADER
