/*******************************************************************************
 * wrs/par_scan_alias2.hpp
 *
 * Parallel 1-alias table construction
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PAR_SCAN_ALIAS2_HEADER
#define WRS_PAR_SCAN_ALIAS2_HEADER

#include <wrs/psa/psa_base.hpp>
#include <wrs/timer.hpp>

#include <tlx/define.hpp>
#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace wrs {

template <typename size_type = int32_t>
class par_scan_alias2 : public psa::psa_base<par_scan_alias2, size_type> {
    using base = psa::psa_base<par_scan_alias2, size_type>;
    // base class types
    using tableitem = typename base::tableitem;
    using splitter = typename base::Splitter;
    using subproblem = typename base::Subproblem;
    // constexpr definitions
    using base::debug;
    using base::empty;
    using base::relative_epsilon;
    using base::time;

    // member variables
    using base::classifier_;
    using base::eps_;
    using base::indices_;
    using base::num_light_;
    using base::prefsum_;
    using base::size_;
    using base::splitters_;
    using base::subproblems_;
    using base::table_;
    using base::timers_;
    using base::W_;
    using base::W_n_;

    // (non-templated) functions
    using base::base_case;
    using base::fixup_boundary;
    using base::get_heavy;
    using base::get_item;
    using base::is_light;
    using base::set_heavy;
    using base::split;
    using base::subproblem_from_splitters;
    using base::verify_subproblem_feasibility;

public:
    using result_type = size_type;
    static constexpr const char* name = "psa2";

    par_scan_alias2() : base() {}

    template <typename Iterator>
    par_scan_alias2(Iterator begin, Iterator end) {
        init(end - begin);
        construct(begin, end);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        timer t;
        if (end - begin != static_cast<ssize_t>(size_)) {
            sLOG1 << "Error: tried to construct alias table of incorrect size!"
                  << "Expected" << size_ << "got" << end - begin;
            return;
        }
        timers_.clear();

        if (size_ <= 100)
            LOG << "input: " << std::vector<double>(begin, end);

        // Step 0: calculate total weight and weight per bucket
        W_ = parallel_accumulate(begin, end, 0.0);
        W_n_ = W_ / size_;

        eps_ = relative_epsilon * W_;
        LOG << "Using epsilon: " << eps_;

        const int num_threads = wrs::get_num_threads();
        sLOG << "have" << num_threads << "threads";

        // try to improve stability for normalized inputs
        if (std::abs(W_n_ - 1.0) < 1e-8) { // need fixed epsilon here
            sLOG << "Setting W/n to 1 - was off by" << W_n_ - 1.0 << "and W to"
                 << size_ << "- was off by" << W_ - static_cast<double>(size_);
            W_n_ = 1.0;
            W_ = static_cast<double>(size_);
        }
        LOG << "W = " << W_ << ", W/n = " << W_n_;
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 0: preprocessing took " << timers_.back() << "ms";

        // Step 1: best-effort process subranges of n/p items
        auto startpos = std::make_unique<size_type[]>(num_threads);
        auto subsizes = std::make_unique<size_type[]>(num_threads + 1);
        auto sequential_assignment = [&](size_type min, size_type max, int thread) {
            auto my_end = begin + max;
            auto [last, unfinished] = best_effort_fill(begin, begin + min, my_end);
            startpos[thread] = last - begin;
            subsizes[thread] = unfinished;
            sLOG << "Thread" << thread << "filled" << min << "to"
                 << startpos[thread] << "->" << unfinished << "remaining";
        };
        parallel_do_range(sequential_assignment, size_);
        LOG_ARR(table_.get(), "initial table");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: best-effort fill took " << timers_.back() << "ms";

        wrs::exclusive_scan(subsizes.get(), subsizes.get() + num_threads,
                            subsizes.get(), std::plus<>(), 0, true);
        LOG_ARRs(startpos.get(), "startpos", num_threads);
        LOG_ARRs(subsizes.get(), "subsizes", num_threads + 1);
        remaining_ = subsizes[num_threads];
        sLOG1 << remaining_ << "of" << size_ << "remaining";

        // We're done here :)
        if (remaining_ <= 1) {
            sLOG << "Base case handled everything, we're done here";
            return;
        }

        // Step 2: classify remainder
        parallel_do_range(
            [&](size_type min, size_type max, int thread) {
                (void)min; // prevent unused parameter warning
                assert(min <= startpos[thread]);
                assert(max >= startpos[thread]);
                size_type pos = startpos[thread];
                sLOG << "Thread" << thread << "classifying range"
                     << startpos[thread] << ".." << max << "output to pos"
                     << subsizes[thread] << "...";
                for (size_type index = subsizes[thread]; pos < max; pos++) {
                    // skip finished buckets
                    if (!is_unfinished(table_[pos])) {
                        sLOG << "classifier skipping finished bucket" << pos
                             << table_[pos];
                        continue;
                    }
                    LOG << "classifier_[" << index << "] = " << table_[pos].p
                        << " < W/n (item " << pos << ")";
                    classifier_[index++] = !is_light(table_[pos]);
                }
            },
            size_);
        LOG_ARRs(classifier_.get(), "classifier", remaining_);
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 2: classification & table init took "
                   << timers_.back() << "ms";

        // Step 3: prefix sum
        parallel_scan(classifier_.get(), classifier_.get() + remaining_,
                      classifier_.get(), std::plus<>());
        num_light_ = remaining_ - classifier_[remaining_ - 1];
        LOG_ARRs(classifier_.get(), "post scan", remaining_);
        timers_.push_back(t.get_and_reset());
        LOG << "num_light_ = " << num_light_ << " remaining_ = " << remaining_
            << " size_ = " << size_;
        LOGC(time) << "Step 3: prefsum took " << timers_.back() << "ms";

        if (num_light_ == 0) {
            // All items are considered heavy due to numerical issues, but
            // actually, we're done here
            sLOG
                << "num_light_ == 0, aborting -- just some numerical weirdness";
            return;
        }

        // Step 4: write indices_ array. But first handle the first item in both
        // possible outcomes (the wrong one will be overwritten)
        indices_[0] = 0;
        indices_[num_light_] = 0;
        // XXX
        // prefsum_[0] = *begin; prefsum_[num_light_] = *begin;

        const size_type light_offset = num_light_ - 1;
        parallel_do_range(
            [&](size_type /* ignored */, size_type max, int thread) {
                size_type min = std::max<size_type>(0, startpos[thread]);
                size_type item_index;
                size_type out_index = subsizes[thread];
                size_type prev = out_index == 0 ? 0 : classifier_[out_index - 1];
                sLOG << "Thread" << thread << "computing indices for item range"
                     << min << ".." << max << "output to indices_[] pos"
                     << out_index << ".." << subsizes[thread + 1]
                     << "prev =" << prev << "light_offset =" << light_offset;

                for (item_index = min; item_index < max; out_index++, item_index++) {
                    // Skip finished items
                    while (item_index < max && !is_unfinished(table_[item_index])) {
                        ++item_index;
                    }
                    if (item_index >= max)
                        break;
                    size_type curr = classifier_[out_index];
                    size_type pos;
                    if (curr == prev) {
                        assert(out_index >= curr && out_index - curr < num_light_);
                        pos = out_index - curr;
                        sLOG << "Thread" << thread << "pos" << item_index
                             << curr << "= curr = prev -> pos =" << pos;
                    } else {
                        assert(light_offset + curr < size_);
                        pos = light_offset + curr;
                        sLOG << "Thread" << thread << "pos" << item_index << curr
                             << "= curr != prev =" << prev << "-> pos =" << pos;
                    }
                    prev = curr;
                    assert(0 <= pos && pos < size_);
                    indices_[pos] = item_index;
                    prefsum_[pos] = table_[item_index].p;
                }
            },
            size_);
        LOG_ARRs(indices_.get(), "indices", remaining_);
        if constexpr (debug) {
            std::vector<tableitem> t;
            double sum = 0.0;
            for (size_type i = 0; i < remaining_; i++) {
                t.push_back(table_[indices_[i]]);
                sum += t.back().p;
            }
            LOG_ARRs(t.data(), "table-subset", remaining_);
            sLOG << "weight sum =" << sum << "=" << sum / W_n_ << "* W/n";
        }
        LOG_ARRs(prefsum_.get(), "L/H weights", remaining_);
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 4: computing indices took " << timers_.back() << "ms";

        // Step 5: prefix sum over prefsum_ weights
        parallel_scan(prefsum_.get(), prefsum_.get() + num_light_,
                      prefsum_.get(), std::plus<>());
        parallel_scan(prefsum_.get() + num_light_, prefsum_.get() + remaining_,
                      prefsum_.get() + num_light_, std::plus<>());
        LOG_ARR(prefsum_.get(), "L/H prefsum");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 5: weight prefix sum took " << timers_.back() << "ms";

        // Step 6: split up
        // clang-format off
        subproblem problem = {0, num_light_, num_light_, remaining_, remaining_,
                              num_light_, remaining_ - num_light_, 0.0, 0.0,
                              false, false};
        // clang-format on

        // 6a: compute splitters
        splitters_[0] = {0, num_light_, 0.0, 0, false, false};
        splitters_[num_threads] = {num_light_, remaining_, 0,
                                   0.0,        false,      false};
        wrs::parallel_do(
            [&](size_type i) {
                double fraction = static_cast<double>(i) / num_threads;
                splitters_[i] = split(problem, fraction);
                sLOG << "splitter" << i << splitters_[i];
            },
            num_threads, 1);
        LOG << "splitting done.\n";

        // 6b: compute subproblems from neighbouring splitters and fix boundaries
        wrs::parallel_do(
            [&](size_type i) {
                subproblems_[i] =
                    subproblem_from_splitters(splitters_[i], splitters_[i + 1]);
                if constexpr (debug)
                    verify_subproblem_feasibility(subproblems_[i]);
                fixup_boundary(subproblems_[i]);
                sLOG << "problem" << i << subproblems_[i];
            },
            num_threads);
        LOG << "splitter to problem conversion done.\n";
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 6: splitting took " << timers_.back() << "ms";

        // Step 7: call base cases
        wrs::parallel_do(
            [&](size_type i) {
                sLOG << "\nCalling base case" << i << subproblems_[i];
                base_case(subproblems_[i]);
            },
            num_threads);

        LOG_ARR(table_.get(), "table");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 7: base cases took " << timers_.back() << "ms";
    }

protected:
    TLX_ATTRIBUTE_ALWAYS_INLINE
    bool is_unfinished(tableitem& i) const {
        return (i.p == 0 || i.p > W_n_ || (i.p < W_n_ && i.a == empty));
    }

    // returns iterator to first unfinished item + number of unfinished
    template <typename Iterator>
    std::pair<Iterator, size_type> best_effort_fill(const Iterator base,
                                                    Iterator begin, Iterator end) {
        LOG << "best_effort_fill(" << begin - base << ", " << end - base << ")";
        Iterator i = begin, j = begin;
        while (i != end && !is_light(*i))
            i++; // light item
        while (j != end && is_light(*j))
            j++; // heavy item
        double w = (j != end) ? *j : 0; // residual weight of current heavy item
        // Can't fill if only light or heavy items
        while (i != end && j != end) {
            if (is_light(w)) {
                // pack a heavy bucket
                Iterator jprime = j;
                size_type jidx = j - base;
                // find next heavy item
                do {
                    jprime++;
                } while (jprime != end && is_light(*jprime));
                if (jprime == end) {
                    // no more heavy items
                    sLOG << "residual heavy item" << jidx << "is light"
                         << "but search for next heavy item exceeded end"
                         << "aborting; i =" << i - base << "j =" << jidx
                         << "w =" << w;
                    break;
                }
                table_[jidx].p = w;
                table_[jidx].a = jprime - base;
                w = (w + *jprime) - W_n_;
                j = jprime;
            } else {
                // pack a light bucket with piece of heavy item
                size_type iidx = i - base;
                sLOG << "packing light bucket" << iidx << "with heavy item"
                     << j - base << "residual weight" << w << "->"
                     << (w + *i) - W_n_;
                table_[iidx].p = *i;
                table_[iidx].a = j - base;
                w = (w + *i) - W_n_;
                do {
                    i++;
                } while (i != end && !is_light(*i));
                if (i == end) {
                    sLOG << "residual heavy item" << j - base
                         << "was still heavy but search for"
                         << "next light item exceeded end, aborting; i =" << i - base
                         << "j =" << j - base << "w =" << w;
                    break;
                }
            }
        }
        assert(begin <= i && i <= end);
        assert(begin <= j && j <= end);
        // OK, we've processed as many items as we can right now. Now we need to
        // do two things:
        // 1) Count how many items are unfinished
        // 2) Initialise their table entries
        size_type unfinished = 0;
        Iterator last = std::min(i, j);
        // Only one of the next two loops will execute anything
        for (Iterator it = i; it < j; ++it) {
            // skip heavy items between i and j (they're finished)
            if (!is_light(*it))
                continue;
            sLOG0 << "setting light item" << it - base << *it << table_[it - base];
            unfinished++;
            table_[it - base].p = *it;
            table_[it - base].a = empty;
        }
        for (Iterator it = j; it < i; ++it) {
            // skip light items between j and i (they're finished)
            if (is_light(*it))
                continue;
            sLOG0 << "setting heavy item" << it - base << *it << table_[it - base];
            unfinished++;
            table_[it - base].p = *it;
            table_[it - base].a = empty;
        }
        // Items after i and j are all unprocessed, no need to explicitly count,
        // but we still need to fill the table entries
        Iterator it = std::max(i, j);
        unfinished += end - it;
        for (Iterator it = std::max(i, j); it < end; ++it) {
            sLOG0 << "last step processing" << it - base << *it;
            table_[it - base].p = *it;
            table_[it - base].a = empty;
        }
        // The last heavy item already had part of its weight assigned, restore
        // that entry
        if (j != end) {
            sLOG << "resetting item" << j - base << "to w =" << w
                 << table_[j - base];
            table_[j - base].p = w;
        }
        return std::make_pair(last, unfinished);
    }


    size_type remaining_;
};


} // namespace wrs

#endif // WRS_PAR_SCAN_ALIAS2_HEADER
