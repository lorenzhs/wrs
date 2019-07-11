/*******************************************************************************
 * wrs/par_scan_alias.hpp
 *
 * Parallel 1-alias table construction
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PAR_SCAN_ALIAS_HEADER
#define WRS_PAR_SCAN_ALIAS_HEADER

#include <wrs/psa/psa_base.hpp>
#include <wrs/timer.hpp>

#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <functional>

namespace wrs {

template <typename size_type = int32_t>
class par_scan_alias : public psa::psa_base<par_scan_alias, size_type> {
    using base = psa::psa_base<par_scan_alias, size_type>;
    // base class types
    using splitter = typename base::Splitter;
    using subproblem = typename base::Subproblem;
    // constexpr definitions
    using base::debug;
    using base::time;
    using base::empty;
    using base::relative_epsilon;

    // member variables
    using base::table_;
    using base::indices_;
    using base::prefsum_;
    using base::classifier_;
    using base::subproblems_;
    using base::splitters_;
    using base::size_;
    using base::num_light_;
    using base::W_;
    using base::W_n_;
    using base::eps_;
    using base::timers_;

    // (non-templated) functions
    using base::is_light;
    using base::get_heavy;
    using base::get_item;
    using base::set_heavy;
    using base::split;
    using base::verify_subproblem_feasibility;
    using base::base_case;
    using base::subproblem_from_splitters;
    using base::fixup_boundary;
public:
    using result_type = size_type;
    static constexpr const char* name = "psa";

    par_scan_alias() : base() {}

    template <typename Iterator>
    par_scan_alias(Iterator begin, Iterator end) {
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

        // Step 0: calculate total weight and weight per bucket
        W_ = parallel_accumulate(begin, end, 0.0);
        W_n_ = W_ / size_;

        eps_ = relative_epsilon * W_;
        LOG << "Using epsilon: " << eps_;

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

        // Step 1: classify
        parallel_do_range(
            [&](size_type min, size_type max, int /* thread */) {
                Iterator it = begin + min;
                for (size_type index = min; index < max; index++, it++) {
                    classifier_[index] = !is_light(*it);
                    table_[index].p = *it;
                    table_[index].a = empty;
                }
            }, size_);
        LOG_ARR(classifier_.get(), "classifier");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 1: classification & table init took " << timers_.back() << "ms";

        // Step 2: prefix sum
        parallel_scan(classifier_.get(), classifier_.get() + size_,
                      classifier_.get(), std::plus<>());
        num_light_ = size_ - classifier_[size_ - 1];
        LOG_ARR(classifier_.get(), "post scan");
        timers_.push_back(t.get_and_reset());
        LOG << "num_light_ = " << num_light_ << " size_ = " << size_;
        LOGC(time) << "Step 2: prefsum took " << timers_.back() << "ms";

        // Step 3: write indices_ array. But first handle the first item in both
        // possible outcomes (the wrong one will be overwritten)
        indices_[0] = 0; indices_[num_light_] = 0;
        prefsum_[0] = *begin; prefsum_[num_light_] = *begin;

        const size_type light_offset = num_light_ - 1;
        parallel_do_range(
            [&](size_type min, size_type max, int /* thread */) {
                size_type prev = classifier_[min - 1];
                size_type loopmax = max - 1, index;
                // manually unrolled because compilers aren't that smart
                for (index = min; index < loopmax; index += 2) {
                    size_type nexti = index + 1,
                        curr = classifier_[index],
                        next = classifier_[nexti],
                        pos1, pos2;
                    if (curr == prev) {
                        assert(index >= curr && index - curr < num_light_);
                        pos1 = index - curr;
                    } else {
                        assert(light_offset + curr < size_);
                        pos1 = light_offset + curr;
                    }
                    if (next == curr) {
                        assert(nexti >= next && nexti - next < num_light_);
                        pos2 = nexti - next;
                    } else {
                        assert(light_offset + next < size_);
                        pos2 = light_offset + next;
                    }
                    assert(0 <= pos1 && pos1 < size_);
                    assert(0 <= pos2 && pos2 < size_);
                    indices_[pos1] = index;
                    prefsum_[pos1] = *(begin + index);
                    indices_[pos2] = nexti;
                    prefsum_[pos2] = *(begin + nexti);
                    prev = next;
                }
                if (index < max) {
                    size_type curr = classifier_[index];
                    size_type pos;
                    if (curr == prev) {
                        assert(index >= curr && index - curr < num_light_);
                        pos = index - curr;
                    } else {
                        assert(light_offset + curr < size_);
                        pos = light_offset + curr;
                    }
                    assert(0 <= pos && pos < size_);
                    indices_[pos] = index;
                    prefsum_[pos] = *(begin + index);
                }
            }, size_, static_cast<size_type>(1));
        LOG_ARR(indices_.get(), "indices");
        LOG_ARR(prefsum_.get(), "L/H weights");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 3: computing indices took " << timers_.back() << "ms";

        // Step 4: prefix sum over prefsum_ weights
        parallel_scan(prefsum_.get(), prefsum_.get() + num_light_,
                      prefsum_.get(), std::plus<>());
        parallel_scan(prefsum_.get() + num_light_, prefsum_.get() + size_,
                      prefsum_.get() + num_light_, std::plus<>());
        LOG_ARR(prefsum_.get(), "L/H prefsum");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 4: weight prefix sum took " << timers_.back() << "ms";

        // Step 5: split up
        LOG_ARR(table_.get(), "initial table");
        int num_threads = wrs::get_num_threads();
        LOG << "have" << num_threads << "threads";
        subproblem problem =
            {0, num_light_, num_light_, size_, size_, num_light_,
             size_ - num_light_, 0.0, 0.0, false, false};

        // 5a: compute splitters
        splitters_[0] = { 0, num_light_, 0.0, 0, false, false };
        splitters_[num_threads] = { num_light_, size_, 0, 0.0, false, false};
        wrs::parallel_do(
            [&](size_type i) {
                double fraction = static_cast<double>(i) / num_threads;
                splitters_[i] = split(problem, fraction);
            }, num_threads, 1);
        LOG_ARRs(splitters_.get(), "splitters", num_threads + 1);
        LOG << "splitting done.\n";

        // 5b: compute subproblems from neighbouring splitters and fix boundaries
        wrs::parallel_do(
            [&](size_type i) {
                subproblems_[i] = subproblem_from_splitters(
                    splitters_[i], splitters_[i + 1]);
                if constexpr(debug)
                    verify_subproblem_feasibility(subproblems_[i]);
                fixup_boundary(subproblems_[i]);
                sLOG << "problem" << i << subproblems_[i];
            }, num_threads);
        LOG << "splitter to problem conversion done.\n";
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 5: splitting took " << timers_.back() << "ms";

        // Step 6: call base cases
        wrs::parallel_do(
            [&](size_type i) {
                sLOG << "\nCalling base case" << i << subproblems_[i];
                base_case(subproblems_[i]);
            }, num_threads);

        LOG_ARR(table_.get(), "table");
        timers_.push_back(t.get_and_reset());
        LOGC(time) << "Step 6: base cases took " << timers_.back() << "ms";
    }

};


} // namespace wrs

#endif // WRS_PAR_SCAN_ALIAS_HEADER
