/*******************************************************************************
 * wrs/outsens.hpp
 *
 * Divide-and-conquer tree for sampling with replacement
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_OUTSENS_HEADER
#define WRS_OUTSENS_HEADER

#include <wrs/dctree.hpp>
#include <wrs/generators/select.hpp>
#include <wrs/generators/stocc.hpp>
#include <wrs/memory.hpp>
#include <wrs/parallel_do.hpp>
#include <wrs/prefix_sum.hpp>
#include <wrs/util.hpp>

#include <tlx/math.hpp>
#include <tlx/logger.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace wrs {

template <typename size_type, bool dedup, size_t bc_size>
class outsens {
    static constexpr bool debug = false;
public:
    static constexpr const char* name = "seqoutsens";
    static constexpr bool yields_single_sample = false;
    static constexpr bool init_with_seed = true;
    using result_type = size_type;
    using subtree = wrs::dctree<size_type, dedup, bc_size>;

    explicit outsens(size_t seed = 0) : size_(-1), num_groups_(-1), W_(0)
                                      , stoc_(seed), rng_(seed + 1) {}

    template <typename Iterator>
    outsens(Iterator begin, Iterator end, size_type id_offset = 0, size_t seed = 0)
        : stoc_(seed), rng_(seed + 1)
    {
        init(end - begin, seed);
        construct(begin, end, id_offset);
    }

    outsens & operator = (const outsens & other) {
        LOG0 << "outsens copy assignment/constructor called";
        size_ = other.size_;
        num_groups_ = other.num_groups_;
        wmin_ = other.wmin_;
        wmax_ = other.wmax_;
        group_tree_ = other.group_tree_;
        if (num_groups_ > 0 && num_groups_ != static_cast<size_t>(-1)) {
            groups_ = std::make_unique<subtree[]>(num_groups_);
            // copy groups
            for (size_t i = 0; i < num_groups_; i++) {
                groups_[i] = other.groups_[i];
            }

            groupsize_ = std::make_unique<size_type[]>(num_groups_);
            memcpy(groupsize_.get(), other.groupsize_.get(), num_groups_ * sizeof(size_type));
        } else {
            groups_ = nullptr;
            groupsize_ = nullptr;
        }
        stoc_ = other.stoc_;
        // seed rng_ from stoc to get a deterministic result
        rng_.seed(stoc_.Random() * std::numeric_limits<size_t>::max());
        return *this;
    }
    outsens(const outsens &other) : stoc_(other.stoc_)
                                    // seed from stoc to get a deterministic result
                                  , rng_(stoc_.Random() * std::numeric_limits<size_t>::max()) {
        *this = other;
    }
    //! delete move-constructor
    outsens(outsens &&) = delete;
    //! delete move-assignment
    outsens & operator = (outsens &&) = delete;

    void init(size_t size, size_t seed) {
        size_ = size;
        stoc_.RandomInit(seed);
        rng_.seed(seed + 1);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end, size_type id_offset = 0) {
        if (end - begin != static_cast<ssize_t>(size_)) {
            sLOG1 << "Error: tried to call construct() with incorrect size!"
                  << "Expected" << size_ << "got" << end - begin;
            return;
        }
        do_init(begin, end);
        if (num_groups_ == 0 || num_groups_ == static_cast<size_t>(-1)) {
            sLOG1 << "Error: cannot construct outsens sampler with"
                  << num_groups_ << "groups! There's probably a zero in the"
                  << "input: wmin =" << wmin_;
            return;
        }

        // Sort into groups. Sequential prototype.
        auto group_for_weight = [&](const double &weight) {
            int group = std::ilogb(weight / wmin_);
            assert(group >= 0 && static_cast<size_t>(group) <= num_groups_);
            return group;
        };
        auto nextleafid = std::make_unique<size_type[]>(num_groups_);
        auto groupweight = std::make_unique<double[]>(num_groups_);
        std::fill(nextleafid.get(), nextleafid.get() + num_groups_, 0);
        std::fill(groupsize_.get(), groupsize_.get() + num_groups_, 0);
        std::fill(groupweight.get(), groupweight.get() + num_groups_, 0);
        for (Iterator it = begin; it != end; ++it) {
            size_t group = group_for_weight(*it);
            groupsize_[group]++;
            groupweight[group] += *it;
        }
        LOG_ARRs(groupsize_.get(), "group sizes", num_groups_);
        LOG_ARRs(groupweight.get(), "group weights", num_groups_);

        LOG << "Initialising groups";
        W_ = 0;
        for (size_t i = 0; i < num_groups_; ++i) {
            groups_[i].init(groupsize_[i]);
            W_ += groupweight[i];
        }

        LOG << "Assigning elements to their groups";
        size_type id = id_offset;
        for (Iterator it = begin; it != end; ++it, ++id) {
            size_t group = group_for_weight(*it);
            groups_[group].set_leaf(nextleafid[group]++, id, *it);
        }

        LOG << "Constructing dc-trees for groups";
        for (size_t i = 0; i < num_groups_; ++i) {
            groups_[i].construct();
            assert(std::abs(groups_[i].total_weight() - groupweight[i]) <
                   1e-9 * std::max((size_type)1, groupsize_[i]));
        }

        LOG << "Constructing dc-tree over groups";
        group_tree_.construct(groupweight.get(), groupweight.get() + num_groups_);
    }

    template <typename Callback>
    void sample(Callback && callback, size_t num_samples) const {
        // assign without early aborting because group weights can be *very*
        // different, rejection sampling could take forever
        group_tree_.assign(
            num_samples, [&](size_t group, size_t samples) {
                sLOG << "Group" << group << "got" << samples << "samples";
                groups_[group].assign_early_abort(samples, callback, stoc_, rng_);
            }, stoc_);
    }

    // Sample with explicit RNG
    template <typename Callback>
    void sample(StochasticLib1 &rng, Callback && callback, size_t num_samples) const {
        // assign without early aborting because group weights can be *very*
        // different, rejection sampling could take forever
        group_tree_.assign(
            num_samples, [&](size_t group, size_t samples) {
                sLOG << "Group" << group << "got" << samples << "samples";
                groups_[group].assign_early_abort(samples, callback, rng, stoc_);
            }, stoc_);
    }

    double total_weight() const {
        return W_;
    }

    size_t size() const {
        return size_;
    }

    // for internal use only
    void verify_helper(std::vector<double> &weights, size_t) const {
        if (size_ == static_cast<size_t>(-1)) return;
        for (size_t i = 0; i < num_groups_; ++i) {
            size_t num_leaves = groupsize_[i];
            for (size_t j = 0; j < num_leaves; ++j) {
                auto [id, w] = groups_[i].get_leaf(j);
                weights[id] += w;
            }
        }
    }

    template <typename Callback>
    void find(size_type item, Callback && callback) const {
        for (size_t group = 0; group < num_groups_; ++group) {
            auto cb = [&](int, size_type idx, double w, auto bucket)
                      { callback(group, idx, w, bucket); };
            groups_[group].find(item, cb);
        }
    }

protected:

    template <typename Iterator>
    void do_init(Iterator begin, Iterator end) {
        auto [minit, maxit] = std::minmax_element(begin, end);
        wmin_ = *minit;
        wmax_ = *maxit;
        assert(wmin_ > 0);
        num_groups_ = std::log2(wmax_ / wmin_) + 1;

        sLOG << "outsens:" << size_ << "elements with" << num_groups_
             << "groups, wmin:" << wmin_ << "wmax:" << wmax_;

        group_tree_.init(num_groups_);

        groups_ = std::make_unique<subtree[]>(num_groups_);
        groupsize_ = std::make_unique<size_type[]>(num_groups_);
    }

    wrs::dctree<size_type, false, 1> group_tree_;
    std::unique_ptr<subtree[]> groups_;
    std::unique_ptr<size_type[]> groupsize_;
    size_t size_;
    size_t num_groups_;
    double wmin_, wmax_, W_;

    mutable StochasticLib1 stoc_;
    mutable generators::select_t rng_;
};


template <typename size_type = int32_t, bool dedup = true, size_t bc_size = 128>
class par_outsens {
    static constexpr bool debug = false;
public:
    static constexpr const char* name = "paroutsens";
    static constexpr bool yields_single_sample = false;
    static constexpr bool init_with_seed = true;
    using result_type = size_type;
    using outsens = wrs::outsens<size_type, dedup, bc_size>;

    par_outsens() {   }

    template <typename Iterator>
    par_outsens(Iterator begin, Iterator end, size_t seed) {
        init(end - begin, seed);
        construct(begin, end);
    }
    par_outsens & operator = (const par_outsens & other) {
        LOG0 << "par_outsens copy assignment/constructor called";
        size_ = other.size_;
        count_ = other.count_;
        samplers_ = std::make_unique<outsens[]>(count_);

        tree_ = other.tree_;
        for (size_t i = 0; i < count_; i++) {
            // copy samplers
            samplers_[i] = other.samplers_[i];
        }

        return *this;
    }
    par_outsens(const par_outsens &other) {
        *this = other;
    }
    //! delete move-constructor
    par_outsens(par_outsens &&) = delete;
    //! delete move-assignment
    par_outsens & operator = (par_outsens &&) = delete;


    void init(size_t size, size_t seed) {
        size_ = size;
        count_ = get_num_threads();
        samplers_ = std::make_unique<outsens[]>(count_);

        tree_.init(count_);

        CRandomMersenne rng(seed);
        seeds_ = std::make_unique<size_t[]>(count_);
        for (size_t i = 0; i < count_; i++) {
            seeds_[i] = rng.Random() * std::numeric_limits<size_t>::max();
        }
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        if (end - begin != static_cast<ssize_t>(size_)) {
            sLOG1 << "Error: tried to call construct() with incorrect size!"
                  << "Expected" << size_ << "got" << end - begin;
            return;
        }

        LOG << "Constructing subsamplers";
        auto subconstruct =
            [&](size_t min, size_t max, int thread) {
                sLOG << "Thread" << thread << "constructing table for"
                     << max - min << "items, range" << min << "to" << max;
                samplers_[thread].init(max - min, seeds_[thread]);
                samplers_[thread].construct(begin + min, begin + max, min);
            };
        parallel_do_range(subconstruct, size_);

        LOG << "Constructing dc-tree over subsamplers";
        for (size_t i = 0; i < count_; i++) {
            sLOG << "subsampler" << i << "has weight" << samplers_[i].total_weight();
            tree_.set_leaf(i, i, samplers_[i].total_weight());
        }
        tree_.construct();
        LOG << "subsampler tree:" << tree_;
    }

    template <typename Callback>
    size_t sample(StochasticLib1 &rng, Callback && callback, size_t num_samples) {
        thread_local size_t local_count; // initialised later
        std::atomic<size_t> global_count = 0;
        // Wrap callback to keep track
        auto cb = [&](auto sample, auto mult) {
            callback(sample, mult);
            ++local_count;
        };
        tree_.assign(
            num_samples, [&](size_t id, size_t samples) {
                sLOG << "Subsampler" << id << "got" << samples << "samples";
                wrs::do_for_thread(
                    [&, id, samples]() {
                        wrs::timer t;
                        local_count = 0;
                        samplers_[id].sample(cb, samples);
                        global_count.fetch_add(local_count);
                        sLOG << "Thread" << id << "took" << t.get() << "ms for"
                             << samples << "samples, got" << local_count
                             << "unique elements";
                    }, id);
            }, rng);
        wrs::wait_all_threads();
        sLOG << "Overall," << global_count << "of" << num_samples << "were unique";
        return global_count;
    }

    std::vector<double> get_timers() const {
        return {};
    }

    size_t size() const {
        return size_;
    }

    double total_weight() const {
        double w = 0;
        for (size_t i = 0; i < count_; i++) {
            w += samplers_[i].total_weight();
        }
        return w;
    }

    void verify_helper(std::vector<double> &weights, size_t) const {
        for (size_t i = 0; i < count_; i++) {
            samplers_[i].verify_helper(weights, /* dummy */ 0);
        }
    }

    template <typename Callback>
    void find(size_type item, Callback && callback) const {
        for (size_t i = 0; i < count_; i++) {
            auto cb = [&](int group, auto idx, double w, auto &bucket)
                      { callback(std::make_pair(i, group), idx, w, bucket); };
            samplers_[i].find(item, cb);
        }
    }

protected:
    std::unique_ptr<outsens[]> samplers_;
    std::unique_ptr<size_t[]> seeds_;
    wrs::dctree<size_type, false, 1> tree_;
    size_t size_, count_;
};

} // namespace wrs

#endif // WRS_OUTSENS_HEADER
