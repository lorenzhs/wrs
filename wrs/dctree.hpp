/*******************************************************************************
 * wrs/dctree.hpp
 *
 * Divide-and-conquer tree for sampling with replacement
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_DCTREE_HEADER
#define WRS_DCTREE_HEADER

#include <wrs/generators/stocc.hpp>
#include <wrs/memory.hpp>
#include <wrs/tinyhashtable.hpp>

#include <tlx/math.hpp>
#include <tlx/logger.hpp>

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace wrs {

template <typename size_type, bool dedup, size_t bc_size>
class dctree {
    static constexpr bool debug = false;
public:
    static constexpr size_t abort_count = bc_size;
    static constexpr size_t cache_size = dedup ? bc_size : 0;

    dctree() : weights_(nullptr), cache_(cache_size)
             , num_leaves_(-1), num_inners_(-1), size_(-1), levels_(-1) {}

    template <typename Iterator>
    dctree(Iterator begin, Iterator end) : cache_(cache_size) {
        init(end - begin);
        construct(begin, end);
    }

    dctree & operator = (const dctree & other) {
        LOG0 << "dctree copy assignment/constructor called";
        max_weight_ = other.max_weight_;
        num_leaves_ = other.num_leaves_;
        num_inners_ = other.num_inners_;
        size_ = other.size_;
        levels_ = other.levels_;
        if (size_ > 0)
            weights_ = copy_alloc_arr<double>(other.weights_.get(), size_);
        else
            weights_ = nullptr;
        if (num_leaves_ > 0)
            ids_ = copy_alloc_arr<size_type>(other.ids_.get(), num_leaves_);
        else
            ids_ = nullptr;
        return *this;
    }
    dctree(const dctree &other) : cache_(abort_count) {
        *this = other;
    }
    //! delete move-constructor
    dctree(dctree &&) = delete;
    //! delete move-assignment
    dctree & operator = (dctree &&) = delete;

    void init(size_type num_leaves) {
        num_leaves_ = num_leaves;
        levels_ = tlx::integer_log2_ceil(num_leaves);
        num_inners_ = (size_type{1} << levels_) - 1;
        // + 1 because 0 is a dummy
        size_ = num_leaves_ + num_inners_ + 1;
        weights_ = make_alloc_arr<double>(size_);
        ids_ = make_alloc_arr<size_type>(num_leaves_);

        sLOG << "dctree with" << num_leaves_ << "leaves," << num_inners_
             << "inner nodes and" << levels_ << "levels -> size" << size_;
    }

    void set_leaf(size_type leaf_id, size_type id, double weight) {
        assert(leaf_id < num_leaves_);
        weights_[num_inners_ + leaf_id + 1] = weight;
        ids_[leaf_id] = id;
    }

    std::pair<size_type, double> get_leaf(size_type leaf_id) const {
        assert(leaf_id < num_leaves_);
        return std::make_pair(ids_[leaf_id],
                              weights_[num_inners_ + leaf_id + 1]);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end, size_type id_offset = 0) {
        assert(end - begin == static_cast<ssize_t>(num_leaves_));
        // Step 1: copy weights
        Iterator in_it = begin;
        auto out_it = weights_.get() + num_inners_ + 1; // +1 for dummy
        while (in_it != end) {
            *out_it++ = *in_it++;
        }
        // Step 2: write IDs
        for (size_type i = 0; i < num_leaves_; i++) {
            ids_[i] = i + id_offset;
        }

        construct();
    }

    // Construct method if leaf weights are already there
    void construct() {
        if (num_leaves_ == 0) {
            weights_[0] = 0;
            // empty
            return;
        }
        // find max weight
        max_weight_ = *std::max_element(weights_.get() + num_inners_ + 1,
                                       weights_.get() + size_);

        // fill lowest level with zeroes where needed. + 1 because for dummy
        for (size_type i = parent(size_); i < num_inners_ + 1; i++) {
            LOG << "filling 0 for childless node " << i;
            weights_[i] = 0;
        }

        size_type min = num_inners_ + 1, max = size_;
        while (max > 1) { // abort once only dummy is left
            sLOG << "processing range" << min << max;
            size_type i;
            for (i = min; i < max - 1; i += 2) {
                weights_[parent(i)] = weights_[i] + weights_[i+1];
                LOG << "weights_[" << parent(i) << "] = weights_[" << i
                    << "] + weights_[" << i+1 << "] = "
                    << weights_[i] + weights_[i+1];
            }
            if (i < max) { // level not full
                sLOG << "level not full, copying" << i << "to" << parent(i)
                     << "weight" << weights_[i];
                weights_[parent(i)] = weights_[i];
                i++;
                assert(i == max);
            }
            min = parent(min);
            // need to round up because num_leaves_ need not be a power of two
            max = parent(tlx::round_up_to_power_of_two(max));
        }

        if (size_ < 100)
            LOG << std::vector<double> { weights_.get(), weights_.get() + size_ };
    }

    size_type leaves_ptr(double** out) {
        // + 1 because of the dummy node (0)
        *out = weights_.get() + num_inners_ + 1;
        return num_leaves_;
    }

    const double& total_weight() const {
        return weights_[0];
    }

    size_type size() const {
        return size_;
    }

    const double& get(size_type index) const {
        return weights_[index];
    }

    template <typename Callback>
    void assign(size_t count, Callback && callback, StochasticLib1 &stoc) const {
        sLOG << "Assigning" << count << "samples in dctree of"
             << num_inners_ << "inner nodes and" << num_leaves_ << "leaves";
        // Pass stoc as dummy base case RNG
        assign_recursive<false>(count, 1, callback, stoc, stoc);
    }

    // Use stoc as base case RNG, too
    template <typename Callback, typename RNG>
    void assign_early_abort(size_t count, Callback && callback,
                            StochasticLib1 &stoc) const
    {
        sLOG << "Assigning" << count << "samples with early abort in dctree of"
             << num_inners_ << "inner nodes and" << num_leaves_ << "leaves";
        assign_recursive<true>(count, 1, callback, stoc, stoc);
    }

    // Use explicit base case RNG
    template <typename Callback, typename RNG>
    void assign_early_abort(size_t count, Callback && callback,
                            StochasticLib1 &stoc, RNG &rng) const
    {
        sLOG << "Assigning" << count << "samples with early abort in dctree of"
             << num_inners_ << "inner nodes and" << num_leaves_ << "leaves";
        assign_recursive<true>(count, 1, callback, stoc, rng);
    }

    friend std::ostream &operator << (std::ostream &os, const dctree &t) {
        os << "dctree(size=" << t.size_ << " inners=" << t.num_inners_
           << " leaves=" << t.num_leaves_ << " => " << t.levels_ << " levels"
           << ", maxw=" << t.max_weight_ << ", weights = ";
        tlx::LoggerFormatter<std::vector<double>>::print(
            os, std::vector<double>(
                t.weights_.get(),
                t.weights_.get() + std::min<size_type>(t.size_, 100)));
        os << ", ids = ";
        tlx::LoggerFormatter<std::vector<size_type>>::print(
            os, std::vector<size_type>(
                t.ids_.get(),
                t.ids_.get() + std::min<size_type>(t.num_leaves_, 50)));
        return os << ")";
    }


    template <typename Callback>
    void find(size_type item, Callback && callback) const {
        if (size_ == static_cast<size_type>(-1)) return;
        for (size_type i = 0; i < num_leaves_; i++) {
            if (ids_[i] == item) {
                callback(0, i, weights_[num_inners_ + i + 1], max_weight_);
            }
        }
    }


protected:
    // Assign a single element with rejection sampling
    template <typename Callback, typename RNG>
    void assign_one(size_type min_id, size_type max_id,
                    Callback && callback, RNG &rng) const
    {
        LOG << "assign_one(" << min_id << ", " << max_id << ") size_=" << size_;
        assert(max_id < size_);
        size_type range = max_id - min_id + 1;
        double rand;
        size_type id;
        do {
            rand = rng.next() * range;
            id = rand;
            rand = (rand - id) * max_weight_;
            id += min_id;
            sLOG << "Rejection sampling: trying item" << id + 1 - min_id << "of" << range
                 << "item weight" << weights_[ id] << "rand=" << rand
                 << "max_weight=" << max_weight_
                 << (rand <= weights_[id] ? "ACCEPTED" : "REJECTED")
                 << "leaf" << leaf_id(id);
        } while (rand > weights_[id]);
        callback(ids_[leaf_id(id)], 1);
    }

    // Assign several elements with rejection sampling
    template <typename Callback, typename RNG>
    void assign_few(size_t count, size_type min_id, size_type max_id,
                    Callback && callback, RNG &rng) const
    {
        LOG << "assign_few(" << count << ", " << min_id << ", " << max_id
            << ") size_=" << size_;
        assert(max_id < size_);
        assert(is_leaf(min_id) && is_leaf(max_id));
        assert(min_id < max_id);
        size_type range = max_id - min_id + 1;

        if constexpr (dedup) {
            cache_.clear();
        }
        for (size_t i = 0; i < count; i++) {
            double rand;
            size_type id;
            do {
                rand = rng.next() * range;
                id = rand;
                rand = (rand - id) * max_weight_;
                id += min_id;
                assert(0 <= rand && rand <= max_weight_);
                sLOG << "Rejection sampling: trying item" << id + 1 - min_id << "of" << range
                     << "item weight" << weights_[id] << "rand=" << rand
                     << "max_weight=" << max_weight_
                     << (rand <= weights_[id] ? "ACCEPTED" : "REJECTED")
                     << "leaf" << leaf_id(id);
            } while (rand > weights_[id]);
            if constexpr (dedup) {
                cache_[ids_[leaf_id(id)]] += 1;
            } else {
                callback(ids_[leaf_id(id)], 1);
            }
        }
        if constexpr (dedup) {
            cache_.foreach(callback);
        }
    }

    template <bool early_abort, typename Callback, typename RNG>
    void assign_recursive(
        size_t count, size_type index, Callback && callback,
        StochasticLib1 &stoc, RNG &rng) const
    {
        LOG << "assign<" << early_abort << ">(" << count << ", " << index << ")";

        if (is_leaf(index)) {
            callback(ids_[leaf_id(index)], count);
            return;
        }

        if constexpr(early_abort) {
            if (count <= abort_count) {
                size_type level = tlx::integer_log2_floor(index);
                size_type maxchildren = size_type(1) << (levels_ - level);
                if (count <= 2 * static_cast<size_t>(maxchildren)) {
                    // Compute leaf range
                    size_type min_id = index << (levels_ - level);
                    size_type max_id = min_id | (maxchildren - 1);
                    if (max_id >= size_) max_id = size_ - 1;
                    sLOG << "Early abort at node" << index << "with leaf range"
                         << min_id << "to" << max_id;
                    if (max_id == min_id) {
                        // we hit a subtree with only one leaf!
                        callback(ids_[leaf_id(min_id)], count);
                    } else if (count > 1) {
                        assign_few(count, min_id, max_id, callback, rng);
                    } else {
                        assign_one(min_id, max_id, callback, rng);
                    }
                    return;
                }
            }
        }

        size_type left_id = left(index), right_id = right(index);
        if (right_id == size_) {
            // Right child does not exist
            assign_recursive<early_abort>(count, left_id, callback, stoc, rng);
            return;
        }
        double left_w = weights_[left_id], right_w = weights_[right_id];

        if (left_w < 1e-10) {
            sLOG << "\tleft subtree of" << index << "is empty, giving all"
                 << count << "samples to" << right_id << "of weight" << right_w;
            if (is_leaf(right_id)) callback(ids_[leaf_id(right_id)], count);
            else assign_recursive<early_abort>(
                count, right_id, callback, stoc, rng);
            return;
        } else if (right_w < 1e-10) {
            sLOG << "\tright subtree of" << index << "is empty, giving all"
                 << count << "samples to" << left_id << "of weight" << left_w;
            if (is_leaf(left_id)) callback(ids_[leaf_id(left_id)], count);
            else assign_recursive<early_abort>(
                count, left_id, callback, stoc, rng);
            return;
        }

        double p_left = left_w / (left_w + right_w);
        size_t num_left = stoc.Binomial(count, p_left);
        sLOG << "Of" << count << "items," << num_left
             << "come from the left subtree of" << index << "left weight"
             << left_w << "vs right" << right_w << "=" << p_left*100 << "%";

        if (num_left > 0)
            assign_recursive<early_abort>(
                num_left, left_id, callback, stoc, rng);
        if (num_left < count)
            assign_recursive<early_abort>(
                count - num_left, right_id, callback, stoc, rng);
    }


    alloc_arr_ptr<double> weights_;
    alloc_arr_ptr<size_type> ids_;
    mutable tinyhashtable<size_type, size_type> cache_;
    double max_weight_;
    size_type num_leaves_;
    size_type num_inners_;
    size_type size_;
    unsigned levels_;

    static constexpr size_type parent(size_type index) { return index/2; }
    static constexpr size_type left(size_type index) { return 2 * index; }
    static constexpr size_type right(size_type index) { return 2 * index + 1; }
    size_type leaf_id(size_type index) const { return index - num_inners_ - 1; }
    bool is_leaf(size_type index) const { return index > num_inners_; }
};

} // namespace wrs

#endif // WRS_DCTREE_HEADER
