/*******************************************************************************
 * wrs/dedup.hpp
 *
 * Sequential de-duplication for sampling with replacement
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_DEDUP_HEADER
#define WRS_DEDUP_HEADER

#include <wrs/generators/select.hpp>

#include <google/dense_hash_map>

#include <vector>

namespace wrs {

template <typename alias_t>
class dedup {
public:
    static constexpr const char* name = "dedup";
    static constexpr bool yields_single_sample = false;
    static constexpr bool init_with_seed = true;
    using result_type = typename alias_t::result_type;
    using size_type = size_t; // for now?

    explicit dedup(size_t seed = 0) : rng_(seed) {}

    template <typename Iterator>
    dedup(Iterator begin, Iterator end, size_t seed = 0) : rng_(seed) {
        alias_table_.init(end - begin);
        alias_table_.construct(begin, end);
    }

    void init(size_t size, size_t seed) {
        rng_.seed(seed);
        if (!flag)
            hash_table_.set_empty_key(static_cast<result_type>(-1));
        flag = true;
        alias_table_.init(size);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        alias_table_.construct(begin, end);
    }

    template <typename Callback>
    size_t sample(Callback&& callback, size_t num_samples) {
        hash_table_.clear_no_resize();
        for (size_t i = 0; i < num_samples; i++) {
            result_type sample = alias_table_.sample(rng_.next());
            hash_table_[sample]++;
        }

        size_t count = hash_table_.size();
        for (auto it = hash_table_.begin(); it != hash_table_.end(); ++it) {
            callback(it->first, it->second);
        }
        return count;
    }

    // interface hack for compatibility with outsens / par_outsens
    template <typename Ignored, typename Callback>
    size_t sample(Ignored, Callback&& callback, size_t num_samples) {
        return sample(callback, num_samples);
    }

    size_t size() const {
        return alias_table_.size();
    }

    double total_weight() const {
        return alias_table_.total_weight();
    }

    std::vector<double> get_timers() const {
        return alias_table_.get_timers();
    }

    void verify_helper(std::vector<double>& weights, size_t offset) const {
        alias_table_.verify_helper(weights, offset);
    }

    template <typename Callback>
    void find(size_type item, Callback&& callback) const {
        alias_table_.find(item, callback);
    }

protected:
    alias_t alias_table_;
    google::dense_hash_map<result_type, size_type> hash_table_;
    wrs::generators::select_t rng_;
    bool flag = false;
};


template <typename alias_t>
class store_vec {
public:
    static constexpr const char* name = "storevec";
    static constexpr bool yields_single_sample = false;
    static constexpr bool init_with_seed = true;
    using result_type = typename alias_t::result_type;
    using size_type = size_t; // for now?

    explicit store_vec(size_t seed = 0) : rng_(seed) {}

    template <typename Iterator>
    store_vec(Iterator begin, Iterator end, size_t seed = 0) : rng_(seed) {
        alias_table_.init(end - begin);
        alias_table_.construct(begin, end);
    }

    void init(size_t size, size_t seed) {
        rng_.seed(seed);
        alias_table_.init(size);
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        alias_table_.construct(begin, end);
    }

    template <typename Callback>
    size_t sample(Callback&& callback, size_t num_samples) {
        samples_.clear();
        samples_.reserve(num_samples);
        for (size_t i = 0; i < num_samples; i++) {
            samples_.push_back(alias_table_.sample(rng_.next()));
        }
        for (size_t i = 0; i < num_samples; i++) {
            callback(samples_[i], 1);
        }
        return samples_.size();
    }

    // interface hack for compatibility with outsens / par_outsens
    template <typename Ignored, typename Callback>
    size_t sample(Ignored, Callback&& callback, size_t num_samples) {
        return sample(callback, num_samples);
    }

    size_t size() const {
        return alias_table_.size();
    }

    double total_weight() const {
        return alias_table_.total_weight();
    }

    std::vector<double> get_timers() const {
        return alias_table_.get_timers();
    }

    void verify_helper(std::vector<double>& weights, size_t offset) const {
        alias_table_.verify_helper(weights, offset);
    }

    template <typename Callback>
    void find(size_type item, Callback&& callback) const {
        alias_table_.find(item, callback);
    }

protected:
    alias_t alias_table_;
    std::vector<result_type> samples_;
    wrs::generators::select_t rng_;
};

} // namespace wrs

#endif // WRS_DEDUP_HEADER
