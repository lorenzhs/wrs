/*******************************************************************************
 * wrs/verify.hpp
 *
 * Alias table verification helpers
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_VERIFY_HEADER
#define WRS_VERIFY_HEADER

#include <tlx/logger/all.hpp>

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace wrs {

// Sum up the weights in a [probability, alias] alias table, adding `offset` to
// every item ID
template <typename size_type, typename tableitem>
void sum_table(const tableitem *table, size_t size, double total_weight,
               double bucket_size, std::vector<double> &weights, size_t offset)
{
    if (size == static_cast<size_t>(-1)) return;
    const  __attribute__((unused)) double max_dev = 1e-10 * total_weight;
    for (size_t i = 0; i < size; i++) {
        auto [p, a] = table[i];
        if (p <= 0) LOG1 << "p <= 0 at position " << i;
        assert(p > 0);

        weights[offset + i] += p;
        if (p < bucket_size) {
            if (a == static_cast<size_type>(-1)) {
                double dev = std::abs(bucket_size - p);
                if (dev < 1e-10) continue; // ignore small deviations
                sLOG1 << "no alias despite p < W/n at position" << i
                      << "diff" << dev;
                assert(a != static_cast<size_type>(-1) || dev < max_dev);
            } else {
                weights[offset + a] += bucket_size - p;
            }
        }
    }

}

// Sum up all weights in a [probability, item, alias] alias table, adding
// `offset` to every item ID.
template <typename size_type, typename tableitem>
void sum_tritable(const tableitem* table, size_t size, double total_weight,
                  double bucket_size, std::vector<double> &weights,
                  size_t offset) {
    if (size == static_cast<size_t>(-1)) return;
    const __attribute__((unused)) double max_dev = 1e-10 * total_weight;
    for (size_t i = 0; i < size; i++) {
        auto [p, i1, i2] = table[i];
        if (p <= 0) LOG1 << "p <= 0 at position " << i;
        assert(p > 0);

        assert(i1 != static_cast<size_type>(-1));
        weights[offset + i1] += p;
        if (p < bucket_size) {
            if (i2 == static_cast<size_type>(-1)) {
                double dev = std::abs(bucket_size - p);
                if (dev < 1e-10) continue; // ignore small deviations
                sLOG1 << "no alias despite p < W/n at position" << i
                      << "diff" << dev;
                assert(i2 != static_cast<size_type>(-1) || dev < max_dev);
            } else {
                weights[offset + i2] += bucket_size - p;
            }
        }
    }
}

// Check that `weights` is close to the weights that the iterator range
// [i_begin, i_end) dereferences to
template <typename Iterator, typename alias_t>
void check_weights(Iterator i_begin,  __attribute__((unused)) Iterator i_end,
                   const std::vector<double> &weights,
                   const alias_t &alias)
{
    assert(i_end - i_begin == static_cast<ssize_t>(weights.size()));
    const double total_weight = alias.total_weight();
    const double max_dev = total_weight * 1e-10;
    size_t i, last = static_cast<size_t>(-1);
    double sum = 0.0; // last and sum are for the find() callback
    auto cb = [&](auto table, auto idx, double weight, const auto & bucket) {
                  if (last == i) { sum += weight; }
                  else { sum = weight; last = i; }
                  sLOG1 << "found item" << i << "in table" << table
                        << "index" << idx << "with weight" << weight
                        << "bucket:" << bucket
                        << "nominal weight:" << weights[i] << "diff:"
                        << weights[i] - sum;
              };
    // GSL has much higher deviations, don't check for low relative error
    const bool skip_rel_err_check = (std::string(alias_t::name) == "gsl");
    for (i = 0; i < weights.size(); ++i) {
        double should = *(i_begin + i);
        double have = weights[i];
        double dev = std::abs(should - have);
        if (have > 1e-16 && !skip_rel_err_check) {
            double rel_err = should / have - 1;
            double max_err = (have > 1e-10) ?
                // large tables are allowed some more errors
                (alias.size() >= 50'000'000 ? 1e-3 : 1e-6) : 1e-4;
            if (std::abs(rel_err) >= max_err) {
                sLOG1 << "Item" << i << "error" << rel_err << ">=" << max_err
                      << "should" << should << "have" << have << "quot"
                      << should / have;
                alias.find(i, cb);
            }
            assert(std::abs(rel_err) < max_err);
        }
        if (dev > max_dev) {
            sLOG1 << "Item" << i << "deviating by" << dev << ">" << max_dev
                  << "should" << should << "have" << have << "diff"
                  << should - have;
            alias.find(i, cb);
        }
        // GSL does weird things for index 0, or maybe we're interpreting its
        // data structures wrong there? not sure
        assert(dev <= max_dev || (skip_rel_err_check && i == 0));
    }
    sLOG1 << alias_t::name << "verification succeeded!";
}

// Verify an alias table (not for keyed versions)
template <typename Iterator, typename alias_t>
void verify(Iterator begin, Iterator end, const alias_t &alias,
            const bool debug = false)
{
        (void) begin;
        (void) end;
        (void) alias;
        (void) debug;
#ifndef NDEBUG
        assert(end - begin == static_cast<ssize_t>(alias.size()));
        std::vector<double> weights(alias.size());
        alias.verify_helper(weights, 0);
        if (alias.size() < 100) {
            LOG << "reconstructed: " << weights;
        }
        check_weights(begin, end, weights, alias);
#endif
}


} // namespace wrs

#endif // WRS_VERIFY_HEADER
