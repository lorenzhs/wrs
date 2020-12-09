/*******************************************************************************
 * tests/accumulate.cpp
 *
 * Accumulation test
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

#include <tests/alias_common.hpp>

#include <wrs/accumulate.hpp>

#include <gtest/gtest.h>

#include <numeric>
#include <vector>

TEST(accumulate, random_double) {
    std::vector<size_t> sizes = {1,    2,       12,      31,          32,
                                 33,   100,     512,     1023,        1024,
                                 1025, 31 * 32, 32 * 33, 32 * 32 * 32};
    for (size_t size : sizes) {
        std::vector<double> weights(size);
        uniform_gen()(weights, size);

        double expect = std::accumulate(weights.begin(), weights.end(), 0.0);
        double have = wrs::accumulate(weights.begin(), weights.end(), 0.0);
        // can't use EXPECT_DOUBLE_EQ because the whole point is to be more
        // numerically stable
        EXPECT_NEAR(have, expect, 1e-9);
    }
}

TEST(accumulate, simple_int) {
    std::vector<size_t> sizes = {1,    2,       12,      31,          32,
                                 33,   100,     512,     1023,        1024,
                                 1025, 31 * 32, 32 * 33, 32 * 32 * 32};
    for (size_t size : sizes) {
        std::vector<double> weights(size);
        for (size_t i = 0; i < weights.size(); i++) {
            weights[i] = i;
        }

        double expect = weights.size() * (weights.size() - 1) / 2;
        double have = wrs::accumulate(weights.begin(), weights.end(), 0.0);
        EXPECT_EQ(have, expect);
    }
}
