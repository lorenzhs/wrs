/*******************************************************************************
 * tests/alias_query.cpp
 *
 * Alias table tests
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

// Make sure NDEBUG isn't set
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <tests/alias_common.hpp>

#include <gtest/gtest.h>

#include <vector>


template <typename alias_t>
struct alias_query {
    static constexpr bool debug = false;

    void operator()() {
        std::vector<double> weights;
        size_t size = 64, samples = 50'000 * size;
        double sum = uniform_gen(42)(weights, size);
        alias_t a;
        construct(a, weights);
        EXPECT_EQ(a.size(), weights.size());
        EXPECT_NEAR(a.total_weight(), sum, 1e-8);
        ASSERT_NO_THROW(wrs::verify(weights.begin(), weights.end(), a));

        std::vector<int> freq(size, 0);
        do_sample(a, freq, samples, 12345);

        for (size_t i = 0; i < size; i++) {
            double expected = samples * weights[i] / sum;
            int min = 0.95 * expected, max = 1.05 * expected;

            sLOG << "Expecting item" << i << "of weight" << weights[i] << "="
                 << weights[i] / sum * 100 << "% to have around" << expected
                 << "but (most likely) between" << min << "and" << max << "samples,"
                 << "have" << freq[i];

            EXPECT_LE(min, freq[i])
                << alias_t::name << " min sample count fail for item " << i
                << " wanted MIN " << min << " exp " << expected << " max "
                << max << " prob[i] = " << weights[i] / sum;
            EXPECT_GE(max, freq[i])
                << alias_t::name << " max sample count fail for item " << i
                << " wanted min " << min << " exp " << expected << " MAX "
                << max << " prob[i] = " << weights[i] / sum;
        }
    }
};

TEST(alias, query) {
    run_for_all<alias_query>()();
}
