/*******************************************************************************
 * tests/alias_edgecase.cpp
 *
 * Alias table tests: edge case distributions
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

// Make sure NDEBUG isn't set
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <tests/alias_common.hpp>

#include <wrs/generators/stl.hpp>

#include <gtest/gtest.h>

#include <vector>


// Run tests for a specific algorithm-input-combination
template <typename alias_t, typename generator_t>
void alias_test_distribution(generator_t&& gen) {
    static constexpr bool debug = false;

    std::vector<size_t> sizes = {10, 42, 100, 1000, 10000, 100000};
    std::vector<double> weights;
    for (size_t size : sizes) {
        weights.resize(size);
        auto sum = gen(weights, size);

        alias_t alias;
        construct(alias, weights);
        EXPECT_NEAR(alias.total_weight(), sum, 1e-9 * size);
        ASSERT_NO_THROW(wrs::verify(weights.begin(), weights.end(), alias))
            << alias_t::name << " verify fail with " << generator_t::name
            << " generator, size " << size;
    }
    sLOG << "tests for" << alias_t::name << "with" << generator_t::name << "okay!";
}


// Run all inputs for a particular construction algorithm
template <typename alias_t>
struct alias_test_edgecases {
    void operator()() {
        alias_test_distribution<alias_t>(one_heavy_gen());
        alias_test_distribution<alias_t>(one_light_gen());
        alias_test_distribution<alias_t>(one_heavy2_gen());
        alias_test_distribution<alias_t>(powerlaw_gen(0.1));
        alias_test_distribution<alias_t>(powerlaw_gen(1.0));
        alias_test_distribution<alias_t>(powerlaw_gen(2.0));
        alias_test_distribution<alias_t>(powerlaw_gen(3.0));
        alias_test_distribution<alias_t>(powerlaw_shuffle_gen(1.0));
        alias_test_distribution<alias_t>(powerlaw_shuffle_gen(2.0));
    }
};

TEST(alias, edgecase) {
    run_for_all<alias_test_edgecases>()();
}
