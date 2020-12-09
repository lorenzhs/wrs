/*******************************************************************************
 * tests/alias_verify.cpp
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
struct alias_verify {
    static constexpr bool debug = false;

    void operator()() {
        std::vector<double> weights = {1, 2, 0.5, 0.1, 0.9, 1.5};
        alias_t a1;
        construct(a1, weights);
        EXPECT_EQ(a1.size(), weights.size()) << alias_t::name << " size fail";
        EXPECT_NEAR(a1.total_weight(), 6, 1e-9) << alias_t::name << " weight fail";
        ASSERT_NO_THROW(wrs::verify(weights.begin(), weights.end(), a1))
            << alias_t::name << " verify fail";

        double sum = uniform_gen(42)(weights, 1000);
        alias_t a2;
        construct(a2, weights);
        EXPECT_EQ(a2.size(), weights.size()) << alias_t::name << " size fail";
        EXPECT_NEAR(a2.total_weight(), sum, 1e-8) << alias_t::name << " weight fail";
        ASSERT_NO_THROW(wrs::verify(weights.begin(), weights.end(), a2))
            << alias_t::name << " verify fail";
    }
};


TEST(alias, verify) {
    run_for_all<alias_verify>()();
}
