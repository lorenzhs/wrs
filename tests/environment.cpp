/*******************************************************************************
 * tests/environment.cpp
 *
 * Set up and tear down global test environment
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

#include <wrs/parallel_do.hpp>

#include <gtest/gtest.h>

// for stocc
#include <mersenne.cpp>
#include <stoc1.cpp>
#include <userintf.cpp>

class WRSEnvironment : public ::testing::Environment {
public:
    void SetUp() {
        wrs::init_threads(4);
    }

    void TearDown() {
        wrs::release_threads();
    }
};


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new WRSEnvironment);
    return RUN_ALL_TESTS();
}
