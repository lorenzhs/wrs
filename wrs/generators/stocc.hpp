/*******************************************************************************
 * wrs/generators_stocc.hpp
 *
 * Wrapper around Agner Fog's stocc non-uniform random number generator
 * available at https://www.agner.org/random/ (GPLv3)
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_GENERATORS_STOCC_HEADER
#define WRS_GENERATORS_STOCC_HEADER

#include <stocc.h>

// don't include the cpps for the tests
#ifndef WRS_IS_TEST

#include <mersenne.cpp>
#include <stoc1.cpp>
#include <userintf.cpp>

#endif

#endif // WRS_GENERATORS_STOCC_HEADER
