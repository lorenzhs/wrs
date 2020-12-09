/*******************************************************************************
 * wrs/rng/stl.hpp
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WRS_GENERATORS_SELECT_HEADER
#define WRS_GENERATORS_SELECT_HEADER

#include "dSFMT.hpp"
#include "stl.hpp"
#ifdef WRS_HAVE_MKL
#include "mkl.hpp"
#endif

namespace wrs {
namespace generators {

struct select {
#ifdef WRS_IS_TEST
    using type = dSFMT;
#else
#ifdef WRS_HAVE_MKL
    // MKL is much faster than anything else, by a factor that's not even funny
    // any more for large block sizes
    using type = mkl;
#else
    // dSFMT is at least twice as fast as std::mt19937_64 for large blocks with
    // gcc, and more using clang. It's practically never slower, so prefer it.
    using type = dSFMT;
#endif // WRS_HAVE_MKL
#endif // WRS_IS_TEST
};

using select_t = select::type;

} // namespace generators
} // namespace wrs

#endif // WRS_GENERATORS_SELECT_HEADER
