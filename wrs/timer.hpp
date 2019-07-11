/*******************************************************************************
 * wrs/timer.hpp
 *
 * Copyright (C) 2018-2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_TIMER_HEADER
#define WRS_TIMER_HEADER

#include <chrono>

namespace wrs {

/**
 * A flexible timer.
 *
 * resolution is the precision of the timing, while scaling_factor
 * is the factor by which the output will be scaled. The default is to
 * return milliseconds with microsecond precision.
 */
template <typename resolution = std::chrono::microseconds,
          int scaling_factor = 1000, typename return_type = double>
struct base_timer {
    base_timer() {
        reset();
    }

    void reset() {
        start = std::chrono::system_clock::now();
    }

    return_type get() const {
        auto duration = std::chrono::duration_cast<resolution>(
            std::chrono::system_clock::now() - start);
        return static_cast<return_type>(duration.count()) / scaling_factor;
    }

    return_type get_and_reset() {
        auto t = get();
        reset();
        return t;
    }

private:
    std::chrono::system_clock::time_point start;
};

/// A timer that is accurate to microseconds, formatted as milliseconds
typedef base_timer<std::chrono::microseconds, 1000, double> timer;
/// A timer that is accurate to milliseconds, formatted as seconds
typedef base_timer<std::chrono::milliseconds, 1000, double> sec_timer;

} // namespace wrs

#endif // WRS_TIMER_HEADER
