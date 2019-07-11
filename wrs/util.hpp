/*******************************************************************************
 * wrs/util.hpp
 *
 * Utilities
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_UTIL_HEADER
#define WRS_UTIL_HEADER

#include <tlx/logger.hpp>

#include <algorithm>
#include <atomic>
#include <string>
#include <vector>

namespace wrs {

template <typename Array>
void log_arr(const bool debug, const Array &arr, size_t size,
             const std::string &desc, size_t max_log_size = 100) {
    using value_type = std::remove_reference_t<decltype(arr[0])>;
    if (size > max_log_size) return;
    LOG << desc << ": "
        << std::vector<value_type>(arr, arr + std::min(max_log_size, size))
        << (max_log_size < size ? " (truncated)" : "");
}

#define LOG_ARR(arr, desc) log_arr(debug, (arr), size_, (desc))
#define LOG_ARRs(arr, desc, size) log_arr(debug, (arr), (size), (desc))

// for types that don't have std::atomic<T>::fetch_add (mostly double)
template <typename T>
T atomic_fetch_add(std::atomic<T> *obj, T difference) {
  T expected = obj->load();
  while(!atomic_compare_exchange_weak(obj, &expected, expected + difference));
  return expected;
}

template <typename T>
T atomic_fetch_sub(std::atomic<T> *obj, T difference) {
  T expected = obj->load();
  while(!atomic_compare_exchange_weak(obj, &expected, expected - difference));
  return expected;
}

// technically, this is an unsafe hack :)
template <typename T>
T atomic_fetch_sub(T * ptr, T difference) {
    auto hack = &reinterpret_cast<std::atomic<T>&>(*ptr);
    return atomic_fetch_sub(hack, difference);
}

}  // namespace wrs

#endif // WRS_UTILI_HEADER
