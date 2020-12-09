/*******************************************************************************
 * tests/alias_common.cpp
 *
 * Alias table tests: common stuff
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 * *****************************************************************************/

#include <tests/alias_common.hpp>

#include <wrs/generators/stl.hpp>


// uniformly random in [0, 1)
double uniform_gen::operator()(std::vector<double> &out, size_t size) const {
    wrs::generators::stl<> gen(seed_);
    gen.generate_block(out, size, /* left_open */ true);
    double sum = 0.0;
    for (const double &w : out)
        sum += w;
    return sum;
}

// [1.1, 1-0.1/(n-1), ..., 1-0.1/(n-1)]
double one_heavy_gen::operator()(std::vector<double> &out, size_t size) const {
    out[0] = 1.1;
    double others = 1 - 0.1 / (size - 1);
    std::fill(out.begin() + 1, out.end(), others);
    return size;
}

// [0.9, 1+0.1/(n-1), ..., 1+0.1/(n-1)]
double one_light_gen::operator()(std::vector<double> &out, size_t size) const {
    out[0] = 0.5;
    double others = 1 + 0.5 / (size - 1);
    std::fill(out.begin() + 1, out.end(), others);
    return size;
}

// [n, 1, ..., 1]
double one_heavy2_gen::operator()(std::vector<double> &out, size_t size) const {
    out[0] = size;
    std::fill(out.begin() + 1, out.end(), 1);
    return 2 * size - 1;
}


// [1^-exp, 2^-exp, 3^-exp, ..., n^-exp]
double powerlaw_gen::operator()(std::vector<double> &out, size_t size) const {
    double sum = 0;
    for (size_t i = 0; i < size; i++) {
        out[i] = std::pow(static_cast<double>(i + 1), exp_);
        sum += out[i];
    }
    return sum;
}


// [1^-exp, 2^-exp, 3^-exp, ..., n^-exp], but shuffled
double powerlaw_shuffle_gen::operator()(std::vector<double> &out, size_t size) const {
    double sum = powerlaw_gen(exp_)(out, size);
    wrs::generators::stl<> generator(42);
    for (size_t i = size - 1; i != 0; i--) {
        size_t j = generator.next_int<size_t>(0, i);
        std::swap(out[i], out[j]);
    }
    return sum;
}
