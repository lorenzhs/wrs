/*******************************************************************************
 * wrs/gsl.hpp
 *
 * Wrapper around the GSL alias method implementation
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#pragma once
#ifndef WRS_GSL_HEADER
#define WRS_GSL_HEADER

#ifdef WRS_HAVE_GSL

#include <wrs/timer.hpp>

#include <tlx/logger.hpp>

#include <gsl/gsl_randist.h>

#include <cassert>
#include <random>
#include <utility>

namespace wrs {
namespace gsl {

extern "C" {

typedef struct {
    size_t N;                      /* max number of elts on stack */
    size_t *v;                     /* array of values on the stack */
    size_t i;                      /* index of top of stack */
} gsl_stack_t;


// Custom gsl_ran_discrete_preproc without the allocations
    gsl_ran_discrete_t * my_gsl_ran_discrete_preproc(
        size_t K, const double *P, gsl_ran_discrete_t *g, double* E,
        gsl_stack_t *bigs, gsl_stack_t *smalls);

} // extern "C"

} // namespace gsl

// Wrapper struct
struct gsl_alias {
    static constexpr bool pass_rand = false;
    static constexpr bool debug = false;
    static constexpr bool time = false;

    using result_type = size_t;

private:
    gsl::gsl_stack_t* new_stack(size_t size) {
        gsl::gsl_stack_t *s = new gsl::gsl_stack_t();
        s->N = size;
        s->i = 0;                  /* indicates stack is empty */
        s->v = new size_t[size];
        return s;
    }

public:
    gsl_alias() : r(nullptr), g(nullptr), E(nullptr), size(-1) {}

    template <typename Iterator>
    gsl_alias(Iterator begin, Iterator end) {
        init(end - begin);

        construct(begin, end);
    }

    void init(size_t size_, size_t seed = 0) {
        size = size_;
        r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, seed != 0 ? seed : std::random_device{}());
        g = new gsl_ran_discrete_t();
        g->K = size;
        g->F = new double[size];
        g->A = new size_t[size];
        E = new double[size];

        bigs = new_stack(size);
        smalls = new_stack(size);
    }

    ~gsl_alias() {
        assert(r != nullptr && g != nullptr && E != nullptr);
        assert(bigs != nullptr && smalls != nullptr);
        gsl_rng_free(r);
        delete[] g->F;
        delete[] g->A;
        delete g;
        delete[] E;
        delete[] bigs->v;
        delete bigs;
        delete[] smalls->v;
        delete smalls;
    }

    template <typename Iterator>
    void construct(Iterator begin, Iterator end) {
        if (end - begin != static_cast<ssize_t>(size)) {
            sLOG1 << "Error: tried to construct alias table of incorrect size!"
                  << "Expected" << size << "got" << end - begin;
            return;
        }
        timers_.clear();
        timer t;

        bigs->i = 0;
        smalls->i = 0;
        gsl::my_gsl_ran_discrete_preproc(
            size, static_cast<double*>(&(*begin)), g, E, bigs, smalls);

        timers_.push_back(t.get_and_reset());
        LOGC(time) << "GSL construction took " << timers_.back() << "ms";
    }

    size_t sample() {
        return gsl_ran_discrete(r, g);
    }

    std::vector<double> get_timers() const {
        return timers_; // copy
    }

    template <typename Iterator>
    void verify(Iterator, Iterator) {}

    std::vector<double> timers_; // for measurements
    gsl_rng *r;
    gsl_ran_discrete_t *g;
    gsl::gsl_stack_t *bigs, *smalls;
    double *E;
    size_t size;
};

} // namespace wrs

#endif // WRS_HAVE_GSL

#endif // WRS_GSL_HEADER
