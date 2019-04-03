/*******************************************************************************
 * benchmark/rmat_graph500.hpp
 *
 * graph500 R-MAT generation stuff
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 *
 * Graph500 code:
 * Copyright (C) 2009-2010 The Trustees of Indiana University.                 *
 *                                                                             *
 * Use, modification and distribution is subject to the Boost Software         *
 * License, Version 1.0. (See http://www.boost.org/LICENSE_1_0.txt)            *
 *                                                                             *
 *  Authors: Jeremiah Willcock                                                 *
 *           Andrew Lumsdaine                                                  *
 ******************************************************************************/

#pragma once
#ifndef BENCHMARK_RMAT_GRAPH500_HEADER
#define BENCHMARK_RMAT_GRAPH500_HEADER


// disable parallelisation
#ifdef _OPENMP
#undef _OPENMP
#endif // _OPENMP

#ifdef __MTA__
#undef __MTA__
#endif // __MTA__

// generate unweighted graphs
#ifdef SSSP
#undef SSSP
#endif // SSSP

#define restrict __restrict__

#include <extlib/graph500/generator/splittable_mrg.c>
#include <extlib/graph500/generator/graph_generator.c>

#include <wrs/aggregate.hpp>

#include <tlx/define.hpp>

#include <cassert>
#include <limits>
#include <stdint.h>

// check graph500 user settings
#ifndef FAST_64BIT_ARITHMETIC
#warning "FAST_64BIT_ARITHMETIC is not defined, have you mucked with graph500's generator/user_settings.h?"
#endif


namespace graph500 {

struct scramble_state {
    uint64_t val0;
    uint64_t val1;
    int lgN;
    scramble_state() = delete;
    scramble_state(uint64_t v0, uint64_t v1, int lgN_) : val0(v0), val1(v1), lgN(lgN_) {}
};

/* adapted from generator/graph_generator.c generate_kronecker_range */
scramble_state init_scramble_helper(const uint_fast32_t seed[5], int lgN) {
    mrg_state state;
    mrg_seed(&state, seed);

    uint64_t val0, val1; /* Values for scrambling */
    {
      mrg_state new_state = state;
      mrg_skip(&new_state, 50, 7, 0);

      val0 = mrg_get_uint_orig(&new_state);
      val0 *= UINT64_C(0xFFFFFFFF);
      val0 += mrg_get_uint_orig(&new_state);
      val1 = mrg_get_uint_orig(&new_state);
      val1 *= UINT64_C(0xFFFFFFFF);
      val1 += mrg_get_uint_orig(&new_state);
    }

    return scramble_state(val0, val1, lgN);
}


// init scramble state, generating seed from RNG
template <typename RNG>
scramble_state init_scramble_state(RNG &rng, int lgN) {
    uint_fast32_t seed[5]; /* All values in [0, 2^31 - 1), not all zero */
    for (int i = 0; i < 5; i++) {
        seed[i] = static_cast<uint_fast32_t>(
            rng.next() * ((static_cast<uint_fast32_t>(1) << 31) - 1));
    }

    return init_scramble_helper(seed, lgN);
}


/* from graph500 generator/graph_generator.c */
/* Reverse bits in a number; this should be optimized for performance
 * (including using bit- or byte-reverse intrinsics if your platform has them).
 * */
constexpr uint64_t bitreverse(uint64_t x) {
  /* 64-bit code */
  x = __builtin_bswap64(x);
  x = ((x >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) |
      ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
  x = ((x >> 2) & UINT64_C(0x3333333333333333)) |
      ((x & UINT64_C(0x3333333333333333)) << 2);
  x = ((x >> 1) & UINT64_C(0x5555555555555555)) |
      ((x & UINT64_C(0x5555555555555555)) << 1);
  return x;
}

constexpr void bitreverse_two(uint64_t &x, uint64_t &y) {
  /* 64-bit code */
  x = __builtin_bswap64(x);
  y = __builtin_bswap64(y);
  x = ((x >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) |
      ((x & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
  y = ((y >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F)) |
      ((y & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4);
  x = ((x >> 2) & UINT64_C(0x3333333333333333)) |
      ((x & UINT64_C(0x3333333333333333)) << 2);
  y = ((y >> 2) & UINT64_C(0x3333333333333333)) |
      ((y & UINT64_C(0x3333333333333333)) << 2);
  x = ((x >> 1) & UINT64_C(0x5555555555555555)) |
      ((x & UINT64_C(0x5555555555555555)) << 1);
  y = ((y >> 1) & UINT64_C(0x5555555555555555)) |
      ((y & UINT64_C(0x5555555555555555)) << 1);
}


TLX_ATTRIBUTE_ALWAYS_INLINE
constexpr void scramble_two(int64_t &v0, int64_t &w0,
                            int lgN, uint64_t val0, uint64_t val1) {
    uint64_t v = (uint64_t)v0;
    uint64_t w = (uint64_t)w0;
    v += val0 + val1;
    w += val0 + val1;
    v *= (val0 | UINT64_C(0x4519840211493211));
    w *= (val0 | UINT64_C(0x4519840211493211));
    /*
    v = (bitreverse(v) >> (64 - lgN));
    w = (bitreverse(w) >> (64 - lgN));
    */
    bitreverse_two(v, w);
    v = (v >> (64 - lgN));
    w = (w >> (64 - lgN));

    assert ((v >> lgN) == 0);
    assert ((w >> lgN) == 0);
    v *= (val1 | UINT64_C(0x3050852102C843A5));
    w *= (val1 | UINT64_C(0x3050852102C843A5));
    v = (bitreverse(v) >> (64 - lgN));
    w = (bitreverse(w) >> (64 - lgN));
    assert ((v >> lgN) == 0);
    assert ((w >> lgN) == 0);
    v0 = (int64_t)v;
    w0 = (int64_t)w;
}

TLX_ATTRIBUTE_ALWAYS_INLINE
constexpr void scramble_two(int64_t &v0, int64_t &w0, const scramble_state &state) {
    scramble_two(v0, w0, state.lgN, state.val0, state.val1);
}


/* Make a single graph edge using a pre-set MRG state. */
template <bool scramble_ids, typename callback>
static
void make_one_edge(int64_t nverts, int level, int lgN, mrg_state* st,
                   callback && cb, uint64_t __attribute__((unused)) val0,
                   uint64_t __attribute__((unused)) val1) {
  int64_t base_src = 0, base_tgt = 0;
  while (nverts > 1) {
    int square = generate_4way_bernoulli(st, level, lgN);
    int src_offset = square / 2;
    int tgt_offset = square % 2;
    assert (base_src <= base_tgt);
#ifdef RMAT_CLIPFLIP
    if (base_src == base_tgt) {
      /* Clip-and-flip for undirected graph */
      if (src_offset > tgt_offset) {
        int temp = src_offset;
        src_offset = tgt_offset;
        tgt_offset = temp;
      }
    }
#endif // RMAT_CLIPFLIP
    nverts /= 2;
    ++level;
    base_src += nverts * src_offset;
    base_tgt += nverts * tgt_offset;
  }
  if constexpr(scramble_ids) {
      scramble_two(base_src, base_tgt, lgN, val0, val1);
  }
  cb(base_src, base_tgt);
}

} // namespace graph500

// Wrapper around the graph500 generator
template <bool Scramble_IDs = true>
class rmat_graph500 {
public:
    static constexpr const char* name = "graph500";

    using node = int64_t;
    static constexpr bool scramble_ids = Scramble_IDs;

    template <typename RNG>
    rmat_graph500(RNG &rng, int log_n_, double, double, double):
        n(static_cast<size_t>(1) << log_n_),
        log_n(log_n_)
    {
        uint_fast32_t seed[5]; /* All values in [0, 2^31 - 1), not all zero */
        for (int i = 0; i < 5; i++) {
            seed[i] = static_cast<uint_fast32_t>(
                rng.next() * ((static_cast<uint_fast32_t>(1) << 31) - 1));
        }
        mrg_seed(&state, seed);

        {
            mrg_state new_state = state;
            mrg_skip(&new_state, 50, 7, 0);
            val0 = mrg_get_uint_orig(&new_state);
            val0 *= UINT64_C(0xFFFFFFFF);
            val0 += mrg_get_uint_orig(&new_state);
            val1 = mrg_get_uint_orig(&new_state);
            val1 *= UINT64_C(0xFFFFFFFF);
            val1 += mrg_get_uint_orig(&new_state);
        }
    }

    template <typename RNG>
    std::pair<node, node> get_edge(size_t ei, RNG&) const {
        mrg_state new_state = state;
        mrg_skip(&new_state, 0, (uint64_t)ei, 0);

        node src, dst;
        auto callback = [&](auto src_, auto dst_) {
            src = src_;
            dst = dst_;
        };
        graph500::make_one_edge<scramble_ids>(n, 0, log_n, &new_state, callback, val0, val1);
        return std::make_pair(src, dst);
    }

    template <typename RNG, typename Callback>
    void get_edges(Callback && callback, size_t min, size_t max, RNG&) const {
        for (size_t ei = min; ei < max; ++ei) {
            mrg_state new_state = state;
            mrg_skip(&new_state, 0, (uint64_t)ei, 0);
            graph500::make_one_edge<scramble_ids>(
                n, 0, log_n, &new_state, callback, val0, val1);
        }
    }

    size_t table_size() const {
        return 0;
    }

    wrs::Aggregate<double> get_depth_stats() const {
        return wrs::Aggregate<double>().Add(log_n);
    }

    wrs::Aggregate<double> get_sample_stats() const {
        return wrs::Aggregate<double>().Add(1);
    }


protected:
    mrg_state state;
    size_t n;
    uint64_t val0, val1;
    int log_n;
};


#endif // BENCHMARK_RMAT_GRAPH500_HEADER
