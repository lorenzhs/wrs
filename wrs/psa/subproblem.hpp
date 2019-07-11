/*******************************************************************************
 * wrs/psa/subproblem.hpp
 *
 * Helper classes for parallel alias table construction
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_PSA_SUBPROBLEM_HEADER
#define WRS_PSA_SUBPROBLEM_HEADER

#include <cassert>
#include <ostream>

namespace wrs {
namespace psa {

template <typename size_type>
struct splitter {
    // first light item of the right side.  Light items are never split.
    size_type l_end;
    // heavy item that is split between both sides, or the first of the right
    // side if no item is split
    size_type h_split;
    // partial weight of the split item, need not add up to its weight if the
    // item is split over more than two subproblems
    double split_weight_left, split_weight_right;
    bool has_split_item, right_owns_split;

    friend std::ostream &operator << (std::ostream &os, const splitter &s) {
        return os << "splitter(l_end=" << s.l_end
                  << " h_split=" << s.h_split
                  << " split=(?" << s.has_split_item
                  << " O" << s.right_owns_split
                  << " lw=" << s.split_weight_left
                  << " rw=" << s.split_weight_right << "))";
    }
};

template <typename size_type>
struct subproblem {
    // light item range, exclusive to this subproblem.
    size_type l_begin, l_end;
    // heavy item range, h_begin and h_end-1 may be shared
    size_type h_begin, h_end;
    // number of *owned* items per group
    size_type size, num_light, num_heavy;
    // Shared items (partial weight and ownership flag for left side)
    // If flag is set, no right boundary item may be present
    double left_boundary_weight, right_boundary_weight;
    bool has_left_boundary, has_right_boundary;

    bool sanity_check() const {
        if ((l_end - l_begin) + (h_end - h_begin) - has_right_boundary != size) {
            sLOG1 << "size sanity check failed:" << l_end - l_begin
                  << "+" << h_end - h_begin << "="
                  << (l_end - l_begin) + (h_end - h_begin) << "!=" << size;
            assert((l_end - l_begin) + (h_end - h_begin) - has_right_boundary == size);
            return false;
        }
        if (has_left_boundary && has_right_boundary)
        {
            assert(!has_left_boundary || !has_right_boundary);
            return false;
        }

        return true;
    }

    friend std::ostream &operator << (std::ostream &os, const subproblem &p) {
        return os << "subproblem(l=" << p.l_begin << ".." << p.l_end
                  << "#" << p.num_light
                  << " h=" << p.h_begin << ".." << p.h_end
                  << "#" << p.num_heavy
                  << " #=" << p.size
                  << " lb=(?" << p.has_left_boundary
                  << " W" << p.left_boundary_weight
                  << ") rb=(?" << p.has_right_boundary
                  << " W" << p.right_boundary_weight << "))";
    }

};

} // namespace psa
} // namespace wrs

#endif // WRS_PSA_SUBPROBLEM_HEADER
