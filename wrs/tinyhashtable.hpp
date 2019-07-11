/*******************************************************************************
 * wrs/tinyhashtable.hpp
 *
 * A tiny hash table
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#pragma once
#ifndef WRS_TINYHASHTABLE_HEADER
#define WRS_TINYHASHTABLE_HEADER

#include <tlx/math/round_to_power_of_two.hpp>

#include <memory>

namespace wrs {

// Hash table with no safety. Insert at most as many elements as you promised
// the constructor, or it will segfault. Clear before use or it will segfault.
// Power-of-two size recommended.
template <typename key_t, typename val_t>
class tinyhashtable {
public:
    static constexpr key_t empty = static_cast<key_t>(-1);
    struct entry {
        key_t key;
        val_t val;
        entry() : key(empty), val(0) {}
        entry(const key_t &k, const val_t &v) : key(k), val(v) {}
        bool operator==(const entry &other) const { return key == other.key; }
    };
    friend std::ostream &operator << (std::ostream &os, const entry &e) {
        return os << '(' << e.key << ',' << e.val << ')';
    }

    // need to call clear before use!!!
    explicit tinyhashtable(size_t max_size)
        : size_(4 * tlx::round_up_to_power_of_two(max_size))
        , mask_((size_ >> 1) - 1)
        , table_(max_size == 0 ? nullptr : std::make_unique<entry[]>(size_))
    {}

    val_t& operator[](const key_t &key) {
        size_t pos = key & mask_;
        entry *e  = table_.get() + pos;
        while (e->key != empty && e->key != key) {
            e++;
            assert(e - table_.get() < static_cast<ssize_t>(size_));
        }
        e->key = key; // doesn't matter if already set
        return e->val;
    }

    void clear() {
        const entry empty_e(empty, 0);
        std::fill(table_.get(), table_.get() + size_, empty_e);
    }

    template <typename Callback>
    void foreach(Callback && callback) const {
        entry* it = table_.get();
        const entry* end = table_.get() + size_;
        while (it != end) {
            if (it->key != empty)
                callback(it->key, it->val);
            it++;
        }
    }

protected:
    size_t size_, mask_;
    std::unique_ptr<entry[]> table_;
};

} // namespace wrs

#endif // WRS_TINYHASHTABLE_HEADER
