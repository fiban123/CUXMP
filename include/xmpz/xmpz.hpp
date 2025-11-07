#pragma once

#include <algorithm>
#include <cstdlib>

#include "cuxmpdef.hpp"

namespace xmp {

struct xmpz_t {
    cuxmp_limb_t* limbs = nullptr;
    cuxmp_len_t n = 0;
    cuxmp_len_t reserved = 0;

    xmpz_t() : limbs(nullptr), n(0), reserved(0) {}

    CUXMP_ALWAYS_INLINE void reserve(cuxmp_len_t _n) {
        if (_n > reserved) {
            reserved = std::max(reserved * 2, _n);
            limbs = (cuxmp_limb_t*)std::realloc(
                limbs, reserved * sizeof(cuxmp_limb_t));
        }
    }
    CUXMP_ALWAYS_INLINE void free() { std::free(limbs); }
};

using xmpz_arr_t = xmpz_t*;

}  // namespace xmp
