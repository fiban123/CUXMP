#pragma once

#include "xmpz.hpp"

namespace xmp {

void xmpz_set_u32(xmpz_t& dst, const uint32_t value);
void xmpz_set_u64(xmpz_t& dst, const uint64_t value);
void xmpz_set_limb(xmpz_t& dst, const cuxmp_limb_t& value);
void xmpz_set(xmpz_t& dst, const xmpz_t& src);
void xmpz_copy(xmpz_t& dst, const xmpz_t& src);

template <typename ui_t>
CUXMP_ALWAYS_INLINE void _xmpz_init_ui(xmpz_t& dst, const ui_t& value) {
    constexpr size_t bits = sizeof(ui_t) * 8;
    constexpr size_t limb_bits = CUXMP_LIMB_BITS;
    constexpr size_t limbs_needed = (bits + limb_bits - 1) / limb_bits;

    dst.reserve(limbs_needed);

    cuxmp_len_t i = 0;
    while (value && i < limbs_needed) {
        dst.limbs[i++] = (cuxmp_limb_t)(value & CUXMP_LIMB_MAX);
        value >>= limb_bits;
    }

    // handle zero value (still needs 1 limb)
    if (i == 0) {
        dst.limbs[0] = 0;
        i = 1;
    }

    dst.n = i;
}

}  // namespace xmp
