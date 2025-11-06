#pragma once

#include "xmpz.hpp"

namespace xmp {

cuxmp_stat_t xmpz_div_limb(xmpz_t& dst, const xmpz_t& left_op,
                           const cuxmp_limb_t& right_op);

cuxmp_stat_t xmpz_div_limb_rem(xmpz_t& dst, const xmpz_t& left_op,
                               const cuxmp_limb_t& right_op,
                               cuxmp_limb_t& rem_dst);

// at least left_op.n must be reserved; left_op must not be 0; right_op must not
// be 0;
CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_div_limb_rem(
    xmpz_t& dst, const xmpz_t& left_op, const cuxmp_limb_t& right_op) {
    dst.n = left_op.n;

    cuxmp_sum_t rem = 0;
    for (cuxmp_len_t i = left_op.n; i > 0;) {
        i--;
        cuxmp_sum_t cur =
            (rem << CUXMP_LIMB_BITS) | (cuxmp_sum_t)left_op.limbs[i];
        dst.limbs[i] = static_cast<cuxmp_limb_t>(cur / right_op);
        rem = cur % right_op;
    }

    // remove trailing zeroes
    while (dst.n > 1 && dst.limbs[dst.n - 1] == 0) dst.n--;

    return static_cast<cuxmp_limb_t>(rem);
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t
_xmpz_div_eq_limb(xmpz_t& dst, const cuxmp_limb_t& right_op) {
    return _xmpz_div_limb_rem(dst, dst, right_op);
}

}  // namespace xmp
