#pragma once

#include <cassert>
#include <cstring>

#include "xmpz.hpp"

namespace xmp {

cuxmp_stat_t xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                         const xmpz_t& right_op);

// left.n + right.n must be reserved in dst
CUXMP_ALWAYS_INLINE void _xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                                      const xmpz_t& right_op) {
    assert(dst.reserved >= left_op.n + right_op.n);

    dst.n = left_op.n + right_op.n;

    // zero all needed limbs of dst
    std::memset(dst.limbs, 0, (left_op.n + right_op.n) * sizeof(cuxmp_limb_t));

    // schoolbook multiply
    for (cuxmp_len_t i = 0; i < left_op.n; i++) {
        cuxmp_sum_t carry = 0;
        cuxmp_limb_t left_limb = left_op.limbs[i];

        for (cuxmp_len_t j = 0; j < right_op.n; j++) {
            cuxmp_len_t index = i + j;
            cuxmp_limb_t right_limb = right_op.limbs[j];
            cuxmp_sum_t acc = (cuxmp_sum_t)dst.limbs[index] +
                              (cuxmp_sum_t)left_limb * (cuxmp_sum_t)right_limb +
                              carry;
            dst.limbs[index] = (cuxmp_limb_t)acc;
            carry = (cuxmp_sum_t)(acc >> CUXMP_LIMB_BITS);
        }

        // final carry
        dst.limbs[i + right_op.n] = (cuxmp_limb_t)carry;
    }

    while (dst.n > 1 && dst.limbs[dst.n - 1] == 0) dst.n--;
}

}  // namespace xmp
