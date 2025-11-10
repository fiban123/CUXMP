#pragma once

#include <cassert>
#include <format>
#include <iostream>

#include "xmpz.hpp"

namespace xmp {

// schoolbook addition
cuxmp_stat_t xmpz_add_sb(xmpz_t& dst, const xmpz_t& bigger,
                         const xmpz_t& smaller);

cuxmp_stat_t xmpz_add(xmpz_t& dst, const xmpz_t& left_op,
                      const xmpz_t& right_op);

// !: dst.n >= offset + 6
CUXMP_ALWAYS_INLINE void _xmpz_coef_add_eq_offset(xmpz_t& dst,
                                                  const cuxmp_crt_coef_t& coef,
                                                  const cuxmp_len_t offset) {
    assert(dst.n >= offset + 6);

    if (coef == 0) return;

    // these loops will hopefully be unrolled

    // add the uint128 to the limbs of dst
    cuxmp_sum_t carry = 0;
    for (cuxmp_static_bit_len_t i = 0;
         i < sizeof(cuxmp_crt_coef_t) / sizeof(cuxmp_limb_t); i++) {
        cuxmp_static_bit_len_t bit_offset = i * CUXMP_LIMB_BITS;
        cuxmp_sum_t sum = (cuxmp_sum_t)dst.limbs[offset] +
                          (cuxmp_limb_t)(coef >> bit_offset) + carry;

        dst.limbs[offset + i] = (cuxmp_limb_t)sum;
        carry = (sum >> CUXMP_LIMB_BITS);
    }

    // add any remaining carry
    if (carry != 0) {
        cuxmp_len_t carry_offset =
            offset + sizeof(cuxmp_crt_coef_t) / sizeof(cuxmp_limb_t);
        for (cuxmp_len_t i = 0; i < sizeof(cuxmp_sum_t) / sizeof(cuxmp_limb_t);
             i++) {
            cuxmp_sum_t sum = (cuxmp_sum_t)dst.limbs[carry_offset + i] + carry;
            dst.limbs[carry_offset + i] = (cuxmp_limb_t)sum;
            carry = sum >> CUXMP_LIMB_BITS;
            if (carry == 0) {
                break;
            }
        }
    }
}

cuxmp_limb_t CUXMP_ALWAYS_INLINE _xmp_add_with_carry(cuxmp_limb_t a,
                                                     cuxmp_limb_t b,
                                                     cuxmp_limb_t& carry) {
    cuxmp_sum_t sum = (cuxmp_sum_t)a + (cuxmp_sum_t)b + (cuxmp_sum_t)carry;

    cuxmp_limb_t result = (cuxmp_limb_t)sum;

    std::cout << std::format("a: {}, b: {}, carry: {}, sum: {}, result: {}", a,
                             b, carry, sum, result);
    carry = (cuxmp_limb_t)(sum >> CUXMP_LIMB_BITS);

    std::cout << "result carry: " << carry << "\n";

    return result;
}

}  // namespace xmp
