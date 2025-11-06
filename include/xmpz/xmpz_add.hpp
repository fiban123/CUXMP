#pragma once

#include <format>
#include <iostream>

#include "xmpz.hpp"

namespace xmp {

// schoolbook addition
cuxmp_stat_t xmpz_add_sb(xmpz_t& dst, const xmpz_t& bigger,
                         const xmpz_t& smaller);

cuxmp_stat_t xmpz_add(xmpz_t& dst, const xmpz_t& left_op,
                      const xmpz_t& right_op);

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
