#pragma once

#include <cmath>
#include <cstdlib>

#include "xmpz.hpp"
#include "xmpz_div.hpp"
#include "xmpz_set.hpp"

namespace xmp {

char* xmpz_get_str(const xmpz_t& src, cuxmp_base_t base, cuxmp_len_t& out_len);

CUXMP_ALWAYS_INLINE char* _xmpz_get_str(const xmpz_t& src, cuxmp_base_t base,
                                        cuxmp_len_t& out_len) {
    // handle 0
    if (src.n == 0 || (src.n == 1 && src.limbs[0] == 0)) {
        char* str = (char*)malloc(sizeof(char));
        str[0] = '0';
        out_len = 1;
        return str;
    }

    xmpz_t temp;
    xmpz_set(temp, src);

    // calculate log_base(sum_max)
    double max_digits_per_limb =
        (double)std::log2(CUXMP_LIMB_MAX) / (double)std::log2(base);
    cuxmp_len_t max_digits =
        (cuxmp_sum_t)(std::ceil(max_digits_per_limb * src.n));

    char* digits = (char*)malloc(max_digits * sizeof(char));

    cuxmp_len_t len = 0;

    while (!(temp.n == 1 && temp.limbs[0] == 0)) {
        cuxmp_limb_t rem = _xmpz_div_eq_limb(temp, base);

        digits[len] = CUXMP_STR_CHARSET[rem % base];
        len++;
    }

    // reverse digits
    for (cuxmp_len_t i = 0; i < len / 2; i++) {
        char t = digits[i];
        digits[i] = digits[len - i - 1];
        digits[len - i - 1] = t;
    }

    out_len = len;
    temp.free();
    return digits;
}

}  // namespace xmp
