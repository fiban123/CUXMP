#pragma once

#include <cmath>
#include <cstdlib>
#include <format>
#include <iostream>

#include "cuxmpdef.hpp"

namespace xmp {

struct xmpz_t {
    cuxmp_limb_t* limbs = nullptr;
    cuxmp_len_t n = 0;
    cuxmp_len_t reserved = 0;

    xmpz_t();
    void reserve(cuxmp_len_t _n);
};

// primitive operations

// schoolbook addition
cuxmp_stat_t xmpz_add_sb(xmpz_t& dst, const xmpz_t& bigger,
                         const xmpz_t& smaller);

// operations

// add
cuxmp_stat_t xmpz_add(xmpz_t& dst, const xmpz_t& left_op,
                      const xmpz_t& right_op);
// sub
// mul
// div
cuxmp_stat_t xmpz_div_limb(xmpz_t& dst, const xmpz_t& left_op,
                           const cuxmp_limb_t& right_op);

// util

// setters
void xmpz_set_u32(xmpz_t& dst, const uint32_t value);
void xmpz_set_u64(xmpz_t& dst, const uint64_t value);
void xmpz_set_limb(xmpz_t& dst, const cuxmp_limb_t& value);
void xmpz_set(xmpz_t& dst, const xmpz_t& value);

// getters
char* xmpz_get_str(const xmpz_t& src, cuxmp_base_t base, cuxmp_len_t& out_len);

// free
void xmpz_free(xmpz_t& value);

// helpers

// add with carry
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

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_div_limb(xmpz_t& dst,
                                                const xmpz_t& left_op,
                                                const cuxmp_limb_t& right_op) {
    if (left_op.n == 0) {
        dst.reserve(1);
        dst.n = 1;
        dst.limbs[0] = 0;
        return 0;
    }

    dst.reserve(left_op.n);
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
    return _xmpz_div_limb(dst, dst, right_op);
}

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
    xmpz_free(temp);
    return digits;
}

}  // namespace xmp
