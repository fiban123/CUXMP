#include "xmp.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>

using namespace xmp;

xmpz_t::xmpz_t() : limbs(nullptr), n(0) {}

void xmpz_t::reserve(cuxmp_len_t _n) {
    if (_n > reserved) {
        reserved = std::max(reserved * 2, _n);
        limbs = (cuxmp_limb_t*)realloc(limbs, reserved * sizeof(cuxmp_limb_t));
    }
}

// schoolbook addition
cuxmp_stat_t xmp::xmpz_add_sb(xmpz_t& dst, const xmpz_t& bigger,
                              const xmpz_t& smaller) {
    cuxmp_len_t i;
    cuxmp_limb_t carry = 0;

    assert(dst.reserved > bigger.n);
    assert(dst.limbs && bigger.limbs && smaller.limbs);

    dst.n = bigger.n;

    for (i = 0; i < smaller.n; i++) {
        dst.limbs[i] =
            _xmp_add_with_carry(smaller.limbs[i], bigger.limbs[i], carry);
    }

    for (; i < bigger.n; i++) {
        dst.limbs[i] = _xmp_add_with_carry(bigger.limbs[i], 0, carry);
    }

    if (carry > 0) {
        std::cout << "carry performed\n";
        dst.n++;
        dst.limbs[i] = carry;
    }

    return CUXMP_STAT_OK;
}

// operations

cuxmp_stat_t xmp::xmpz_add(xmpz_t& dst, const xmpz_t& left_op,
                           const xmpz_t& right_op) {
    const bool left_bigger = left_op.n > right_op.n;
    const xmpz_t& bigger = left_bigger ? left_op : right_op;
    const xmpz_t& smaller = left_bigger ? right_op : left_op;

    dst.reserve(bigger.n + 1);

    if (bigger.n > CUXMP_ADD_SB_MIN) {
        return xmpz_add_sb(dst, bigger, smaller);
    }

    // ...
    return 0;
}

// setters
void xmp::xmpz_set_u32(xmpz_t& dst, const uint32_t value) {
    if constexpr (sizeof(uint32_t) == sizeof(cuxmp_limb_t)) {
        xmpz_set_limb(dst, (cuxmp_limb_t)value);
    } else {
        _xmpz_init_ui<uint32_t>(dst, value);
    }
}

void xmp::xmpz_set_u64(xmpz_t& dst, const uint64_t value) {
    if constexpr (sizeof(uint64_t) == sizeof(cuxmp_limb_t) * 2) {
        dst.reserve(2);
        dst.limbs[0] = (cuxmp_limb_t)(value & CUXMP_LIMB_MAX);
        dst.limbs[1] =
            (cuxmp_limb_t)((value >> CUXMP_LIMB_BITS) & CUXMP_LIMB_MAX);
        dst.n = dst.limbs[1] == 0 ? 1 : 2;
    } else {
        _xmpz_init_ui<uint64_t>(dst, value);
    }
}

void xmp::xmpz_set_limb(xmpz_t& dst, const cuxmp_limb_t& value) {
    dst.reserve(1);
    dst.limbs[0] = value;
    dst.n = 1;
}

void xmp::xmpz_set(xmpz_t& result, const xmpz_t& value) {
    result.reserve(value.n);
    memcpy(result.limbs, value.limbs, value.n * sizeof(cuxmp_limb_t));

    result.n = value.n;
}

// getters
char* xmp::xmpz_get_str(const xmpz_t& src, cuxmp_base_t base,
                        cuxmp_len_t& out_len) {
    return _xmpz_get_str(src, base, out_len);
}

// misc

void xmp::xmpz_free(xmpz_t& value) { free(value.limbs); }
