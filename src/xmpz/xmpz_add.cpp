#include "xmpz/xmpz_add.hpp"

#include <cassert>

using namespace xmp;

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
