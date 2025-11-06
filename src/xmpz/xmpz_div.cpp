#include "xmpz/xmpz_div.hpp"

#include "xmpz/xmpz.hpp"

using namespace xmp;

cuxmp_stat_t xmp::xmpz_div_limb(xmpz_t& dst, const xmpz_t& left_op,
                                const cuxmp_limb_t& right_op) {
    if (left_op.n == 0) {
        dst.reserve(1);
        dst.n = 1;
        dst.limbs[0] = 0;
        return 0;
    }

    if (right_op == 0) {
        return CUXMP_STAT_DIV_BY_ZERO;
    }

    dst.reserve(left_op.n);

    (void)_xmpz_div_limb_rem(dst, left_op, right_op);

    return CUXMP_STAT_OK;
}

cuxmp_stat_t xmp::xmpz_div_limb_rem(xmpz_t& dst, const xmpz_t& left_op,
                                    const cuxmp_limb_t& right_op,
                                    cuxmp_limb_t& rem_dst) {
    if (left_op.n == 0) {
        dst.reserve(1);
        dst.n = 1;
        dst.limbs[0] = 0;
        return 0;
    }

    if (right_op == 0) {
        return CUXMP_STAT_DIV_BY_ZERO;
    }

    dst.reserve(left_op.n);

    rem_dst = _xmpz_div_limb_rem(dst, left_op, right_op);

    return CUXMP_STAT_OK;
}
