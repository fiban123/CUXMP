#include "xmpz/xmpz_mul.hpp"

#include "xmpz/xmpz.hpp"

using namespace xmp;

cuxmp_stat_t xmp::xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                              const xmpz_t& right_op) {
    dst.reserve(left_op.n + right_op.n);

    (void)_xmpz_mul_sb(dst, left_op, right_op);

    return CUXMP_STAT_OK;
}
