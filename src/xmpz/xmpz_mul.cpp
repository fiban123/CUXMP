#include "xmpz/xmpz_mul.hpp"

#include "xmpz/xmpz.hpp"

using namespace xmp;

cuxmp_stat_t xmp::xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                              const xmpz_t& right_op) {
    dst.reserve(left_op.n + right_op.n);

    (void)_xmpz_mul_sb(dst, left_op, right_op);

    return CUXMP_STAT_OK;
}

cuxmp_stat_t xmp::xmpz_mul_nttcrt(xmpz_t& dst, const xmpz_t left_op,
                                  const xmpz_t& right_op) {
    dst.reserve(left_op.n + right_op.n);

    cuxmp_len_t n_primes;
    cuxmp_limb_t* primes = _xmpz_ntt_find_primes(dst.reserved, n_primes);

    (void)_xmpz_mul_nttcrt(dst, left_op, right_op, primes, n_primes);

    free(primes);

    return CUXMP_STAT_OK;
}

cuxmp_stat_t xmp::xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                              const xmpz_t& right_op) {
    _xmpz_mul_ss(dst, left_op, right_op);

    return CUXMP_STAT_OK;
}

cuxmp_stat_t xmp::xmpz_mul(xmpz_t& dst, const xmpz_t& left_op,
                           const xmpz_t& right_op) {
    cuxmp_len_t n_result = left_op.n + right_op.n;

    if (n_result < CUXMP_MUL_SB_MAX) {
        return xmpz_mul_sb(dst, left_op, right_op);
    } else if (n_result < CUXMP_MUL_NTTCRT_MAX) {
        return xmpz_mul_nttcrt(dst, left_op, right_op);
    }
    return xmpz_mul_ss(dst, left_op, right_op);
}
