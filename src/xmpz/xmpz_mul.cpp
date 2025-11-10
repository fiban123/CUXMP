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
    // find common N
    cuxmp_len_t n = 1;
    while (n < left_op.n + right_op.n - 1) {
        n <<= 1;
    }

    cuxmp_len_t n_primes;
    xmp_ntt_prime* primes = _xmpz_ntt_find_primes(n, n_primes);

    xmpz_t left_padded, right_padded;
    _xmpz_pad(left_padded, left_op, n);
    _xmpz_pad(right_padded, right_op, n);

    // pad dst and zero it
    dst.reserve(n);
    dst.n = n;
    memset(dst.limbs, 0, n * sizeof(cuxmp_limb_t));

    (void)_xmpz_mul_nttcrt(dst, left_padded, right_padded, primes, n_primes);

    free(primes);

    // remove trailing zeroes here

    return CUXMP_STAT_OK;
}

cuxmp_stat_t xmp::xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                              const xmpz_t& right_op) {
    dst.reserve(left_op.n + right_op.n + 1);

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
