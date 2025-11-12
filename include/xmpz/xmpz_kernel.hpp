#pragma once

#include <vector>

#include "xmpz.hpp"
#include "xmpz_kernel.hpp"

namespace xmp {

// kernels for NTT and CRT

struct xmp_ntt_prime {
    cuxmp_limb_t p;     // the prime
    cuxmp_limb_t root;  // root of unity
    cuxmp_sum_t m;      // the magic coefficient
};

constexpr static cuxmp_limb_t max_ntt_len = 1u << 24;  // 2^24

constexpr static xmp_ntt_prime ntt_primes[] = {
    {4194304001u, 250, 4398046510ull},
    {4076863489u, 243, 4524739208ull},
    {3942645761u, 235, 4678772883ull},
    {3892314113u, 232, 4739274256ull},
    {3489660929u, 208, 5286113594ull}};

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_ntt_add(cuxmp_limb_t a, cuxmp_limb_t b,
                                               cuxmp_limb_t p) {
    cuxmp_limb_t res = a + b;

    // i think this does (a + b) mod p but without an expensive division. but
    // limited to 2 * p.
    return (res >= p) ? (res - p) : res;
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_ntt_sub(cuxmp_limb_t a, cuxmp_limb_t b,
                                               cuxmp_limb_t p) {
    cuxmp_sum_t res = (cuxmp_sum_t)a - (cuxmp_sum_t)b;

    // not really sure what this does ??
    return (cuxmp_limb_t)(res + (p & -((cuxmp_limb_t)res > a)));
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_ntt_barret_mul(cuxmp_limb_t a,
                                                      cuxmp_limb_t b,
                                                      const xmp_ntt_prime& p) {
    cuxmp_sum_t prod = (cuxmp_sum_t)a * b;

    cuxmp_crt_coef_t p_prod = (cuxmp_crt_coef_t)p.m * prod;
    cuxmp_sum_t q_est = (cuxmp_sum_t)(p_prod >> CUXMP_SUM_BITS);

    cuxmp_sum_t r = prod - q_est * p.p;

    return (cuxmp_limb_t)((r >= p.p) ? (r - p.p) : r);
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_ntt_shoup_mul(cuxmp_limb_t a,
                                                     cuxmp_limb_t w_shoup,
                                                     cuxmp_limb_t w,
                                                     cuxmp_limb_t p) {
    cuxmp_limb_t q =
        (cuxmp_limb_t)(((cuxmp_sum_t)a * w_shoup) >> CUXMP_LIMB_BITS);

    int64_t r = (int64_t)a * w - (int64_t)q * p;

    return (cuxmp_limb_t)(r + (p & -(r >> 63)));
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _xmpz_ntt_power(cuxmp_limb_t base,
                                                 cuxmp_limb_t exp,
                                                 const xmp_ntt_prime& p) {
    cuxmp_limb_t res = 1;
    base = _xmpz_ntt_barret_mul(base, 1, p);

    while (exp > 0) {
        if (exp & 1) {
            res = _xmpz_ntt_barret_mul(res, base, p);
        }
        base = _xmpz_ntt_barret_mul(base, base, p);
        exp >>= 1;
    }

    return res;
}

CUXMP_ALWAYS_INLINE void _xmpz_ntt_bit_reverse(cuxmp_limb_t* limbs,
                                               cuxmp_len_t n) {
    for (cuxmp_len_t i = 1, j = 0; i < n; i++) {
        cuxmp_len_t bit = n >> 1;

        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }

        j ^= bit;

        if (i < j) {
            std::swap(limbs[i], limbs[j]);
        }
    }
}

// performs forward NTT
CUXMP_ALWAYS_INLINE void _ui32ntt_fwd_kernel(cuxmp_limb_t* limbs, cuxmp_len_t n,
                                             const xmp_ntt_prime& prime) {
    std::vector<cuxmp_limb_t> twiddles(n / 2);
    std::vector<cuxmp_limb_t> twiddles_shoup(n / 2);

    cuxmp_limb_t n_root = prime.root;
    cuxmp_len_t n_check = max_ntt_len;

    while (n_check > n) {
        n_root = _xmpz_ntt_barret_mul(n_root, n_root, prime);
        n_check >>= 1;
    }

    twiddles[0] = 1;

    twiddles_shoup[0] = (cuxmp_limb_t)(((__uint128_t)1 << 32) / prime.p);

    for (cuxmp_len_t i = 1; i < n / 2; i++) {
        twiddles[i] = _xmpz_ntt_barret_mul(twiddles[i - 1], n_root, prime);
        twiddles_shoup[i] =
            (cuxmp_limb_t)(((__uint128_t)twiddles[i] << 32) / prime.p);
    }

    // butterfly loops

    for (cuxmp_len_t len = n; len > 1; len >>= 1) {
        cuxmp_len_t half_len = len >> 1;

        cuxmp_len_t twiddle_step = (n / 2) / half_len;

        for (cuxmp_len_t i = 0; i < n; i += len) {
            cuxmp_limb_t* w = twiddles.data();
            cuxmp_limb_t* w_s = twiddles_shoup.data();

            for (cuxmp_len_t j = 0; j < half_len; j++) {
                cuxmp_limb_t u = limbs[i + j];
                cuxmp_limb_t v = limbs[i + j + half_len];

                limbs[i + j] = _xmpz_ntt_add(u, v, prime.p);

                cuxmp_limb_t t = _xmpz_ntt_sub(u, v, prime.p);

                limbs[i + j + half_len] =
                    _xmpz_ntt_shoup_mul(t, *w_s, *w, prime.p);

                w += twiddle_step;
                w_s += twiddle_step;
            }
        }
    }
}

// performs inverse NTT
CUXMP_ALWAYS_INLINE void _ui32ntt_inv_kernel(cuxmp_limb_t* limbs, cuxmp_len_t n,
                                             const xmp_ntt_prime& prime) {
    std::vector<cuxmp_limb_t> twiddles(n / 2);
    std::vector<cuxmp_limb_t> twiddles_shoup(n / 2);

    cuxmp_limb_t n_root = prime.root;
    cuxmp_len_t n_check = max_ntt_len;

    while (n_check > n) {
        n_root = _xmpz_ntt_barret_mul(n_root, n_root, prime);
        n_check >>= 1;
    }

    cuxmp_limb_t n_root_inv = _xmpz_ntt_power(n_root, prime.p - 2, prime);

    twiddles[0] = 1;

    twiddles_shoup[0] = (cuxmp_limb_t)(((__uint128_t)1 << 32) / prime.p);

    for (cuxmp_len_t i = 1; i < n / 2; i++) {
        twiddles[i] = _xmpz_ntt_barret_mul(twiddles[i - 1], n_root_inv, prime);
        twiddles_shoup[i] =
            (cuxmp_limb_t)(((__uint128_t)twiddles[i] << 32) / prime.p);
    }

    // butterly loops

    for (cuxmp_len_t len = 2; len <= n; len <<= 1) {
        cuxmp_len_t half_len = len >> 1;
        cuxmp_len_t twiddle_step = (n / 2) / half_len;

        for (cuxmp_len_t i = 0; i < n; i += len) {
            cuxmp_limb_t* w = twiddles.data();
            cuxmp_limb_t* w_s = twiddles_shoup.data();

            for (cuxmp_len_t j = 0; j < half_len; j++) {
                cuxmp_limb_t v = _xmpz_ntt_shoup_mul(limbs[i + j + half_len],
                                                     *w_s, *w, prime.p);
                cuxmp_limb_t u = limbs[i + j];

                limbs[i + j] = _xmpz_ntt_add(u, v, prime.p);
                limbs[i + j + half_len] = _xmpz_ntt_sub(u, v, prime.p);

                w += twiddle_step;
                w_s += twiddle_step;
            }
        }
    }

    cuxmp_limb_t n_inv = _xmpz_ntt_power(n, prime.p - 2, prime);

    for (cuxmp_len_t i = 0; i < n; i++) {
        limbs[i] = _xmpz_ntt_barret_mul(limbs[i], n_inv, prime);
    }
}

CUXMP_ALWAYS_INLINE cuxmp_crt_coef_t
_xmpz_crt_solve(const cuxmp_limb_t* prime_limbs, const xmp_ntt_prime* primes,
                const cuxmp_len_t n_primes) {
    // last thing left for nttcrt
}

}  // namespace xmp
