#pragma once

#include <cstdint>
#include <iostream>
#include <vector>

#include "xmpz.hpp"
#include "xmpz_kernel.hpp"

namespace xmp {

// kernels for NTT and CRT

struct xmp_ntt_prime {
    cuxmp_limb_t p;     // the prime
    cuxmp_limb_t root;  // root of unity
    cuxmp_sum_t m;      // the magic coefficient (2^64 / p)
};

constexpr static cuxmp_limb_t max_ntt_len = 1u << 23;  // 2^24

// biggest primes that fit in 31 bits such that k * 2^23 + 1 = p
constexpr static xmp_ntt_prime ntt_primes[] = {
    {2130706433u, 1791270792, 8657571868ull},    // 250 * 2^24 + 1 = p
    {2113929217u, 1722264568, 8726282755ull},    // 243 * 2^24 + 1 = p
    {2088763393u, 256941453, 8831418692ull},     // 235 * 2^24 + 1 = p
    {2013265921u, 1003846038, 9162596893ull},    // 232 * 2^24 + 1 = p
    {1811939329u, 1762019879, 10180663214ull}};  // 208 * 2^24 + 1 = p

CUXMP_ALWAYS_INLINE cuxmp_limb_t _ntt_addmod(cuxmp_limb_t a, cuxmp_limb_t b,
                                             cuxmp_limb_t p) {
    cuxmp_sum_t res = (cuxmp_sum_t)a + b;

    // i think this does (a + b) mod p but without an expensive division. but
    // limited to 2 * p.
    return (res >= p) ? (res - p) : res;
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _ntt_submod(cuxmp_limb_t a, cuxmp_limb_t b,
                                             cuxmp_limb_t p) {
    uint32_t res = a - b;
    // mask = 0xFFFFFFFF if a < b else 0
    uint32_t mask = (uint32_t)-(int32_t)(a < b);
    return res + (p & mask);
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _ntt_barret_mulmod(cuxmp_limb_t a,
                                                    cuxmp_limb_t b,
                                                    const xmp_ntt_prime& p) {
    cuxmp_sum_t prod = (cuxmp_sum_t)a * b;

    cuxmp_crt_coef_t p_prod = (cuxmp_crt_coef_t)p.m * prod;
    cuxmp_sum_t q_est = (cuxmp_sum_t)(p_prod >> CUXMP_SUM_BITS);

    cuxmp_sum_t r = prod - q_est * p.p;

    return (cuxmp_limb_t)((r >= p.p) ? (r - p.p) : r);
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _ntt_shoup_mulmod(cuxmp_limb_t a,
                                                   cuxmp_limb_t w_shoup,
                                                   cuxmp_limb_t w,
                                                   cuxmp_limb_t p) {
    cuxmp_limb_t q =
        (cuxmp_limb_t)(((cuxmp_sum_t)a * w_shoup) >> CUXMP_LIMB_BITS);

    int64_t r = (int64_t)a * (int64_t)w - (int64_t)q * (int64_t)p;
    if (r < 0)
        r += p;
    else if (r >= p)
        r -= p;

    return (cuxmp_limb_t)(r + (p & -(r >> 63)));
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _ntt_powmod(cuxmp_limb_t base,
                                             cuxmp_limb_t exp,
                                             const xmp_ntt_prime& p) {
    cuxmp_limb_t res = 1;
    base = _ntt_barret_mulmod(base, 1, p);

    while (exp > 0) {
        if (exp & 1) {
            res = _ntt_barret_mulmod(res, base, p);
        }
        base = _ntt_barret_mulmod(base, base, p);
        exp >>= 1;
    }

    return res;
}

CUXMP_ALWAYS_INLINE void _ntt_bit_reverse(cuxmp_limb_t* limbs, cuxmp_len_t n) {
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
        n_root = _ntt_barret_mulmod(n_root, n_root, prime);
        n_check >>= 1;
    }

    twiddles[0] = 1;

    twiddles_shoup[0] = (cuxmp_limb_t)(((__uint128_t)1 << 32) / prime.p);

    for (cuxmp_len_t i = 1; i < n / 2; i++) {
        twiddles[i] = _ntt_barret_mulmod(twiddles[i - 1], n_root, prime);
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

                limbs[i + j] = _ntt_addmod(u, v, prime.p);

                cuxmp_limb_t t = _ntt_submod(u, v, prime.p);

                limbs[i + j + half_len] =
                    _ntt_shoup_mulmod(t, *w_s, *w, prime.p);

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
        n_root = _ntt_barret_mulmod(n_root, n_root, prime);
        n_check >>= 1;
    }

    cuxmp_limb_t n_root_inv = _ntt_powmod(n_root, prime.p - 2, prime);

    twiddles[0] = 1;

    twiddles_shoup[0] = (cuxmp_limb_t)(((__uint128_t)1 << 32) / prime.p);

    for (cuxmp_len_t i = 1; i < n / 2; i++) {
        twiddles[i] = _ntt_barret_mulmod(twiddles[i - 1], n_root_inv, prime);
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
                cuxmp_limb_t v = _ntt_shoup_mulmod(limbs[i + j + half_len],
                                                   *w_s, *w, prime.p);
                cuxmp_limb_t u = limbs[i + j];

                limbs[i + j] = _ntt_addmod(u, v, prime.p);
                limbs[i + j + half_len] = _ntt_submod(u, v, prime.p);

                w += twiddle_step;
                w_s += twiddle_step;
            }
        }
    }

    cuxmp_limb_t n_inv = _ntt_powmod(n, prime.p - 2, prime);

    for (cuxmp_len_t i = 0; i < n; i++) {
        limbs[i] = _ntt_barret_mulmod(limbs[i], n_inv, prime);
    }
}

// CRT
inline cuxmp_limb_t _crt_const_inv[3][3] = {
    {0u, 2113929091u, 1253257986u}, {0u, 0u, 2088763310u}, {0u, 0u, 0u}};
inline cuxmp_crt_coef_t _crt_const_p_prod[3] = {1ull, 2130706433ull,
                                                4504162581568552961ull};

CUXMP_ALWAYS_INLINE cuxmp_limb_t _crt_modpow(cuxmp_limb_t base,
                                             cuxmp_limb_t exp,
                                             const xmp_ntt_prime& p) {
    __uint128_t res = 1;
    base %= p.p;

    while (exp > 0) {
        if (exp & 1) res = (res * base) % p.p;
        base = ((__uint128_t)base * base) % p.p;
        exp >>= 1;
    }

    return (cuxmp_limb_t)res;
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _crt_modsub(const cuxmp_limb_t a,
                                             const cuxmp_limb_t b,
                                             const xmp_ntt_prime& p) {
    int64_t res = (int64_t)a - (int64_t)b;

    int64_t mask = -(res < 0);  // -1 if res<0, else 0
    res += p.p & mask;          // add p if negative
    mask = -(res < 0);          // second check
    res += p.p & mask;

    return (cuxmp_limb_t)res;
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t _crt_modmul(const int64_t a,
                                             const cuxmp_limb_t b,
                                             const xmp_ntt_prime& p) {
    return (cuxmp_limb_t)(((cuxmp_sum_t)a * b) % p.p);
}

CUXMP_ALWAYS_INLINE cuxmp_crt_coef_t
_crt_solve_kernel(const cuxmp_limb_t* prime_limbs, const xmp_ntt_prime* primes,
                  const cuxmp_len_t n_primes) {
    cuxmp_limb_t v[3];

    v[0] = prime_limbs[0];

    // find v[i] for i > 0
    for (cuxmp_len_t i = 1; i < n_primes; i++) {
        v[i] = prime_limbs[i];

        for (cuxmp_len_t j = 0; j < i; j++) {
            cuxmp_limb_t tmp = _crt_modsub(v[i], v[j], primes[i]);

            v[i] = _crt_modmul(tmp, _crt_const_inv[j][i], primes[i]);

            // cuxmp_limb_t tmp = _crt_modsub(v[i], v[j], primes[i]);

            // v[i] = _crt_modmul(tmp, _crt_const_inv[j][i], primes[i]);
        }
    }

    cuxmp_crt_coef_t c = 0;

    for (cuxmp_len_t i = 0; i < n_primes; i++) {
        c = c + (cuxmp_crt_coef_t)v[i] * _crt_const_p_prod[i];
    }
    return c;
}

}  // namespace xmp
