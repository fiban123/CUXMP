#pragma once

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

#include "xmpz.hpp"
#include "xmpz_add.hpp"
#include "xmpz_set.hpp"

namespace xmp {

// TEST SS and NTTC-CRT implementation
// performance is not important rn. just want a working SS and NTTCRT to
// optimize later.

// helpers for debugging

// boilerplate i copied online
inline std::string to_string_u128(__uint128_t value) {
    if (value == 0) return "0";
    std::string result;
    while (value > 0) {
        int digit = value % 10;
        result.push_back('0' + digit);
        value /= 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

inline void _d_xmpz_print_limb_arr(const cuxmp_limb_t* limbs, cuxmp_len_t n) {
    std::cout << "N = " << n << ": ";
    for (cuxmp_len_t i = 0; i < n; i++) {
        std::cout << limbs[i] << ", ";
    }
    std::cout << '\n';
}

inline void _d_xmpz_print_limbs(const xmpz_t& src) {
    _d_xmpz_print_limb_arr(src.limbs, src.n);
}

inline void _d_xmpz_print_xmpz_arr(const xmpz_arr_t src, const cuxmp_len_t n) {
    std::cout << "arr N = " << n << "\n";
    for (cuxmp_len_t i = 0; i < n; i++) {
        _d_xmpz_print_limbs(src[i]);
    }
}

inline void _d_xmpz_print_coef_arr(const cuxmp_crt_coef_t* arr,
                                   const cuxmp_len_t n) {
    std::cout << "coef N = " << n << ": ";

    for (cuxmp_len_t i = 0; i < n; i++) {
        std::cout << to_string_u128(arr[i]) << ", ";
    }

    std::cout << "\n";
}

// public
cuxmp_stat_t xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                         const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul_nttcrt(xmpz_t& dst, const xmpz_t left_op,
                             const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                         const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul(xmpz_t& dst, const xmpz_t& left_op,
                      const xmpz_t& right_op);

// helpers
// pad a number with 0s to length n
CUXMP_ALWAYS_INLINE void _xmpz_pad(xmpz_t& dst, const xmpz_t& src,
                                   const cuxmp_len_t n) {
    // finally something easy
    if (dst.n >= n) {
        return;
    }
    dst.reserve(n);
    xmpz_set(dst, src);
    memset(&dst.limbs[dst.n], 0, (n - dst.n) * sizeof(cuxmp_limb_t));
    dst.n = n;
}

CUXMP_ALWAYS_INLINE void _xmpz_pad_eq(xmpz_t& dst, cuxmp_len_t n) {
    if (dst.n >= n) {
        return;
    }
    dst.reserve(n);
    memset(&dst.limbs[dst.n], 0, (n - dst.n) * sizeof(cuxmp_limb_t));
    dst.n = n;
}

// NTT-CRT

CUXMP_ALWAYS_INLINE void _xmpz_ui32ntt(xmpz_t& dst, const xmpz_t& src,
                                       const cuxmp_limb_t prime) {
    // this and inv NTT probably most difficult and crucial for performance
}

CUXMP_ALWAYS_INLINE void _xmpz_ui32invntt(xmpz_t& dst, const xmpz_t& src,
                                          const cuxmp_limb_t prime) {
    // ...
}

CUXMP_ALWAYS_INLINE cuxmp_crt_coef_t
_xmpz_crt_solve(const cuxmp_limb_t* prime_limbs, const cuxmp_limb_t* primes,
                const cuxmp_len_t n_primes) {
    // does not sound easy
}

CUXMP_ALWAYS_INLINE void _xmpz_construct_coef_crt(xmpz_t& dst,
                                                  const xmpz_arr_t dst_parts,
                                                  const cuxmp_len_t dst_part_n,
                                                  const cuxmp_limb_t* primes,
                                                  const cuxmp_len_t n_primes) {
    // cuxmp_crt_coef_t is a __uint128_t on the CPU
    // on the GPU, this will be a custom uint128 made of 2 uint64s with add,
    // sub, ... functionality. way faster than a dynamically sized xmpz

    // precompute CRT constants
    // ...

    std::vector<cuxmp_limb_t> prime_limbs(n_primes);

    // for each coefficient
    for (cuxmp_len_t c_i = 0; c_i < dst_part_n; c_i++) {
        // for each prime
        for (cuxmp_len_t p_i = 0; p_i < n_primes; p_i++) {
            prime_limbs[p_i] = dst_parts[p_i].limbs[c_i];
        }

        // get true coeffficient
        cuxmp_crt_coef_t coef = _xmpz_crt_solve(prime_limbs.data(), primes,
                                                n_primes /* CRT constants */);

        // add coefficient
        _xmpz_coef_add_eq_offset(dst, coef, c_i);
    }
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t* _xmpz_ntt_find_primes(
    const cuxmp_len_t n, cuxmp_len_t& out_prime_n) {
    // just some if statements that return different predetermined primes
}

CUXMP_ALWAYS_INLINE void _xmpz_pointwise_mul_mod(xmpz_t& dst,
                                                 const xmpz_t& left_op,
                                                 const xmpz_t& right_op,
                                                 const cuxmp_limb_t m) {
    assert(left_op.n == right_op.n && left_op.n <= dst.n);

    // do dst = left * right mod m for each limb.
    // not sure if this is correct
    for (cuxmp_len_t i = 0; i < left_op.n; i++) {
        dst.limbs[i] = (cuxmp_sum_t)(left_op.limbs[i] * right_op.limbs[i]) %
                       (cuxmp_limb_t)m;
    }
}

// right_op.n must = left_op.n, N must be a power of 2
CUXMP_ALWAYS_INLINE void _xmpz_mul_nttcrt(xmpz_t& dst, const xmpz_t& left_op,
                                          const xmpz_t& right_op,
                                          const cuxmp_limb_t* primes,
                                          const cuxmp_len_t n_primes) {
    assert(left_op.n == right_op.n);
    assert(std::exp2(std::log2(left_op.n)) ==
           left_op.n);  // make sure N is a power of 2

    // use vector bc it automatically default constructs each number
    // and wont be included in any C header files.
    std::vector<xmpz_t> dst_parts(n_primes);

    xmpz_t left_ntt;
    xmpz_t right_ntt;
    xmpz_t ntt_product;

    // for each prime
    for (cuxmp_len_t p_i = 0; p_i < n_primes; p_i++) {
        cuxmp_limb_t p = primes[p_i];
        // take NTT of left and right
        _xmpz_ui32ntt(left_ntt, left_op, p);
        _xmpz_ui32ntt(right_ntt, right_op, p);
        // perform the pointwise multiplication mod p
        _xmpz_pointwise_mul_mod(ntt_product, left_ntt, right_ntt, p);
        // do the inverse ntt and put it in the destination
        _xmpz_ui32invntt(dst_parts[p_i], ntt_product, p);
    }
    // i think its possible for dst to not be large enough when using more than
    // one prime.
    _xmpz_pad_eq(dst,
                 left_op.n + (sizeof(cuxmp_crt_coef_t) + sizeof(cuxmp_sum_t)) /
                                 sizeof(cuxmp_limb_t));

    // reconstruct using crt
    _xmpz_construct_coef_crt(dst, dst_parts.data(), left_op.n, primes,
                             n_primes);
}

// Schönhagen-Strassen

CUXMP_ALWAYS_INLINE void _xmpz_xmpzntt(xmpz_arr_t dst, const xmpz_arr_t src,
                                       const xmpz_t& m) {
    // idk
}

CUXMP_ALWAYS_INLINE void _xmpz_xmpzinvntt(xmpz_arr_t dst, const xmpz_arr_t src,
                                          const xmpz_t& m) {}

CUXMP_ALWAYS_INLINE void _xmpz_ss_concat_chunks(xmpz_t& dst,
                                                const xmpz_arr_t src) {
    // probably not too difficult
}

CUXMP_ALWAYS_INLINE void _xmpz_ss_find_m(xmpz_t& dst, const cuxmp_len_t n) {
    // ?
}

CUXMP_ALWAYS_INLINE cuxmp_len_t _xmpz_ss_find_k(const cuxmp_len_t chunk_size,
                                                const cuxmp_len_t left_n,
                                                const cuxmp_len_t right_n) {
    // must return a power of 2

    cuxmp_len_t k_left = (left_n + chunk_size - 1) / chunk_size;
    cuxmp_len_t k_right = (right_n + chunk_size - 1) / chunk_size;

    cuxmp_len_t k = 1;
    while (k < k_left + k_right - 1) {
        k <<= 1;
    }

    return k;
}

// dst must be an array of initialized xmpz_ts of length k
CUXMP_ALWAYS_INLINE void _xmpz_ss_make_chunks(xmpz_arr_t dst, const xmpz_t& src,
                                              const cuxmp_len_t k) {
    // make k even-sized chunks of limbs out of limbs of src
    assert(src.n % k == 0);
    cuxmp_len_t chunk_size = src.n / k;

    for (cuxmp_len_t i = 0; i < k; i++) {
        dst[i].reserve(chunk_size);
        memcpy(dst[i].limbs, &src.limbs[i * chunk_size],
               k * sizeof(cuxmp_limb_t));
        dst[i].n = chunk_size;
    }
}

CUXMP_ALWAYS_INLINE void _xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                                      const xmpz_t& right_op) {
    // _n is N before padding
    cuxmp_len_t _n = std::max(left_op.n, right_op.n);

    // find chunk_size
    cuxmp_len_t chunk_size = (cuxmp_len_t)std::ceil(std::sqrt(_n));
    // find number of chunks k
    cuxmp_len_t k = _xmpz_ss_find_k(chunk_size, left_op.n, right_op.n);
    cuxmp_len_t n = k * chunk_size;

    // create copies of left and right and pad them to n
    xmpz_t left_padded;
    xmpz_t right_padded;

    _xmpz_pad(left_padded, left_op, n);
    _xmpz_pad(right_padded, right_op, n);

    // find modulo m
    xmpz_t m;
    _xmpz_ss_find_m(m, n);

    // create chunks
    std::vector<xmpz_t> left_chunks(k);
    std::vector<xmpz_t> right_chunks(k);

    _xmpz_ss_make_chunks(left_chunks.data(), left_padded, k);
    _xmpz_ss_make_chunks(right_chunks.data(), right_padded, k);

    // take xmpz NTT of chunks
    std::vector<xmpz_t> left_chunks_ntt(k);
    std::vector<xmpz_t> right_chunks_ntt(k);

    _xmpz_xmpzntt(left_chunks_ntt.data(), left_chunks.data(), m);
    _xmpz_xmpzntt(right_chunks_ntt.data(), right_chunks.data(), m);

    std::vector<xmpz_t> dst_chunks_ntt(k);

    // recursive multiplication
    for (cuxmp_len_t chunk = 0; chunk < k; chunk++) {
        xmpz_mul(dst_chunks_ntt[chunk], left_chunks_ntt[chunk],
                 right_chunks_ntt[chunk]);
        // have not implemented any modulo
        // xmpz_mod_eq(dst_chunks[chunk], m)
    }

    // do inverse xmpz NTT
    std::vector<xmpz_t> dst_chunks(k);

    _xmpz_xmpzinvntt(dst_chunks.data(), dst_chunks_ntt.data(), m);

    // generate final product

    _xmpz_ss_concat_chunks(dst, dst_chunks.data());
}

// left.n + right.n must be reserved in dst
CUXMP_ALWAYS_INLINE void _xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                                      const xmpz_t& right_op) {
    assert(dst.reserved >= left_op.n + right_op.n);

    dst.n = left_op.n + right_op.n;

    // zero all needed limbs of dst
    std::memset(dst.limbs, 0, (left_op.n + right_op.n) * sizeof(cuxmp_limb_t));

    // schoolbook multiply
    for (cuxmp_len_t i = 0; i < left_op.n; i++) {
        cuxmp_sum_t carry = 0;
        cuxmp_limb_t left_limb = left_op.limbs[i];

        for (cuxmp_len_t j = 0; j < right_op.n; j++) {
            cuxmp_len_t index = i + j;
            cuxmp_limb_t right_limb = right_op.limbs[j];
            cuxmp_sum_t acc = (cuxmp_sum_t)dst.limbs[index] +
                              (cuxmp_sum_t)left_limb * (cuxmp_sum_t)right_limb +
                              carry;
            dst.limbs[index] = (cuxmp_limb_t)acc;
            carry = (cuxmp_sum_t)(acc >> CUXMP_LIMB_BITS);
        }

        // final carry
        dst.limbs[i + right_op.n] = (cuxmp_limb_t)carry;
    }

    while (dst.n > 1 && dst.limbs[dst.n - 1] == 0) dst.n--;
}

}  // namespace xmp
