#pragma once

#include <cassert>
#include <cstring>
#include <vector>

#include "xmpz.hpp"

namespace xmp {

cuxmp_stat_t xmpz_mul_sb(xmpz_t& dst, const xmpz_t& left_op,
                         const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul_nttcrt(xmpz_t& dst, const xmpz_t left_op,
                             const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                         const xmpz_t& right_op);

cuxmp_stat_t xmpz_mul(xmpz_t& dst, const xmpz_t& left_op,
                      const xmpz_t& right_op);

// NTT-CRT

CUXMP_ALWAYS_INLINE void _xmpz_ui32ntt(xmpz_t& dst, const xmpz_t& src,
                                       const cuxmp_limb_t prime) {
    // no idea
}

CUXMP_ALWAYS_INLINE void _xmpz_ui32invntt(xmpz_t& dst, const xmpz_t& src,
                                          const cuxmp_limb_t prime) {
    // ...
}

// requires at least parts[0].n be reserved in dst. all parts must have same N.
CUXMP_ALWAYS_INLINE void _xmpz_crt(xmpz_t& dst, const xmpz_t* parts,
                                   const cuxmp_limb_t* primes,
                                   const cuxmp_len_t n_primes) {
    // shouldnt be too difficult ??
}

// requires at least coef_src.n be reserved in dst
CUXMP_ALWAYS_INLINE void _xmpz_construct_from_coef(xmpz_t& dst,
                                                   const xmpz_t& coef_src) {
    // probably not too difficult either ?
}

CUXMP_ALWAYS_INLINE cuxmp_limb_t* _xmpz_ntt_find_primes(
    const cuxmp_len_t n, cuxmp_len_t& out_prime_n) {
    // perhaps just some if statements that return different predetermined
    // primes up to N = 2^64 ??
}

CUXMP_ALWAYS_INLINE void _xmpz_pointwise_mul_mod(xmpz_t& dst,
                                                 const xmpz_t& left_op,
                                                 const xmpz_t& right_op,
                                                 const cuxmp_limb_t m) {
    // probably not too difficult
}

CUXMP_ALWAYS_INLINE void _xmpz_mul_nttcrt(xmpz_t& dst, const xmpz_t& left_op,
                                          const xmpz_t& right_op,
                                          const cuxmp_limb_t* primes,
                                          const cuxmp_len_t primes_n) {
    assert(left_op.n == right_op.n);

    // use vector bc it automatically default constructs each number
    // and wont be included in any C header files.
    std::vector<xmpz_t> dst_parts(primes_n);

    xmpz_t left_ntt;
    xmpz_t right_ntt;
    xmpz_t ntt_product;

    for (cuxmp_len_t p_i = 0; p_i < primes_n; p_i++) {
        cuxmp_limb_t p = primes[p_i];
        _xmpz_ui32ntt(left_ntt, left_op, p);
        _xmpz_ui32ntt(right_ntt, right_op, p);
        _xmpz_pointwise_mul_mod(ntt_product, left_ntt, right_ntt, p);
        _xmpz_ui32invntt(dst_parts[p_i], ntt_product, p);
    }

    xmpz_t dst_coef;
    dst_coef.reserve(left_op.n);
    _xmpz_crt(dst_coef, dst_parts.data(), primes, primes_n);
    dst.reserve(left_op.n);
    _xmpz_construct_from_coef(dst, dst_coef);
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
    // ...
}

CUXMP_ALWAYS_INLINE cuxmp_len_t _xmpz_ss_find_k(const cuxmp_len_t n) {
    // probably return sqrt(n) for now
}

CUXMP_ALWAYS_INLINE void _xmpz_ss_make_chunks(xmpz_arr_t dst, const xmpz_t& src,
                                              const cuxmp_len_t k) {
    // make k even-sized chunks of limbs out of limbs of src
}

CUXMP_ALWAYS_INLINE void _xmpz_mul_ss(xmpz_t& dst, const xmpz_t& left_op,
                                      const xmpz_t& right_op) {
    cuxmp_len_t n = std::max(left_op.n, right_op.n);

    cuxmp_len_t k = _xmpz_ss_find_k(n);

    // i think i need to pad left and right op such that left_op.n % k == 0. so
    // i can make k even-sized chunks

    xmpz_t m;  // modulo
    _xmpz_ss_find_m(m, n);

    std::vector<xmpz_t> left_chunks(k);
    std::vector<xmpz_t> right_chunks(k);

    _xmpz_ss_make_chunks(left_chunks.data(), left_op, k);
    _xmpz_ss_make_chunks(right_chunks.data(), right_op, k);

    std::vector<xmpz_t> left_chunks_ntt(k);
    std::vector<xmpz_t> right_chunks_ntt(k);

    _xmpz_xmpzntt(left_chunks_ntt.data(), left_chunks.data(), m);
    _xmpz_xmpzntt(right_chunks_ntt.data(), right_chunks.data(), m);

    std::vector<xmpz_t> dst_chunks_ntt(k);

    for (cuxmp_len_t chunk = 0; chunk < k; chunk++) {
        xmpz_mul(dst_chunks_ntt[chunk], left_chunks_ntt[chunk],
                 right_chunks_ntt[chunk]);
        // have not implemented any modulo
        // xmpz_mod_eq(dst_chunks[chunk], m)
    }

    std::vector<xmpz_t> dst_chunks(k);

    _xmpz_xmpzinvntt(dst_chunks.data(), dst_chunks_ntt.data(), m);

    dst.reserve(n);

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
