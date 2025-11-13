// TEST stuff
#pragma once

#include <stddef.h>

#include <cstdint>
#include <format>
#include <iostream>
#include <iterator>
#include <random>

#include "xmpz/xmpz.hpp"
#include "xmpz/xmpz_kernel.hpp"
#include "xmpz/xmpz_mul.hpp"

using namespace xmp;

#define RED "\e[0;31m"
#define GREEN "\e[0;32m"
#define WHITE "\e[0;97m"

inline int n_primes = 3;
inline int check_n = 1000;
inline std::string corr = GREEN "correct" WHITE;
inline std::string wrng = RED "WRONG !!!" WHITE;

template <typename T, typename U>
inline bool check(const T& value, const U& expected,
                  std::string_view msg = "") {
    if (value == expected) {
        // return true;
    }
    if (!msg.empty()) std::cout << msg << " -> ";
    std::cout << std::format("{} ({}, {})\n", (value == expected ? corr : wrng),
                             value, expected);

    return false;
}

inline std::mt19937_64 rng{67};

template <typename T>
inline T rand_val()
    requires std::is_trivially_copyable_v<T>
{
    T value;
    std::uniform_int_distribution<uint64_t> dist(0, UINT64_MAX);

    // Fill sizeof(T) bytes with random data
    unsigned char* p = reinterpret_cast<unsigned char*>(&value);
    for (size_t i = 0; i < sizeof(T); i += sizeof(uint64_t)) {
        uint64_t chunk = dist(rng);
        std::memcpy(p + i, &chunk, std::min(sizeof(uint64_t), sizeof(T) - i));
    }

    return value;
}

template <typename T>
inline std::vector<T> rand_arr(size_t n) {
    std::vector<T> vec(n);

    for (size_t i = 0; i < n; i++) {
        vec[i] = rand_val<T>();
    }

    return vec;
}

template <typename T>
inline void print_vec(std::vector<T> vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << ", ";
    }
    std::cout << "\n";
}

template <typename T>
inline void print_arr(T* ptr, size_t n) {
    for (size_t i = 0; i < n; i++) {
        std::cout << ptr[i] << ", ";
    }
    std::cout << "\n";
}

inline void test_const() {
    std::cout << "prime root test: \nall these should be exactly 1\n";

    for (int i = 0; i < n_primes; i++) {
        cuxmp_limb_t res =
            _ntt_powmod(ntt_primes[0].root, 1u << 24, ntt_primes[0]);

        check(res, 1, std::format("prime {}", i));
    }

    std::cout << "\ncrt inv constants check:\n";

    for (int i = 0; i < n_primes; i++) {
        for (int j = 0; j < n_primes; j++) {
            if (i < j) {
                cuxmp_limb_t res = _crt_modmul(_crt_const_inv[i][j],
                                               ntt_primes[i].p, ntt_primes[j]);

                check(res, 1, std::format("const {}, {}", i, j));
            } else {
                std::cout << std::format("skipping {}, {}\n", i, j);
            }
        }
    }

    std::cout << "\ncrt prod constants check:\n";

    check(_crt_const_p_prod[0], 1, "1");
    check(_crt_const_p_prod[1], ntt_primes[0].p, "2");
    check(_crt_const_p_prod[2], (__uint128_t)ntt_primes[0].p * ntt_primes[1].p,
          "3");

    std::cout << "\nmicro-kernel check:\n";

    xmp_ntt_prime tp = ntt_primes[0];
    uint32_t tpp = tp.p;

    auto nums1 = rand_arr<cuxmp_limb_t>(check_n);
    auto nums2 = rand_arr<cuxmp_limb_t>(check_n);

    for (int i = 0; i < check_n; i++) {
        nums1[i] %= tpp;
        nums2[i] %= tpp;

        // std::cout << nums1[i] << " " << nums2[i] << "\n";

        check(
            //
            _ntt_addmod(nums1[i], nums2[i], tpp),
            ((cuxmp_sum_t)nums1[i] + nums2[i]) % tpp,
            std::format("ntt addmod {}", i)
            //
        );

        check(
            //
            _ntt_submod(nums1[i], nums2[i], tpp),
            ((uint64_t)nums1[i] + tpp - (uint64_t)nums2[i]) % tpp,
            std::format("ntt submod {}", i)
            //
        );

        check(
            //
            _ntt_barret_mulmod(nums1[i], nums2[i], tp),
            ((uint64_t)nums1[i] * (uint64_t)nums2[i]) % tpp,
            std::format("ntt barret mulmod {}", i)
            //
        );

        cuxmp_limb_t a = nums1[i];
        cuxmp_limb_t w = nums2[i];

        cuxmp_limb_t w_s =
            (cuxmp_limb_t)(((__uint128_t)w << CUXMP_LIMB_BITS) / tpp);

        check(
            //
            _ntt_shoup_mulmod(a, w_s, w, tpp),
            (cuxmp_limb_t)(((__uint128_t)a * w) % tpp),
            std::format("ntt shoup mulmod {}", i)
            //
        );

        cuxmp_limb_t res = _ntt_powmod(nums1[i], nums2[i], tp);

        // std::cout << std::format("{}^{} mod {} = {}\n", nums1[i], nums2[i],
        //                         tp.p, res);
    }

    std::cout << "\nntt check \n";

    for (int i = 0; i < check_n; i++) {
        xmpz_t a, b, c;

        a.reserve(1024);
        a.n = 1024;

        std::vector<cuxmp_limb_t> limbs = rand_arr<cuxmp_limb_t>(1024);

        for (cuxmp_limb_t& limb : limbs) {
            limb %= tpp;
        }

        memcpy(a.limbs, limbs.data(), sizeof(cuxmp_limb_t) * 1024);

        // print_arr(a.limbs, a.n);
        _xmpz_ui32ntt_fwd(b, a, tp);

        //
        // print_arr(b.limbs, b.n);

        _xmpz_ui32ntt_inv(c, b, tp);

        // print_arr(c.limbs, c.n);

        check(
            //
            memcmp(limbs.data(), c.limbs, 1024 * sizeof(cuxmp_limb_t)), 0,
            std::format("ntt {}", i));
    }

    std::cout << "\nCRT test\n";

    for (int i = 0; i < check_n; i++) {
        auto residues = rand_arr<cuxmp_limb_t>(n_primes);
        // std::vector<cuxmp_limb_t> residues = {2128530233, 1661722638,
        // 28532098};
        //  std::cout << "residues: ";
        //  print_vec(residues);

        for (int j = 0; j < n_primes; j++) {
            residues[j] %= ntt_primes[j].p;
        }

        cuxmp_crt_coef_t x =
            _crt_solve_kernel(residues.data(), ntt_primes, n_primes);

        // verify
        for (int j = 0; j < n_primes; j++) {
            if (!check(
                    //
                    (cuxmp_limb_t)(x % ntt_primes[j].p), residues[j],
                    std::format("CRT solver {} {}", i, j)
                    //
                    )) {
                print_vec(residues);
            }

            cuxmp_limb_t c = (cuxmp_limb_t)(x % ntt_primes[i].p);
        }
    }
}
