#include <cstdint>
#include <iostream>
#include <string>

#include "xmp.hpp"

using namespace xmp;

int find_primes();

int main() {
    xmpz_t test;
    xmpz_t test2;
    xmpz_t result;

    xmpz_set_u64(test, UINT64_MAX);
    xmpz_set_u64(test2, UINT64_MAX);

    xmpz_mul_sb(result, test, test2);

    std::cout << test.n << " " << test.limbs[0] << " " << test.limbs[1] << "\n"
              << test2.n << " " << test2.limbs[0] << " " << test2.limbs[1]
              << "\n"
              << result.n << " " << result.reserved << " " << result.limbs[0]
              << " " << result.limbs[1] << "\n";

    cuxmp_len_t len;
    char* str = xmpz_get_str(result, 10, len);

    std::string s(str, str + len);

    std::cout << s << std::endl;

    find_primes();
}

using namespace std;

// Simple 64-bit modular multiplication (avoids overflow)
static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t mod) {
    __uint128_t res = (__uint128_t)a * b;
    return (uint64_t)(res % mod);
}

// Fast modular exponentiation
uint64_t powmod(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    while (exp) {
        if (exp & 1) result = mulmod(result, base, mod);
        base = mulmod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// Miller-Rabin primality test for 64-bit integers
bool is_prime(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
        if (n % p == 0) return n == p;

    uint64_t d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        s++;
    }

    // Deterministic bases for 64-bit range
    for (uint64_t a : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {
        if (a % n == 0) continue;
        uint64_t x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (uint64_t r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }
    return true;
}

int find_primes() {
    constexpr uint64_t MAX_BITS = 32;
    constexpr uint64_t N = 1ULL << 24;  // 2^20
    const uint64_t MAX_P = 1ULL << MAX_BITS;

    cout << "Searching for 30-bit NTT primes of the form p = m*2^20 + 1...\n";

    for (uint64_t m = 1;; m++) {
        uint64_t p = m * N + 1;
        if (p >= MAX_P) break;
        if (is_prime(p)) {
            cout << "m = " << m << "  =>  p = " << p << "\n";
        }
    }

    return 0;
}
