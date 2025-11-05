#include <iostream>
#include <string>

#include "xmp.hpp"

using namespace xmp;

int main() {
    xmpz_t test;
    xmpz_t test2;
    xmpz_t result;

    xmpz_set_u64(test, UINT64_MAX);
    xmpz_set_u64(test2, UINT64_MAX);

    xmpz_add(result, test, test2);
    xmpz_add(result, result, test);

    std::cout << test.n << " " << test.limbs[0] << " " << test.limbs[1] << "\n"
              << test2.n << " " << test2.limbs[0] << " " << test2.limbs[1]
              << "\n"
              << result.n << " " << result.reserved << " " << result.limbs[0]
              << " " << result.limbs[1] << "\n";

    cuxmp_len_t len;
    char* str = xmpz_get_str(result, 10, len);

    std::string s(str, str + len);

    std::cout << s << std::endl;
}
