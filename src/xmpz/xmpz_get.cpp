#include "xmpz/xmpz_get.hpp"

using namespace xmp;

// getters
char* xmp::xmpz_get_str(const xmpz_t& src, cuxmp_base_t base,
                        cuxmp_len_t& out_len) {
    return _xmpz_get_str(src, base, out_len);
}
