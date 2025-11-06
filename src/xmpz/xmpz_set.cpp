#include "xmpz/xmpz_set.hpp"

#include <cstring>

using namespace xmp;

// setters
void xmp::xmpz_set_u32(xmpz_t& dst, const uint32_t value) {
    if constexpr (sizeof(uint32_t) == sizeof(cuxmp_limb_t)) {
        xmpz_set_limb(dst, (cuxmp_limb_t)value);
    } else {
        _xmpz_init_ui<uint32_t>(dst, value);
    }
}

void xmp::xmpz_set_u64(xmpz_t& dst, const uint64_t value) {
    if constexpr (sizeof(uint64_t) == sizeof(cuxmp_limb_t) * 2) {
        dst.reserve(2);
        dst.limbs[0] = (cuxmp_limb_t)(value & CUXMP_LIMB_MAX);
        dst.limbs[1] =
            (cuxmp_limb_t)((value >> CUXMP_LIMB_BITS) & CUXMP_LIMB_MAX);
        dst.n = dst.limbs[1] == 0 ? 1 : 2;
    } else {
        _xmpz_init_ui<uint64_t>(dst, value);
    }
}

void xmp::xmpz_set_limb(xmpz_t& dst, const cuxmp_limb_t& src) {
    dst.reserve(1);
    dst.limbs[0] = src;
    dst.n = 1;
}

void xmp::xmpz_set(xmpz_t& result, const xmpz_t& src) {
    result.reserve(src.n);
    memcpy(result.limbs, src.limbs, src.n * sizeof(cuxmp_limb_t));

    result.n = src.n;
}

void xmp::xmpz_copy(xmpz_t& dst, const xmpz_t& src) {
    dst.limbs = (cuxmp_limb_t*)malloc(src.reserved * sizeof(cuxmp_limb_t));
    memcpy(dst.limbs, src.limbs, src.reserved * sizeof(cuxmp_limb_t));
    dst.n = src.n;
    dst.reserved = src.reserved;
}
