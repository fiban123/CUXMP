#include <cstdint>

using cuxmp_limb_t = uint32_t;
using cuxmp_len_t = size_t;
using cuxmp_stat_t = uint8_t;
using cuxmp_sum_t = uint64_t;
using cuxmp_base_t = uint8_t;

// constants
constexpr cuxmp_sum_t CUXMP_LIMB_BITS = 8ULL * sizeof(cuxmp_limb_t);
constexpr cuxmp_len_t CUXMP_ADD_SB_MIN = 0;
constexpr cuxmp_limb_t CUXMP_LIMB_MAX = UINT32_MAX;
constexpr cuxmp_sum_t CUXMP_SUM_MAX = UINT64_MAX;
constexpr char CUXMP_STR_CHARSET[] =
    "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

// statuses
constexpr cuxmp_stat_t CUXMP_STAT_OK = 0;
constexpr cuxmp_stat_t CUXMP_STAT_OVERFLOW = 0;

// shortcut macros
#define CUXMP_ALWAYS_INLINE inline __attribute__((always_inline))
