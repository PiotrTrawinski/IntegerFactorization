#pragma once

#include "compilerMacros.h"
#include <cstdint>

#if defined(COMPILER_MSVC)
#define NOMINMAX
#include <Windows.h>
#endif

constexpr uint64_t MaxU64 = ~0ull;

uint64_t sizeInBits(uint64_t a) {
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    return 64ull - __builtin_clzll(a);
#elif defined(COMPILER_MSVC)
    DWORD index;
    _BitScanReverse64(&index, a);
    return index + 1ull;
#else
    static_assert(false, "unsupported compiler. Cannot compute sizeInBits");
#endif
}

uint64_t mostSignificantBit(uint64_t a) {
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    return 1ull << (63 - __builtin_clzll(a));
#elif defined(COMPILER_MSVC)
    DWORD index;
    _BitScanReverse64(&index, a);
    return 1ull << index;
#else
    static_assert(false, "unsupported compiler. Cannot compute mostSignificantBit");
#endif
}

uint64_t leadingZeroBitsCount(uint64_t a) {
    return 64 - sizeInBits(a);
}

uint64_t leadingSetBitCount(uint64_t a) {
    return leadingZeroBitsCount(~a);
}

uint64_t bitSetCount(uint64_t a) {
#if defined(_WIN32)
    return __popcnt64(a);
#else
    return __builtin_popcountll(a);
#endif
}

uint64_t trailingZeroBitCount(uint64_t a) {
#if defined(_WIN32)
    return __popcnt64((a & -a) - 1);
#else
    return __builtin_ctzll(a);
#endif
}

