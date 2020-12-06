#pragma once

#include "../Utility/compilerMacros.h"
#include "../Utility/debugAssert.h"

#if defined(_WIN32)
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif

/*
    Fast 64-bit * 64-bit multiplication functions returning 128 bit (high:low) result
*/
void mul128(uint64_t& highRes, uint64_t& lowRes, const uint64_t a, const uint64_t b) {
#ifdef COMPILER_MSVC
    lowRes = _umul128(a, b, &highRes);
#else
    __uint128_t product = (__uint128_t)a * (__uint128_t)b;
    highRes = (uint64_t)(product >> 64);
    lowRes = (uint64_t)(product & 0xFFFFFFFFFFFFFFFF);
#endif
}


/*
    128-bit / 64-bit division function
*/
uint64_t div128(uint64_t highDiv, uint64_t lowDiv, uint64_t divisor) {
#ifdef COMPILER_MSVC
    uint64_t remainder;
    return _udiv128(highDiv, lowDiv, divisor, &remainder);
#elif defined(__x86_64__)
    uint64_t result;
    uint64_t remainder;
    __asm__(
        "divq %[divisor]"
        : "=a"(result), "=d"(remainder)
        : [divisor] "r"(divisor), "a"(lowDiv), "d"(highDiv)
    );
    return result;
#endif
}


/*
    (r + out_carry) = (a + b + in_carry). in_carry can only be 0 or 1.
*/
uint8_t addCarry(uint8_t carry, uint64_t& r, uint64_t a, uint64_t b) {
    debugAssert(carry == 0 || carry == 1);
#if defined(_WIN32)
    return _addcarry_u64(carry, a, b, &r);
#elif defined(COMPILER_CLANG)
    unsigned long long carryOut;
    r = __builtin_addcll(a, b, carry, &carryOut);
    return carryOut;
#else
    r = a + b + carry;
    if (carry)
        carry = (r <= a);
    else
        carry = (r < a);
    return carry;
#endif
}

/*
    (r - out_borrow) = (a - b - in_borrow). in_borrow can only be 0 or 1.
*/
uint8_t subBorrow(uint8_t borrow, uint64_t& r, uint64_t a, uint64_t b) {
    debugAssert(borrow == 0 || borrow == 1);
#if defined(_WIN32)
    return _subborrow_u64(borrow, a, b, &r);
#elif defined(COMPILER_CLANG)
    unsigned long long borrowOut;
    r = __builtin_subcll(a, b, borrow, &borrowOut);
    return borrowOut;
#else
    r = a - b - borrow;
    if (borrow)
        borrow = (r >= a);
    else
        borrow = (r > a);
    return borrow;
#endif
}
