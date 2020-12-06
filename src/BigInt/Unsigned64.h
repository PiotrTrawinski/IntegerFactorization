#pragma once
#include <cstdint>
#include <string>
#include <numeric>
#include "common.h"

template<> struct BigIntParseImpl<uint64_t> {
    static uint64_t parse(const std::string& str) {
        return std::stoull(str);
    }
};
std::string toString(uint64_t n) {
    return std::to_string(n);
}

void assign(uint64_t& dest, const int& src) {
    dest = src;
}
void assign(uint64_t& dest, uint64_t src) {
    dest = src;
}
void add(uint64_t& r, uint64_t a, uint64_t b) {
    r = a + b;
}
void sub(uint64_t& r, uint64_t a, uint64_t b) {
    r = a - b;
}
void absSub(uint64_t& r, uint64_t a, uint64_t b) {
    if (a > b)
        r = a - b;
    else
        r = b - a;
}
void modAdd(uint64_t& r, uint64_t a, uint64_t b, uint64_t m) {
    r = a + b;
    if (r >= m) {
        r -= m;
    }
}
void modSub(uint64_t& r, uint64_t a, uint64_t b, uint64_t m) {
    auto tmp = a;
    r = a - b;
    if (r > tmp) {
        r += m;
    }
}
void modMul(uint64_t& r, uint64_t a, uint64_t b, uint64_t m) {
    r = (a * b) % m;
}
void modSqr(uint64_t& r, uint64_t a, uint64_t m) {
    r = (a * a) % m;
}
void modNeg(uint64_t& r, uint64_t a, uint64_t m) {
    modSub(r, 0, a, m);
}
void modDbl(uint64_t& r, uint64_t a, uint64_t m) {
    modAdd(r, a, a, m);
}
void modInv(int64_t& t, int64_t a, int64_t m) {
    t = 0;
    int64_t newT = 1;
    int64_t r = m;
    int64_t newR = a;

    while (newR != 0) {
        int64_t q = r / newR;
        int64_t tmp = t;
        t = newT;
        newT = tmp - q * newT;

        tmp = r;
        r = newR;
        newR = tmp - q * newR;
    }

    if (r > 1) {
        abort();
    }
    if (t < 0) {
        t += m;
    }
}
void gcd(uint64_t& r, const uint64_t& a, const uint64_t& b) {
    r = std::gcd(a, b);
}
void pow2Mod(uint64_t& r, uint64_t e, uint64_t mod) {
    r = 1;
    uint64_t a = 2;
    while (e > 0) {
        if (e & 1) r = r * a % mod;
        a = a * a % mod;
        e >>= 1;
    }
}
uint64_t convertToValue(uint64_t a, [[maybe_unused]] const uint64_t& dummy_mod) {
    return a;
}
uint64_t convertToValue(uint64_t a, [[maybe_unused]] const BarretReductionMod<uint64_t>& dummy_mod) {
    return a;
}
uint64_t convertToValue([[maybe_unused]] uint64_t a, [[maybe_unused]] const MontgomeryReductionMod<uint64_t>& mod) {
    debugAssert(false);
    return 0;
}
uint64_t getConstant(int a, [[maybe_unused]] const uint64_t& dummy_mod) {
    return static_cast<uint64_t>(a);
}
uint64_t getConstant(int a, [[maybe_unused]] const BarretReductionMod<uint64_t>& dummy_mod) {
    return static_cast<uint64_t>(a);
}
uint64_t getConstant([[maybe_unused]] int a, [[maybe_unused]] const MontgomeryReductionMod<uint64_t>& mod) {
    debugAssert(false);
    return 0;
    //return getValueInMontgomeryForm(static_cast<uint64_t>(a), mod);
}
