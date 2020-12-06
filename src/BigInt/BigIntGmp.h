#pragma once
#include <iostream>
#include <string>
#include <gmp.h>
#include "common.h"
#include "../Utility/bitManipulation.h"

struct BigIntGmp {
    mpz_t data;

    BigIntGmp() {
        mpz_init(data);
    }
    BigIntGmp(const std::string& str) : BigIntGmp() {
        mpz_set_str(data, str.c_str(), 10);
    }
    BigIntGmp(uint64_t value) : BigIntGmp(std::to_string(value)) {
    }
    BigIntGmp(const BigIntGmp& other) {
        mpz_init_set(data, other.data);
    }
    ~BigIntGmp() {
        mpz_clear(data);
    }
    BigIntGmp& operator=(const std::string& str) {
        mpz_set_str(data, str.c_str(), 10);
        return *this;
    }
    BigIntGmp& operator=(const BigIntGmp& other) {
        mpz_set(data, other.data);
        return *this;
    }
    uint64_t operator[](int index) const {
        return mpz_getlimbn(data, index);
    }
};
template<> struct BigIntParseImpl<BigIntGmp> {
    static BigIntGmp parse(const std::string& str) {
        return BigIntGmp{ str };
    }
};
BigIntGmp operator ""_gmp(const char* str) {
    return BigIntGmp(std::string(str));
}
int bitSize(const BigIntGmp& a) {
    return (int)mpz_sizeinbase(a.data, 2);
}
int sizeInLimbs(const BigIntGmp& a) {
    return (int)mpz_size(a.data);
}
int sizeInLimbs(const BarretReductionMod<BigIntGmp>& a) {
    return (int)mpz_size(a.mod.data);
}
int sizeInLimbs(const MontgomeryReductionMod<BigIntGmp>& a) {
    return (int)mpz_size(a.mod.data);
}
void swap(BigIntGmp& a, BigIntGmp& b) {
    std::swap(a.data, b.data);
}
void add(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_add(r.data, a.data, b.data);
}
void sub(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_sub(r.data, a.data, b.data);
}
void absSub(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    sub(r, a, b);
    mpz_abs(r.data, r.data);
}
void mul(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_mul(r.data, a.data, b.data);
}
void sqr(BigIntGmp& r, const BigIntGmp& a) {
    mpz_mul(r.data, a.data, a.data);
}
void div(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_div(r.data, a.data, b.data);
}
void mod(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_mod(r.data, a.data, b.data);
}
void shl(BigIntGmp& r, const BigIntGmp& a, uint32_t c) {
    mpz_mul_2exp(r.data, a.data, c);
}
void shr(BigIntGmp& r, const BigIntGmp& a, uint32_t c) {
    mpz_div_2exp(r.data, a.data, c);
}
int cmp(const BigIntGmp& a, const BigIntGmp& b) {
    return mpz_cmp(a.data, b.data);
}
int sign(const BigIntGmp& a) {
    return mpz_sgn(a.data);
}
void assign(BigIntGmp& dest, const BigIntGmp& src) {
    mpz_set(dest.data, src.data);
}
void assign(BigIntGmp& dest, const uint32_t& src) {
    mpz_set_ui(dest.data, src);
}
void assign(BigIntGmp& dest, const int& src) {
    mpz_set_ui(dest.data, (unsigned)src);
}
void modAdd(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BigIntGmp& m) {
    add(r, a, b);
    if (cmp(r, m) >= 0) {
        sub(r, r, m);
    }
}
void modSub(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BigIntGmp& m) {
    sub(r, a, b);
    if (sign(r) < 0) {
        add(r, r, m);
    }
}
void modMul(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BigIntGmp& m) {
    mul(r, a, b);
    mod(r, r, m);
}
void modSqr(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& m) {
    sqr(r, a);
    mod(r, r, m);
}
void modNeg(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& m) {
    if (sign(a) == 0) {
        assign(r, 0);
    } else {
        sub(r, m, a);
    }
}
void modDbl(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& m) {
    modAdd(r, a, a, m);
}
void modInv(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& m) {
    mpz_invert(r.data, a.data, m.data);
}
void gcd(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b) {
    mpz_gcd(r.data, a.data, b.data);
}
void gcd(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& b) {
    mpz_gcd(r.data, a.data, b.mod.data);
}
void gcd(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    mpz_gcd(r.data, a.data, m.mod.data);
}
bool isZero(const BigIntGmp& a) {
    return mpz_sgn(a.data) == 0;
}
void operator+=(BigIntGmp& r, int a) {
    mpz_add_ui(r.data, r.data, (unsigned long)a);
}
bool operator>(const BigIntGmp& a, int b) {
    return mpz_cmp_si(a.data, b) > 0;
}
bool operator==(const BigIntGmp& a, const BigIntGmp& b) {
    return mpz_cmp(a.data, b.data) == 0;
}
bool operator==(const BarretReductionMod<BigIntGmp>& a, const BarretReductionMod<BigIntGmp>& b) {
    return a.mod == b.mod;
}
bool operator==(const MontgomeryReductionMod<BigIntGmp>& a, const MontgomeryReductionMod<BigIntGmp>& b) {
    return a.mod == b.mod;
}
bool operator!=(const BigIntGmp& a, const BigIntGmp& b) {
    return mpz_cmp(a.data, b.data) != 0;
}
bool operator<(const BigIntGmp& a, const BigIntGmp& b) {
    return mpz_cmp(a.data, b.data) < 0;
}
bool operator!=(const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& b) {
    return a != b.mod;
}
bool operator!=(const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& b) {
    return a != b.mod;
}
BigIntGmp operator*(const BigIntGmp& a, const BigIntGmp& b) {
    BigIntGmp r;
    mul(r, a, b);
    return r;
}

void mod(BigIntGmp& u, const BigIntGmp& x, const MontgomeryReductionMod<BigIntGmp>& m) {
    mpz_mod_2exp(m.s.data, x.data, m.b * 64);
    mul(m.t, m.s, m.k);
    mpz_mod_2exp(m.s.data, m.t.data, m.b * 64);
    mul(m.t, m.s, m.mod);
    add(m.t, m.t, x);
    mpz_div_2exp(u.data, m.t.data, m.b * 64);
    if (mpz_cmp(u.data, m.mod.data) >= 0) {
        sub(u, u, m.mod);
    }
}
MontgomeryReductionMod<BigIntGmp> getMontgomeryReductionMod(const BigIntGmp& m) {
    MontgomeryReductionMod<BigIntGmp> result;
    result.mod = m;
    result.b = (uint32_t)mpz_size(m.data);
    BigIntGmp r, rInv;
    BigIntGmp one = 1;
    mpz_mul_2exp(r.data, one.data, result.b * 64);
    modInv(rInv, r, m);
    mpz_mul_2exp(r.data, rInv.data, result.b * 64);
    sub(r, r, one);
    div(result.k, r, m);
    return result;
}
void convertToMontgomeryForm(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    BigIntGmp mulRes;
    shl(mulRes, a, m.b * 64);
    mod(r, mulRes, m.mod);
}
BigIntGmp getValueInMontgomeryForm(const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    BigIntGmp r;
    convertToMontgomeryForm(r, a, m);
    return r;
}
BigIntGmp operator%(const BigIntGmp& a, const BigIntGmp& m) {
    BigIntGmp r;
    mod(r, a, m);
    return r;
}
BigIntGmp operator%(const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    BigIntGmp r;
    mod(r, a, m);
    return r;
}
void modAdd(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const MontgomeryReductionMod<BigIntGmp>& m) {
    modAdd(r, a, b, m.mod);
}
void modSub(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const MontgomeryReductionMod<BigIntGmp>& m) {
    modSub(r, a, b, m.mod);
}
void modMul(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const MontgomeryReductionMod<BigIntGmp>& m) {
    mul(r, a, b);
    mod(r, r, m);
}
void modSqr(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    sqr(r, a);
    mod(r, r, m);
}
void modNeg(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    modNeg(r, a, m.mod);
}
void modDbl(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    modAdd(r, a, a, m);
}
void modInv(BigIntGmp& r, const BigIntGmp& a, const MontgomeryReductionMod<BigIntGmp>& m) {
    modInv(r, a, m.mod);
}


BarretReductionMod<BigIntGmp> getBarretReductionMod(const BigIntGmp& n) {
    BarretReductionMod<BigIntGmp> result;
    result.mod = n;
    auto nSize = static_cast<int>(mpz_size(n.data));
    result.k = (static_cast<int>(sizeInBits(n[nSize-1])) + 64*(nSize-1)) * 2;

    BigIntGmp powK = 1;
    mpz_mul_2exp(powK.data, powK.data, result.k);
    result.R = 1;
    mpz_div(result.R.data, powK.data, n.data);

    return result;
}
void mod(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& b) {
    if (a < b.mod) {
        r = a;
        return;
    }
    mul(b.q_tmp, a, b.R);
    mpz_div_2exp(b.q_tmp.data, b.q_tmp.data, b.k);
    mul(b.q_tmp, b.q_tmp, b.mod);
    sub(r, a, b.q_tmp);
    if (cmp(r, b.mod) >= 0) {
        sub(r, r, b.mod);
    }
}
void modAdd(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BarretReductionMod<BigIntGmp>& m) {
    modAdd(r, a, b, m.mod);
}
void modSub(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BarretReductionMod<BigIntGmp>& m) {
    modSub(r, a, b, m.mod);
}
void modMul(BigIntGmp& r, const BigIntGmp& a, const BigIntGmp& b, const BarretReductionMod<BigIntGmp>& m) {
    mul(r, a, b);
    mod(r, r, m);
}
void modSqr(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& m) {
    sqr(r, a);
    mod(r, r, m);
}
void modNeg(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& m) {
    modNeg(r, a, m.mod);
}
void modDbl(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& m) {
    modAdd(r, a, a, m);
}
void modInv(BigIntGmp& r, const BigIntGmp& a, const BarretReductionMod<BigIntGmp>& m) {
    modInv(r, a, m.mod);
}
void powMod(BigIntGmp& r, const BigIntGmp& base, const BigIntGmp& exp, const BigIntGmp& mod) {
    mpz_powm(r.data, base.data, exp.data, mod.data);
}
void pow2Mod(BigIntGmp& r, const BigIntGmp& exp, const BigIntGmp& mod) {
    const static BigIntGmp v2 = 2;
    mpz_powm(r.data, v2.data, exp.data, mod.data);
}

BigIntGmp getModValue(const BigIntGmp& m) {
    return m;
}
BigIntGmp getModValue(const BarretReductionMod<BigIntGmp>& m) {
    return m.mod;
}
BigIntGmp getModValue(const MontgomeryReductionMod<BigIntGmp>& m) {
    return m.mod;
}
BigIntGmp convertToValue(const BigIntGmp& value, [[maybe_unused]] const BigIntGmp& dummy_mod) {
    return value;
}
BigIntGmp convertToValue(const BigIntGmp& value, [[maybe_unused]] const BarretReductionMod<BigIntGmp>& dummy_mod) {
    return value;
}
BigIntGmp convertToValue(const BigIntGmp& value, const MontgomeryReductionMod<BigIntGmp>& mod) {
    return value % mod;
}
BigIntGmp getConstant(uint64_t a, const BigIntGmp& m) {
    return BigIntGmp{ a } % m;
}
BigIntGmp getConstant(uint64_t a, const BarretReductionMod<BigIntGmp>& m) {
    return BigIntGmp{ a } % m.mod;
}
BigIntGmp getConstant(uint64_t a, const MontgomeryReductionMod<BigIntGmp>& m) {
    return getValueInMontgomeryForm(BigIntGmp{ a } % m.mod, m);
}

uint32_t mod(const BigIntGmp& a, uint32_t m, [[maybe_unused]] int tableEntryId) {
    return (uint32_t)mpn_mod_1(a.data->_mp_d, sizeInLimbs(a), m);
}


std::string toString(const BigIntGmp& a) {
    char buf[1000] = { 0 };
    mpz_get_str(buf, 10, a.data);
    return std::string(buf);
}
std::ostream& operator<<(std::ostream& out, const BigIntGmp& a) {
    return out << toString(a);
}

namespace std {
    template <> struct hash<BigIntGmp> {
        std::size_t operator()(const BigIntGmp& k) const {
            return hash<uint64_t>()(mpz_getlimbn(k.data, 0));
        }
    };
    template <> struct hash<BarretReductionMod<BigIntGmp>> {
        std::size_t operator()(const BarretReductionMod<BigIntGmp>& k) const {
            return hash<uint64_t>()(mpz_getlimbn(k.mod.data, 0));
        }
    };
    template <> struct hash<MontgomeryReductionMod<BigIntGmp>> {
        std::size_t operator()(const MontgomeryReductionMod<BigIntGmp>& k) const {
            return hash<uint64_t>()(mpz_getlimbn(k.mod.data, 0));
        }
    };
}
