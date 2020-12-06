#pragma once
#include "common.h"
#include "../Utility/compilerMacros.h"
#include "../Utility/debugAssert.h"
#include "../Utility/bitManipulation.h"
#include "kernels.h"
#include <iostream>
#include <string>
#include <cassert>
#include <cstring>
#include <array>
#include <algorithm>


/*
    Big Integer class with fixed stack storage of values.
    The storage capacity is fixed, but the actual size of values dynamically changes between [1, capacity].
    It's the user's responsibility to not overflow the capacity.

    It stores natural + {0} numbers only. 
    If the result of subtraction is negative it stores the absolute value and returns information about the sign.
*/
template<int MaxCap=15> struct BigIntMaxCap {
    // data:
    std::array<uint64_t, MaxCap> data;
    int32_t size_;
    
    // methods:
    BigIntMaxCap(uint64_t value) : size_(1) { data[0] = value; }
    BigIntMaxCap() : BigIntMaxCap(0) {}
    BigIntMaxCap(const std::string& str) { parseString(str); }

    template<int OtherMaxCap> BigIntMaxCap(const BigIntMaxCap<OtherMaxCap>& other) {
        debugAssert(MaxCap >= other.size_);
        memcpy(&data[0], &other.data[0], other.size_ * sizeof(uint64_t));
        size_ = other.size_;
    }
    template<int OtherMaxCap> BigIntMaxCap& operator=(const BigIntMaxCap<OtherMaxCap>& other) {
        debugAssert(MaxCap >= other.size_);
        memcpy(&data[0], &other.data[0], other.size_ * sizeof(uint64_t));
        size_ = other.size_;
        return *this;
    }
    template<int OtherMaxCap> BigIntMaxCap& operator=(uint64_t value) {
        size_ = 1;
        data[0] = value;
        return *this;
    }
    BigIntMaxCap& operator=(const std::string& str) {
        parseString(str);
        return *this;
    }

    uint64_t* ptr()             { return &data[0]; }
    const uint64_t* ptr() const { return &data[0]; }

    uint64_t& operator[](int index)             { return data[index]; }
    const uint64_t& operator[](int index) const { return data[index]; }

    BigIntBit bit(int index)                  { return BigIntBit(ptr() + index / 64, index % 64); }
    const BigIntConstBit bit(int index) const { return BigIntConstBit(ptr() + index / 64, index % 64); }

    int32_t size() const        { return size_; }
    uint32_t sizeInBits() const { return (size_ - 1) * 64 + static_cast<uint32_t>(::sizeInBits(data[size_ -1])); }

private:
    void parseString(const std::string& str) {
        size_ = 1;
        data[0] = 0;
        for (auto c : str) {
            mulAssign(*this, 10);
            add(*this, *this, static_cast<uint64_t>(c) - '0');
        }
    }
};


template<int MaxCap> struct BigIntParseImpl<BigIntMaxCap<MaxCap>> {
    static BigIntMaxCap<MaxCap> parse(const std::string& str) {
        return BigIntMaxCap<MaxCap> { str };
    }
};


// -------------------------------------
// Internal functions
// -------------------------------------
namespace details {
    template<int C> void setZero(BigIntMaxCap<C>& r, uint32_t count) {
        for (uint32_t i = 0; i < count; ++i) {
            r[i] = 0;
        }
    }
    template<int C> void shrinkSizeToFit(BigIntMaxCap<C>& r) {
        while (r.size() >= 1 && r[r.size() - 1] == 0) {
            r.size_ -= 1;
        }
    }
}


// -------------------------------------
// Public functions
// -------------------------------------
template<int C1, int C2> void assign(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a) {
    r = a;
}
template<int C> void assign(BigIntMaxCap<C>& r, const int& a) {
    r = a;
}
template<int C1, int C2> int cmp(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    //return bigIntDetails::cmp(a, b);
    return bigIntKernels::cmp(a.ptr(), b.ptr(), a.size(), b.size());
}
template<int C1, int C2, int C3> void add(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b) {
    r.size_ = std::max(a.size(), b.size());
    r[r.size_] = bigIntKernels::add(r.ptr(), a.ptr(), b.ptr(), a.size(), b.size());
    if (r[r.size_]) {
        r.size_ += 1;
    }
}


/*
    Returns 'true' if negative result.
    Returns 'false if positive result.
*/
template<int C1, int C2, int C3> bool sub(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b) {
    auto isNegative = bigIntKernels::sub(r.ptr(), a.ptr(), b.ptr(), a.size(), b.size());
    r.size_ = std::max(a.size(), b.size());
    details::shrinkSizeToFit(r);
    return isNegative;
}
template<int C1, int C2, int C3> void mul(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b) {
    if (r.ptr() == a.ptr()) {
        auto aCopy = a;
        if (r.ptr() == b.ptr()) {
            auto bCopy = b;
            bigIntKernels::mul(r.ptr(), aCopy.ptr(), bCopy.ptr(), a.size(), b.size());
        } else {
            bigIntKernels::mul(r.ptr(), aCopy.ptr(), b.ptr(), a.size(), b.size());
        }
    } else if (r.ptr() == b.ptr()) {
        auto bCopy = b;
        bigIntKernels::mul(r.ptr(), a.ptr(), bCopy.ptr(), a.size(), b.size());
    } else {
        bigIntKernels::mul(r.ptr(), a.ptr(), b.ptr(), a.size(), b.size());
    }
    r.size_ = a.size() + b.size();
    if (r[r.size_ - 1] == 0) {
        r.size_ -= 1;
    }
}
template<int C1, int C2> void sqr(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a) {
    uint64_t tmpStorage[C2 * 2];
    if (r.ptr() == a.ptr()) {
        auto aCopy = a;
        bigIntKernels::sqr(r.ptr(), aCopy.ptr(), a.size(), tmpStorage);
    } else {
        bigIntKernels::sqr(r.ptr(), a.ptr(), a.size(), tmpStorage);
    }
    r.size_ = a.size()*2;
    if (r[r.size_ - 1] == 0) {
        r.size_ -= 1;
    }
}
template<int C1, int C2> void shr(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint32_t c) {
    r.size_ = bigIntKernels::shr(r.ptr(), a.ptr(), c, r.size(), a.size());
    if (r[r.size_ - 1] == 0) {
        r.size_ -= 1;
    }
}
template<int C1, int C2> void shl(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint32_t c) {
    r.size_ = bigIntKernels::shl(r.ptr(), a.ptr(), c, a.size());
    if (r[r.size_ - 1] == 0) {
        r.size_ -= 1;
    }
}



/*
    Egyptian division algorithm (slow)
*/
template<int C1, int C2, int C3, int C4> void div(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, BigIntMaxCap<C4>& rem) {
    std::vector<BigIntMaxCap<C2>> doublings;
    BigIntMaxCap<C2> doubling = b;

    int i = 0;
    for (; ; ++i) {
        doublings.emplace_back(doubling);
        doubling <<= 1;
        if (doubling > a) {
            break;
        }
    }

    r = doublings.back();
    details::setZero(r, r.size());
    BigIntMaxCap<C2> accumulator = 0;
    BigIntMaxCap<C2> tmp = 0;
    for (; i >= 0; --i) {
        add(tmp, accumulator, doublings[i]);
        //tmp = accumulator + doublings[i];
        if (tmp <= a) {
            accumulator = tmp;
            r.bit(i) = 1;
        }
    }
    sub(rem, a, accumulator);
    //rem = a - accumulator;
}

template<int C1, int C2, int C3> void div(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b) {
    if (b.size() > a.size()) {
        r.data[0] = 0;
        r.size_ = 1;
        return;
    }
    BigIntMaxCap<C1> dummyRem;
    div(r, a, b, dummyRem);
    r.size_ = a.size() - b.size() + 1 - (a.size()-b.size() > 0 && r[a.size() - b.size()] == 0);
}


/*
    Modular operations using standard modulo
*/
template<int C1, int C2, int C3> void mod(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b) {
    BigIntMaxCap<C1> dummyDivResult;
    div(dummyDivResult, a, b, r);
}
template<int C1, int C2, int C3, int C4> void modAdd(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BigIntMaxCap<C4>& m) {
    add(r, a, b);
    if (r >= m) {
        sub(r, r, m);
    }
}
template<int C1, int C2, int C3, int C4> void modSub(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BigIntMaxCap<C4>& m) {
    auto isNegative = sub(r, a, b);
    if (isNegative) {
        sub(r, m, r);
    }
}
template<int C1, int C2, int C3, int C4> void modMul(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BigIntMaxCap<C4>& m) {
    mul(r, a, b);
    mod(r, r, m);
}
template<int C1, int C2, int C3> void modSqr(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& m) {
    sqr(r, a);
    mod(r, r, m);
}
template<int C1, int C2, int C3> void modNeg(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& m) {
    sub(r, m, a);
}
template<int C1, int C2, int C3> void modDbl(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& m) {
    shl(r, a, 1);
    if (r >= m) {
        sub(r, r, m);
    }
}
template<int C> void gcd(BigIntMaxCap<C>& r, const BigIntMaxCap<C>& a, const BigIntMaxCap<C>& b) {
    mpn_gcd(r.ptr(), (uint64_t*)a.ptr(), a.size(), (uint64_t*)b.ptr(), b.size());
    r.size_ = a.size();
    details::shrinkSizeToFit(r);
}
template<int C> void gcd(BigIntMaxCap<C>& r, const BigIntMaxCap<C>& a, const BarretReductionMod<BigIntMaxCap<C>>& b) {
    mpn_gcd(r.ptr(), (uint64_t*)a.ptr(), a.size(), (uint64_t*)b.mod.ptr(), b.mod.size());
    r.size_ = a.size();
    details::shrinkSizeToFit(r);
}



/*
    Modular operations using Barrett reduction. See 'BarretReductionMod' struct for more info.
*/
template<int C> BarretReductionMod<BigIntMaxCap<C>> getBarretReductionMod(const BigIntMaxCap<C>& n) {
    BarretReductionMod<BigIntMaxCap<C>> result;
    result.mod = n;
    result.k = (static_cast<int>(sizeInBits(n[n.size()-1])) + 64*(n.size()-1)) * 2;

    BigIntMaxCap<C> powK = 1;
    shl(powK, powK, result.k);
    result.R = 1;
    div(result.R, powK, n);

    return result;
}
template<int C1, int C2, int C3> void mod(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BarretReductionMod<BigIntMaxCap<C3>>& b) {
    if (a < b.mod) {
        r = a;
        return;
    }
    mul(b.q_tmp, a, b.R);
    shr(b.q_tmp, b.q_tmp, b.k);
    mul(b.q_tmp, b.q_tmp, b.mod);
    sub(r, a, b.q_tmp);
    if (r >= b.mod) {
        sub(r, r, b.mod);
    }
}
template<int C1, int C2, int C3, int C4> void modAdd(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BarretReductionMod<BigIntMaxCap<C4>>& m) {
    modAdd(r, a, b, m.mod);
}
template<int C1, int C2, int C3, int C4> void modSub(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BarretReductionMod<BigIntMaxCap<C4>>& m) {
    modSub(r, a, b, m.mod);
}
template<int C1, int C2, int C3, int C4> void modMul(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const BarretReductionMod<BigIntMaxCap<C4>>& m) {
    mul(r, a, b);
    mod(r, r, m);
}
template<int C1, int C2, int C3> void modSqr(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BarretReductionMod<BigIntMaxCap<C3>>& m) {
    sqr(r, a);
    mod(r, r, m);
}
template<int C1, int C2, int C3> void modNeg(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BarretReductionMod<BigIntMaxCap<C3>>& m) {
    modNeg(r, a, m.mod);
}
template<int C1, int C2, int C3> void modDbl(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BarretReductionMod<BigIntMaxCap<C3>>& m) {
    modDbl(r, a, m.mod);
}


template<int C> MontgomeryReductionMod<BigIntMaxCap<C>> getMontgomeryReductionMod(const BigIntMaxCap<C>& n) {
    MontgomeryReductionMod<BigIntMaxCap<C>> result;
    int32_t s = n.size();
    result.mod = n;
    result.b = s;
    BigIntMaxCap<C> rInv;
    BigIntMaxCap<C> r;
    details::setZero(r, s);
    r[s] = 1;
    r.size_ = s + 1;
    rInv.size_ = s + 1;

    auto mCopy = n;
    mCopy[s] = 0;
    uint64_t tmpStorage[C * 11];
    bigIntKernels::modInv(rInv.ptr(), r.ptr(), mCopy.ptr(), s + 1, tmpStorage);

    BigIntMaxCap<2 * C> r_rInv;
    details::setZero(r_rInv, s);
    bigIntKernels::copy(r_rInv.ptr() + s, rInv.ptr(), s + 1, s + 1);
    r_rInv.size_ = 2 * s + 1;

    BigIntMaxCap<C> one = 1;
    sub(r_rInv, r_rInv, one);
    div(result.k, r_rInv, n);

    return result;
}
template<int C> void convertToMontgomeryForm(BigIntMaxCap<C>& r, const BigIntMaxCap<C>& a, const MontgomeryReductionMod<BigIntMaxCap<C>>& m) {
    mod(r, a << m.b, m.mod);
}
template<int C> BigIntMaxCap<C> getValueInMontgomeryForm(const BigIntMaxCap<C>& a, const MontgomeryReductionMod<BigIntMaxCap<C>>& m) {
    BigIntMaxCap<C> r;
    convertToMontgomeryForm(r, a, m);
    return r;
}
template<int C1, int C2, int C3> void mod(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const MontgomeryReductionMod<BigIntMaxCap<C3>>& m) {
    mul(m.t, a, m.k);
    m.s = m.t;
    mul(m.t, m.s, m.mod);
    add(m.t, m.t, a);
    bigIntKernels::copy(r.ptr(), m.t.ptr() + m.b, m.b, m.b);
    r.size_ = m.t.size() - m.b;
    if (r >= m.mod) {
        sub(r, r, m.mod);
    }

    //Limb s[S];
    //Limb t[2 * S];

    //mul<2 * S, S, S>(t, x, k);
    //copy<S>(s, t);
    //mul<2 * S, S, S>(t, s, m);
    //add<2 * S>(t, t, x);
    //copy<S>(u, t + S);
    //if (cmp<S>(u, m) >= 0) {
    //    sub<S>(u, u, m);
    //}
}

template<int C1, int C2, int C3, int C4> void modAdd(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const MontgomeryReductionMod<BigIntMaxCap<C4>>& m) {
    modAdd(r, a, b, m.mod);
}
template<int C1, int C2, int C3, int C4> void modSub(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const MontgomeryReductionMod<BigIntMaxCap<C4>>& m) {
    modSub(r, a, b, m.mod);
}
template<int C1, int C2, int C3, int C4> void modMul(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const BigIntMaxCap<C3>& b, const MontgomeryReductionMod<BigIntMaxCap<C4>>& m) {
    BigIntMaxCap<C1 * 2> tmp;
    mul(tmp, a, b);
    mod(r, tmp, m);
}
template<int C1, int C2, int C3> void modSqr(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const MontgomeryReductionMod<BigIntMaxCap<C3>>& m) {
    BigIntMaxCap<C1 * 2> tmp;
    sqr(tmp, a);
    mod(r, tmp, m);
}
template<int C1, int C2, int C3> void modNeg(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const MontgomeryReductionMod<BigIntMaxCap<C3>>& m) {
    modNeg(r, a, m.mod);
}
template<int C1, int C2, int C3> void modDbl(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const MontgomeryReductionMod<BigIntMaxCap<C3>>& m) {
    modDbl(r, a, m.mod);
}


template<int C> BigIntMaxCap<C> getModValue(const BigIntMaxCap<C>& m) {
    return m;
}
template<int C> BigIntMaxCap<C> getModValue(const BarretReductionMod<BigIntMaxCap<C>>& m) {
    return m.mod;
}
template<int C> BigIntMaxCap<C> getModValue(const MontgomeryReductionMod<BigIntMaxCap<C>>& m) {
    return m.mod;
}
template<int C> BigIntMaxCap<C> convertToValue(const BigIntMaxCap<C>& value, const BigIntMaxCap<C>& dummy_mod) {
    return value;
}
template<int C> BigIntMaxCap<C> convertToValue(const BigIntMaxCap<C>& value, const BarretReductionMod<BigIntMaxCap<C>>& dummy_mod) {
    return value;
}
template<int C> BigIntMaxCap<C> convertToValue(const BigIntMaxCap<C>& value, const MontgomeryReductionMod<BigIntMaxCap<C>>& mod) {
    return value % mod;
}
template<int C> BigIntMaxCap<C> getConstant(uint64_t a, const BigIntMaxCap<C>& m) {
    return BigIntMaxCap<C> { a } % m;
}
template<int C> BigIntMaxCap<C> getConstant(uint64_t a, const BarretReductionMod<BigIntMaxCap<C>>& m) {
    return BigIntMaxCap<C> { a } % m.mod;
}
template<int C> BigIntMaxCap<C> getConstant(uint64_t a, const MontgomeryReductionMod<BigIntMaxCap<C>>& m) {
    return getValueInMontgomeryForm(BigIntMaxCap<C> { a } % m.mod, m);
}


/*
    Utility
*/
template<int C> bool isZero(const BigIntMaxCap<C>& a) {
    return a.size() == 1 && a[0] == 0;
}
template<int C> std::string toString(const BigIntMaxCap<C>& a) {
    if (isZero(a)) {
        return "0";
    }
    auto aCopy = a;
    aCopy.size_ = abs(aCopy.size_);
    char buf[C*22] = { 0 };
    int i = 0;
    while (!isZero(aCopy)) {
        buf[i++] = (char)('0' + (aCopy % 10));
        aCopy = aCopy / 10;
    }
    std::reverse(buf, buf + i);
    return (a.size() < 0 ? "-" : "") + std::string(buf);
}
template<int C> std::ostream& operator<<(std::ostream& out, const BigIntMaxCap<C>& value) {
    return out << toString(value);
}
BigIntMaxCap<15> operator ""_cap(const char* str) {
    return BigIntMaxCap<15>(std::string(str));
}

template<int C1, int C2> void mul(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint64_t b) {
    if (&r == &a) {
        auto aCopy = a;
        bigIntKernels::mulLimb(r.ptr(), aCopy.ptr(), b, a.size());
    } else {
        bigIntKernels::mulLimb(r.ptr(), a.ptr(), b, a.size());
    }
    r.size_ = a.size() + (r[a.size()] != 0);
}
template<int C> void mulAssign(BigIntMaxCap<C>& r, const uint64_t b) {
    mul(r, r, b);
}
template<int C1, int C2> void add(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint64_t b) {
    auto carry = bigIntKernels::add(r.ptr(), a.ptr(), &b, a.size(), 1);
    if (carry) {
        r.size_ = a.size()+1;
        r[r.size_] = carry;
    } else {
        r.size_ = a.size();
    }
}
template<int C1, int C2> void div(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint64_t b) {
    details::setZero(r, a.size());
    BigIntMaxCap<C2> bCopy = b;
    div(r, a, bCopy);
    r.size_ = a.size() - ((a.size() > 1) && (r[a.size() - 1] == 0));
}
template<int C1, int C2> void mod(BigIntMaxCap<C1>& r, const BigIntMaxCap<C2>& a, const uint64_t b) {
    BigIntMaxCap<C2> bCopy = b;
    mod(r, a, bCopy);
}


// -------------------------------------
// Operators
// -------------------------------------

// TODO: probably hand-written version for each comparison is better
template<int C1, int C2> bool operator>(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) > 0;
}
template<int C1, int C2> bool operator<(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) < 0;
}
template<int C1, int C2> bool operator>=(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) >= 0;
}
template<int C1, int C2> bool operator<=(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) <= 0;
}
template<int C1, int C2> bool operator==(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) == 0;
}
template<int C1, int C2> bool operator!=(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    return cmp(a, b) != 0;
}
template<int C1, int C2> bool operator<(const BigIntMaxCap<C1>& a, const BarretReductionMod<BigIntMaxCap<C2>>& m) {
    return a < m.mod;
}
template<int C> bool operator>(const BigIntMaxCap<C>& a, int b) {
    return cmp(a, BigIntMaxCap<C>((uint64_t)b)) > 0;
}

template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator+(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; add(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator-(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; sub(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator*(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; mul(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator/(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; div(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator%(const BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; mod(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator%(const BigIntMaxCap<C1>& a, const BarretReductionMod<BigIntMaxCap<C2>>& b) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; mod(r, a, b); return r;
}
template<int C1, int C2> BigIntMaxCap<(C1 > C2) ? C1 : C2> operator%(const BigIntMaxCap<C1>& a, const MontgomeryReductionMod<BigIntMaxCap<C2>>& m) {
    BigIntMaxCap<(C1 > C2) ? C1 : C2> r; mod(r, a, m); return r;
}

template<int C> BigIntMaxCap<C> operator+(const BigIntMaxCap<C>& a, const uint64_t b) {
    BigIntMaxCap<C> r; add(r, a, b); return r;
}
//template<int C> BigIntMaxCap<C> operator-(const BigIntMaxCap<C>& a, const uint64_t b) {
//    BigIntMaxCap<C> r; sub(r, a, b); return r;
//}
template<int C> BigIntMaxCap<C> operator*(const BigIntMaxCap<C>& a, const uint64_t b) {
    BigIntMaxCap<C> r; mul(r, a, b); return r;
}
template<int C> BigIntMaxCap<C> operator/(const BigIntMaxCap<C>& a, const uint64_t b) {
    BigIntMaxCap<C> r; div(r, a, b); return r;
}
template<int C> uint64_t operator%(const BigIntMaxCap<C>& a, const uint64_t b) {
    BigIntMaxCap<C> r; mod(r, a, b); return r[0];
}
template<int C> BigIntMaxCap<C> operator<<(const BigIntMaxCap<C>& a, const uint32_t b) {
    BigIntMaxCap<C> r; shl(r, a, b); return r;
}
template<int C> BigIntMaxCap<C> operator>>(const BigIntMaxCap<C>& a, const uint32_t b) {
    BigIntMaxCap<C> r; shr(r, a, b); return r;
}
template<int C> BigIntMaxCap<C> sqr(const BigIntMaxCap<C>& a) {
    BigIntMaxCap<C> r; sqr(r, a); return r;
}

template<int C1, int C2> BigIntMaxCap<C1>& operator+=(BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    add(a, a, b); return a;
}
template<int C1, int C2> BigIntMaxCap<C1>& operator-=(BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    sub(a, a, b); return a;
}
template<int C1, int C2> BigIntMaxCap<C1>& operator*=(BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    mul(a, a, b); return a;
}
template<int C1, int C2> BigIntMaxCap<C1>& operator/=(BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    div(a, a, b); return a;
}
template<int C1, int C2> BigIntMaxCap<C1>& operator%=(BigIntMaxCap<C1>& a, const BigIntMaxCap<C2>& b) {
    mod(a, a, b); return a;
}

template<int C> BigIntMaxCap<C>& operator+=(BigIntMaxCap<C>& a, const uint64_t b) {
    add(a, a, b); return a;
}
//template<int C> BigIntMaxCap<C>& operator-=(BigIntMaxCap<C>& a, const uint64_t b) {
//    sub(a, a, b); return a;
//}
template<int C> BigIntMaxCap<C>& operator*=(BigIntMaxCap<C>& a, const uint64_t b) {
    mul(a, a, b); return a;
}
template<int C> BigIntMaxCap<C>& operator/=(BigIntMaxCap<C>& a, const uint64_t b) {
    div(a, a, b); return a;
}
template<int C> BigIntMaxCap<C>& operator<<=(BigIntMaxCap<C>& a, const uint32_t b) {
    shl(a, a, b); return a;
}
template<int C> BigIntMaxCap<C>& operator>>=(BigIntMaxCap<C>& a, const uint32_t b) {
    shr(a, a, b); return a;
}
