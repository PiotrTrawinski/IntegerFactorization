#pragma once
#include "../Utility/compilerMacros.h"
#include "../Utility/debugAssert.h"
#include "../Utility/bitManipulation.h"
#include "kernels.h"
#include "common.h"
#include "Unsigned64.h"
#include <iostream>
#include <string>
#include <cstring>
#include <array>
#include <algorithm>


/*
    Big Integer class with fixed size.
    All operations are always done on the whole size, so bigger size means more operations even if you store only small numbers.
    It's the user's responsibility to not overflow the size.

    It stores natural + {0} numbers only. 
    If the result of subtraction is negative it stores the absolute value and returns information about the sign.
*/
template<int Size=4> struct BigIntFixedSize { 
    // data:
    std::array<uint64_t, Size> data;
    
    // methods:
    BigIntFixedSize(uint64_t value) { data[0] = value; bigIntKernels::clear<Size - 1>(ptr()+1); }
    BigIntFixedSize() : BigIntFixedSize(0) {}
    BigIntFixedSize(const std::string& str) { parseString(str); }

    template<int OtherSize> BigIntFixedSize(const BigIntFixedSize<OtherSize>& other) {
        bigIntKernels::copy<Size, OtherSize>(ptr(), other.ptr());
    }
    template<int OtherSize> BigIntFixedSize& operator=(const BigIntFixedSize<OtherSize>& other) {
        bigIntKernels::copy<Size, OtherSize>(ptr(), other.ptr());
        return *this;
    }
    BigIntFixedSize& operator=(uint64_t value) {
        data[0] = value;
        bigIntKernels::clear<Size - 1>(ptr() + 1);
        return *this;
    }
    BigIntFixedSize& operator=(const std::string& str) {
        parseString(str);
        return *this;
    }

    uint64_t* ptr()             { return &data[0]; }
    const uint64_t* ptr() const { return &data[0]; }

    uint64_t& operator[](int index) {
        debugAssert(index < Size);
        return data[index];
    }
    const uint64_t& operator[](int index) const {
        debugAssert(index < Size);
        return data[index]; 
    }

    BigIntBit bit(int index)                  { return BigIntBit(ptr() + index / 64, index % 64); }
    const BigIntConstBit bit(int index) const { return BigIntConstBit(ptr() + index / 64, index % 64); }

    int32_t size() const { 
        return Size;
    }
    uint32_t sizeInBits() const {
        auto realS = realSize();
        return (realS-1) * 64 + static_cast<uint32_t>(::sizeInBits(data[realS - 1]));
    }
    BigIntFixedSize<Size> mostSignificantBit() const {
        auto realS = realSize();
        BigIntFixedSize<Size> res;
        res.data[realS - 1] = ::mostSignificantBit(data[realS - 1]);
        return res;
    }

    int32_t realSize() const {
        for (int i = Size - 1; i >= 0; --i) {
            if (data[i] != 0) {
                return i + 1;
            }
        }
        return 1;
    }

private:
    void parseString(const std::string& str) {
        bigIntKernels::clear<Size>(ptr());
        for (auto c : str) {
            mul(*this, *this, 10);
            add(*this, *this, static_cast<uint64_t>(c) - '0');
        }
    }
};


template<int Size> struct BigIntParseImpl<BigIntFixedSize<Size>> {
    static BigIntFixedSize<Size> parse(const std::string& str) {
        return BigIntFixedSize<Size> { str };
    }
};

template<int Size> int sizeInLimbs(const BigIntFixedSize<Size>& a) {
    return Size;
}
template<int Size> int sizeInLimbs(const BarretReductionMod<BigIntFixedSize<Size>>& a) {
    return Size;
}
template<int Size> int sizeInLimbs(const MontgomeryReductionMod<BigIntFixedSize<Size>>& a) {
    return Size;
}
template<int S1, int S2> int cmp(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return bigIntKernels::cmp<S1, S2>(a.ptr(), b.ptr());
}
template<int S> void swap(BigIntFixedSize<S>& a, BigIntFixedSize<S>& b) {
    std::swap(a.data, b.data);
}
template<int Size> void assign(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a) {
    r = a;
}
template<int Size> void assign(BigIntFixedSize<Size>& r, const int& a) {
    r = a;
}
template<int Size> BigIntFixedSize<Size> getModValue(const BigIntFixedSize<Size>& m) {
    return m;
}
template<int Size> BigIntFixedSize<Size> getModValue(const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    return m.mod;
}
template<int Size> BigIntFixedSize<Size> getModValue(const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    return m.mod;
}
template<int Size> BigIntFixedSize<Size> convertToValue(const BigIntFixedSize<Size>& value, const BigIntFixedSize<Size>& dummy_mod) {
    return value;
}
template<int Size> BigIntFixedSize<Size> convertToValue(const BigIntFixedSize<Size>& value, const BarretReductionMod<BigIntFixedSize<Size>>& dummy_mod) {
    return value;
}
template<int Size> BigIntFixedSize<Size> convertToValue(const BigIntFixedSize<Size>& value, const MontgomeryReductionMod<BigIntFixedSize<Size>>& mod) {
    return value % mod;
}
template<int Size> BigIntFixedSize<Size> getConstant(uint64_t a, const BigIntFixedSize<Size>& m) {
    return BigIntFixedSize<Size> { a } % m;
}
template<int Size> BigIntFixedSize<Size> getConstant(uint64_t a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    return BigIntFixedSize<Size> { a } % m.mod;
}
template<int Size> BigIntFixedSize<Size> getConstant(uint64_t a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    return getValueInMontgomeryForm(BigIntFixedSize<Size> { a } % m.mod, m);
}
template<int Size> void add(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    bigIntKernels::add<Size>(r.ptr(), a.ptr(), b.ptr());
}
template<int Size> bool sub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    return bigIntKernels::sub<Size>(r.ptr(), a.ptr(), b.ptr());
}
template<int Size> void absSub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    sub(r, a, b);
}
template<int S1, int S2> void mul(BigIntFixedSize<S1+S2>& r, const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    bigIntKernels::mul<S1, S2>(r.ptr(), a.ptr(), b.ptr());
}
template<int S> void mul(BigIntFixedSize<2*S>& r, const BigIntFixedSize<S>& a, const BigIntFixedSize<S>& b) {
    bigIntKernels::mul<S>(r.ptr(), a.ptr(), b.ptr());
}
template<int S> void mul(BigIntFixedSize<S>& r, const BigIntFixedSize<S>& a, const BigIntFixedSize<S>& b) {
    BigIntFixedSize<2 * S> dummyRes;
    mul<S, S>(dummyRes, a, b);
    r = dummyRes;
}
template<int S> void sqr(BigIntFixedSize<S*2>& r, const BigIntFixedSize<S>& a) {
    bigIntKernels::sqr<S>(r.ptr(), a.ptr());
}
template<int S> void sqr(BigIntFixedSize<S>& r, const BigIntFixedSize<S>& a) {
    BigIntFixedSize<2 * S> res; sqr(res, a); r = res;
}
template<int S1, int S2, int S3> void div(BigIntFixedSize<S1>& r, const BigIntFixedSize<S2>& a, const BigIntFixedSize<S3>& b) {
    bigIntKernels::div<S1, S2, S3>(r.ptr(), a.ptr(), b.ptr());
}
template<int SizeMod, int Size> void mod(BigIntFixedSize<SizeMod>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<SizeMod>& b) {
    bigIntKernels::mod<SizeMod, Size, SizeMod>(r.ptr(), a.ptr(), b.ptr());
}
template<int Size> void shl(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const uint32_t c) {
    bigIntKernels::shl<Size, Size>(r.ptr(), a.ptr(), c);
}
template<int Size> void shr(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const uint32_t c) {
    bigIntKernels::shr<Size, Size>(r.ptr(), a.ptr(), c);
}


template<int Size> void modAdd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BigIntFixedSize<Size>& m) {
    if constexpr (Size == 1) {
        modAdd(r[0], a[0], b[0], m[0]);
        return;
    }
    add(r, a, b);
    if (r >= m) {
        sub(r, r, m);
    }
}
template<int Size> void modSub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BigIntFixedSize<Size>& m) {
    if constexpr (Size == 1) {
        modSub(r[0], a[0], b[0], m[0]);
        return;
    }
    auto isNegative = sub(r, a, b);
    if (isNegative) {
        sub(r, m, r);
    }
}
template<int Size> void modMul(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BigIntFixedSize<Size>& m) {
    uint64_t tmp[2 * Size];
    bigIntKernels::mul<Size>(tmp, a.ptr(), b.ptr());
    bigIntKernels::mod<Size, 2 * Size, Size>(r.ptr(), tmp, m.ptr());
}
template<int Size> void modSqr(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& m) {
    uint64_t tmp[2 * Size];
    bigIntKernels::sqr<Size>(tmp, a.ptr());
    bigIntKernels::mod<Size, 2 * Size, Size>(r.ptr(), tmp, m.ptr());
}
template<int Size> void modNeg(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& m) {
    if constexpr (Size == 1) {
        modNeg(r[0], a[0], m[0]);
        return;
    }
    sub(r, m, a);
}
template<int Size> void modDbl(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& m) {
    if constexpr (Size == 1) {
        modDbl(r[0], a[0], m[0]);
        return;
    }
    shl(r, a, 1);
    if (r >= m) {
        sub(r, r, m);
    }
}


template<int Size> MontgomeryReductionMod<BigIntFixedSize<Size>> getMontgomeryReductionMod(const BigIntFixedSize<Size>& n) {
    MontgomeryReductionMod<BigIntFixedSize<Size>> result;
    bigIntKernels::createMontgomeryReductionMod<Size>(result.k.ptr(), n.ptr(), result.b);
    result.mod = n;
    return result;
}
template<int Size> void convertToMontgomeryForm(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::convertToMontgomeryForm<Size>(r.ptr(), a.ptr(), m.mod.ptr(), m.b);
}
template<int Size> BigIntFixedSize<Size> getValueInMontgomeryForm(const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    BigIntFixedSize<Size> r; 
    convertToMontgomeryForm(r, a, m);
    return r;
}
template<int Size> void mod(BigIntFixedSize<Size>& r, const BigIntFixedSize<2*Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::montgomeryReduction<Size>(r.ptr(), a.ptr(), m.k.ptr(), m.mod.ptr(), m.b);
}
template<int Size> void mod(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    BigIntFixedSize<2 * Size> aa = a;
    bigIntKernels::montgomeryReduction<Size>(r.ptr(), aa.ptr(), m.k.ptr(), m.mod.ptr(), m.b);
}
template<int Size> void modAdd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    modAdd(r, a, b, m.mod);
}
template<int Size> void modSub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    modSub(r, a, b, m.mod);
}
template<int Size> void modMul(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::montgomeryMult<Size>(r.ptr(), a.ptr(), b.ptr(), m.k.ptr(), m.mod.ptr(), m.b);
}
template<int Size> void modSqr(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::montgomerySqr<Size>(r.ptr(), a.ptr(), m.k.ptr(), m.mod.ptr(), m.b);
}
template<int Size> void modNeg(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    modNeg(r, a, m.mod);
}
template<int Size> void modDbl(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    modDbl(r, a, m.mod);
}


template<int Size> BarretReductionMod<BigIntFixedSize<Size>> getBarretReductionMod(const BigIntFixedSize<Size>& n) {
    BarretReductionMod<BigIntFixedSize<Size>> result;
    result.mod = n;
    result.k = bigIntKernels::sizeInBits<Size>(n.ptr()) * 2;

    BigIntFixedSize<2*Size> powK = 1;
    shl(powK, powK, result.k);
    result.R = 1;
    div(result.R, powK, n);

    return result;
}
template<int Size> void mod(BigIntFixedSize<Size>& r, const BigIntFixedSize<2*Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& b) {
    bigIntKernels::modBarret<Size>(r.ptr(), a.ptr(), b.mod.ptr(), b.R.ptr(), b.k);
}
template<int Size> void modAdd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    modAdd(r, a, b, m.mod);
}
template<int Size> void modSub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    modSub(r, a, b, m.mod);
}
template<int Size> void modMul(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    uint64_t tmp[2 * Size];
    bigIntKernels::mul<Size>(tmp, a.ptr(), b.ptr());
    bigIntKernels::modBarret<Size>(r.ptr(), tmp, m.mod.ptr(), m.R.ptr(), m.k);
}
template<int Size> void modSqr(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    uint64_t tmp[2 * Size];
    bigIntKernels::sqr<Size>(tmp, a.ptr());
    bigIntKernels::modBarret<Size>(r.ptr(), tmp, m.mod.ptr(), m.R.ptr(), m.k);
}
template<int Size> void modNeg(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    if constexpr (Size == 1) {
        modNeg(r, a, m.mod);
        return;
    }
    sub(r, m.mod, a);
}
template<int Size> void modDbl(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    if constexpr (Size == 1) {
        modDbl(r, a, m.mod);
        return;
    }
    add(r, a, a);
    if (r >= m.mod) {
        sub(r, r, m.mod);
    }
}
template<int Size> void modInv(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& m) {
    bigIntKernels::modInv<Size>(r.ptr(), a.ptr(), m.ptr());
}
template<int Size> void modInv(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::modInv<Size>(r.ptr(), a.ptr(), m.mod.ptr()); // TODO: is it ok? probably not
}
template<int Size> void modInv(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::modInv<Size>(r.ptr(), a.ptr(), m.mod.ptr()); // TODO: is it ok? probably not
}
template<int Size> void gcd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    bigIntKernels::gcd<Size>(r.ptr(), a.ptr(), b.ptr());
}
template<int Size> void gcd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const BarretReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::gcd<Size>(r.ptr(), a.ptr(), m.mod.ptr());
}
template<int Size> void gcd(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const MontgomeryReductionMod<BigIntFixedSize<Size>>& m) {
    bigIntKernels::gcd<Size>(r.ptr(), a.ptr(), m.mod.ptr());
}


template<int Size> void add(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, uint64_t b) {
    bigIntKernels::add<Size, Size, 1>(r.ptr(), a.ptr(), &b);
}
template<int Size> void sub(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, uint64_t b) {
    bigIntKernels::sub<Size, Size, 1>(r.ptr(), a.ptr(), &b);
}
template<int Size> void mul(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, uint64_t b) {
    bigIntKernels::mul<Size, Size, 1>(r.ptr(), a.ptr(), &b);
}

template<int Size> void div(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const uint64_t b) {
    bigIntKernels::div<Size, Size, 1>(r.ptr(), a.ptr(), &b);
}
template<int Size> void mod(BigIntFixedSize<Size>& r, const BigIntFixedSize<Size>& a, const uint64_t b) {
    bigIntKernels::mod<Size, Size, 1>(r.ptr(), a.ptr(), &b);
}
template<int Size> uint32_t mod(const BigIntFixedSize<Size>& a, uint32_t m, int tableEntryId) {
    return bigIntKernels::mod<Size>(a.ptr(), m, tableEntryId);
}

template<int Size> bool isZero(const BigIntFixedSize<Size>& a) {
    for (int i = 0; i < Size; ++i) {
        if (a[i] != 0)
            return false;
    }
    return true;
}
template<int Size> std::string toString(const BigIntFixedSize<Size>& a) {
    if (isZero(a)) {
        return "0";
    }
    auto aCopy = a;
    char buf[Size*22] = { 0 };
    int i = 0;
    while (!isZero(aCopy)) {
        buf[i++] = (char)('0' + (aCopy % 10));
        aCopy = aCopy / 10;
    }
    std::reverse(buf, buf + i);
    return std::string(buf);
}
template<int Size> std::ostream& operator<<(std::ostream& out, const BigIntFixedSize<Size>& a) {
    return out << toString(a);
}


// TODO: probably hand-written version for each comparison is better
template<int S1, int S2> bool operator>(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) > 0;
}
template<int S1, int S2> bool operator<(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) < 0;
}
template<int S1, int S2> bool operator>=(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) >= 0;
}
template<int S1, int S2> bool operator<=(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) <= 0;
}
template<int S1, int S2> bool operator==(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) == 0;
}
template<int S> bool operator==(const BigIntFixedSize<S>& a, int b) {
    if constexpr (S == 1) {
        return a[0] == b;
    } else {
        return a[0] == b && a.realSize() == 1;
    }
}
template<int S1, int S2> bool operator==(const BarretReductionMod<BigIntFixedSize<S1>>& a, const BarretReductionMod<BigIntFixedSize<S2>>& b) {
    return cmp(a.mod, b.mod) == 0;
}
template<int S1, int S2> bool operator==(const MontgomeryReductionMod<BigIntFixedSize<S1>>& a, const MontgomeryReductionMod<BigIntFixedSize<S2>>& b) {
    return cmp(a.mod, b.mod) == 0;
}
template<int S1, int S2> bool operator!=(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    return cmp(a, b) != 0;
}
template<int S1, int S2> bool operator!=(const BigIntFixedSize<S1>& a, const BarretReductionMod<BigIntFixedSize<S2>>& m) {
    return a != m.mod;
}
template<int S1, int S2> bool operator!=(const BigIntFixedSize<S1>& a, const MontgomeryReductionMod<BigIntFixedSize<S2>>& m) {
    return a != m.mod;
}
template<int S> bool operator>(const BigIntFixedSize<S>& a, int b) {
    uint64_t bb = b;
    return bigIntKernels::cmp<S, 1>(a.ptr(), &bb) > 0;
}
template<int S> bool operator<=(const BigIntFixedSize<S>& a, uint64_t b) {
    return bigIntKernels::cmp<S, 1>(a.ptr(), &b) <= 0;
}
template<int S> bool operator>=(const BigIntFixedSize<S>& a, uint64_t b) {
    return bigIntKernels::cmp<S, 1>(a.ptr(), &b) >= 0;
}
template<int S> BigIntFixedSize<S>& operator+=(BigIntFixedSize<S>& r, int a) {
    add(r, r, a);
    return r;
}
template<int S> BigIntFixedSize<S>& operator-=(BigIntFixedSize<S>& r, int a) {
    if (a >= 0) {
        sub(r, r, a);
    } else {
        add(r, r, -a);
    }
    return r;
}
template<int S> BigIntFixedSize<S>& operator*=(BigIntFixedSize<S>& r, uint64_t a) {
    mul(r, r, a);
    return r;
}
template<int S1, int S2> BigIntFixedSize<S1>& operator*=(BigIntFixedSize<S1>& r, const BigIntFixedSize<S2>& b) {
    BigIntFixedSize<S1 + S2> res;
    mul(res, r, b);
    r = res;
    return r;
}

template<int Size> BigIntFixedSize<Size> operator+(const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    BigIntFixedSize<Size> r; add(r, a, b); return r;
}
template<int Size> BigIntFixedSize<Size> operator-(const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    BigIntFixedSize<Size> r; sub(r, a, b); return r;
}
template<int Size> BigIntFixedSize<2*Size> operator*(const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    BigIntFixedSize<2*Size> r; mul<Size, Size>(r, a, b); return r;
}
template<int Size> BigIntFixedSize<Size> operator/(const BigIntFixedSize<Size>& a, const BigIntFixedSize<Size>& b) {
    BigIntFixedSize<Size> r; div(r, a, b); return r;
}
template<int S1, int S2> BigIntFixedSize<S2> operator%(const BigIntFixedSize<S1>& a, const BigIntFixedSize<S2>& b) {
    BigIntFixedSize<S2> r; mod(r, a, b); return r;
}
template<int S1, int S2> BigIntFixedSize<S2> operator%(const BigIntFixedSize<S1>& a, const BarretReductionMod<BigIntFixedSize<S2>>& b) {
    BigIntFixedSize<S2> r; mod(r, a, b); return r;
}
template<int S1, int S2> BigIntFixedSize<S2> operator%(const BigIntFixedSize<S1>& a, const MontgomeryReductionMod<BigIntFixedSize<S2>>& b) {
    BigIntFixedSize<S2> r; mod(r, a, b); return r;
}
template<int Size> BigIntFixedSize<2*Size> sqr(const BigIntFixedSize<Size>& a) {
    BigIntFixedSize<2*Size> r; sqr(r, a); return r;
}
template<int Size> BigIntFixedSize<Size> operator/(const BigIntFixedSize<Size>& a, const uint64_t b) {
    BigIntFixedSize<Size> r; div(r, a, b); return r;
}
template<int Size> uint64_t operator%(const BigIntFixedSize<Size>& a, const uint64_t b) {
    BigIntFixedSize<Size> r; mod(r, a, b); return r[0];
}
template<int Size> BigIntFixedSize<Size> operator<<(const BigIntFixedSize<Size>& a, const uint32_t c) {
    BigIntFixedSize<Size> r; shl(r, a, c); return r;
}
template<int Size> BigIntFixedSize<Size> operator>>(const BigIntFixedSize<Size>& a, const uint32_t c) {
    BigIntFixedSize<Size> r; shr(r, a, c); return r;
}

namespace std {
    template<int Size> struct hash<BigIntFixedSize<Size>> {
        std::size_t operator()(const BigIntFixedSize<Size>& k) const {
            return hash<uint64_t>()(k[0]);
        }
    };
    template<int Size> struct hash<BarretReductionMod<BigIntFixedSize<Size>>> {
        std::size_t operator()(const BarretReductionMod<BigIntFixedSize<Size>>& k) const {
            return hash<uint64_t>()(k.mod[0]);
        }
    };
    template<int Size> struct hash<MontgomeryReductionMod<BigIntFixedSize<Size>>> {
        std::size_t operator()(const MontgomeryReductionMod<BigIntFixedSize<Size>>& k) const {
            return hash<uint64_t>()(k.mod[0]);
        }
    };
}
