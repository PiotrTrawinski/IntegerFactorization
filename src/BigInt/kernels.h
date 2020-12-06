#pragma once
#include "../Utility/bitManipulation.h"
#include "../Utility/debugAssert.h"
#include "../Utility/alwaysInline.h"
#include "../PrecomputedTables/modInverseTable.h"
#include "common.h"
#include "64bitIntrinsics.h"
#include <cstdint>
#include <gmp.h>

namespace bigIntKernels {
    using Limb     = uint64_t;
    using Int      = Limb*;
    using ConstInt = const Limb*;


    template<int A, int B> constexpr int min() {
        return A < B ? A : B;
    }
    template<int Size> int32_t nonZeroLimbCount(ConstInt a) {
        for (int i = Size - 1; i >= 0; --i) {
            if (a[i] != 0) {
                return i + 1;
            }
        }
        return 1;
    }
    template<int Size> int32_t sizeInBits(ConstInt a) {
        auto limbCount = nonZeroLimbCount<Size>(a);
        return (limbCount - 1) * 64 + static_cast<uint32_t>(::sizeInBits(a[limbCount - 1]));
    }
    template<int Size> void clear(Int a) {
        for (int i = 0; i < Size; ++i) {
            a[i] = 0;
        }
    }
    template<int Size> int32_t realSize(ConstInt a) {
        for (int i = Size - 1; i >= 0; --i) {
            if (a[i] != 0) {
                return i + 1;
            }
        }
        return 1;
    }


    // r = a
    // if size(a) < size(r) the higher limbs are lost
    // if size(r) > size(a) it fills remaining limbs with 0's
    template<int SR, int SA> void copy(Int r, ConstInt a) {
        memcpy(r, a, min<SR, SA>() * sizeof(Limb));
        if constexpr (SR > SA) {
            memset(r + SA, 0, (SR - SA) * sizeof(Limb));
        }
    }
    template<int S> void copy(Int r, ConstInt a) {
        copy<S, S>(r, a);
    }
    void copy(Int r, ConstInt a, int32_t sr, int32_t sa) {
        if (sr <= sa) {
            memcpy(r, a, sr * sizeof(Limb));
        } else {
            memcpy(r, a, sa * sizeof(Limb));
            memset(r + sa, 0, (sr - sa) * sizeof(Limb));
        }
    }


    // Compares values of 'a' and 'b'
    // if a < b: return negative
    // if a > b: return positive
    // if a = b: return 0
    template<int SA, int SB> int cmp(ConstInt a, ConstInt b) {
        if constexpr (SA == 1 && SB == 1) {
            if (a[0] < b[0]) return -1;
            else             return a[0] > b[0];
        } else if constexpr (SA < SB) {
            return -cmp<SB, SA>(b, a);
        } else {
            for (int i = SA - 1; i >= SB; --i) {
                if (a[i] > 0) {
                    return 1;
                }
            }
            for (int i = SB - 1; i >= 0; --i) {
                if      (a[i] < b[i]) return -1;
                else if (a[i] > b[i]) return 1;
            }
            return 0;
        }
    }
    template<int S> int cmp(ConstInt a, ConstInt b) {
        return cmp<S, S>(a, b);
    }
    int cmp(ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        if (sa < sb) {
            return -cmp(b, a, sb, sa);
        }
        for (int i = sa - 1; i >= sb; --i) {
            if (a[i] > 0) {
                return 1;
            }
        }
        for (int i = sb - 1; i >= 0; --i) {
            if (a[i] < b[i]) return -1;
            else if (a[i] > b[i]) return 1;
        }
        return 0;
    }


    // r = a + b
    template<int SR, int SA, int SB> void add(Int r, ConstInt a, ConstInt b) {
        static_assert(SR >= SA && SR >= SB, "result needs to be at least as big as the operands");
        if constexpr (SR == 1) {
            r[0] = a[0] + b[0];
        } else if constexpr (SB > SA) {
            add<SR, SB, SA>(r, b, a);
        } else {
            auto carry = addCarry(0, r[0], a[0], b[0]);
            for (int i = 1; i < SB; ++i) {
                carry = addCarry(carry, r[i], a[i], b[i]);
            }
            for (int i = SB; i < SA; ++i) {
                carry = addCarry(carry, r[i], a[i], 0);
            }
            if constexpr (SR > SA) {
                r[SA] = carry;
                clear<SR-SA-1>(&r[SA + 1]);
            }
        }
    }
    template<int S> void add(Int r, ConstInt a, ConstInt b) {
        add<S, S, S>(r, a, b);
    }
    Limb add(Int r, ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        if (sb > sa) {
            std::swap(a, b);
            std::swap(sa, sb);
        }
        auto carry = addCarry(0, r[0], a[0], b[0]);
        for (int i = 1; i < sb; ++i) {
            carry = addCarry(carry, r[i], a[i], b[i]);
        }
        for (int i = sb; i < sa; ++i) {
            carry = addCarry(carry, r[i], a[i], 0);
        }
        return carry;
    }


    // r = a - b, with a >= b
    template<int SR, int SA, int SB> void sub_firstOperandIsBigger(Int r, ConstInt a, ConstInt b) {
        debugAssert((cmp<SA, SB>)(a, b) >= 0);
        if constexpr (SR == 1) {
            r[0] = a[0] - b[0];
        } else {
            auto borrow = subBorrow(0, r[0], a[0], b[0]);
            for (int i = 1; i < min<SR,SB>(); ++i) {
                borrow = subBorrow(borrow, r[i], a[i], b[i]);
            }
            for (int i = SB; i < min<SR,SA>(); ++i) {
                borrow = subBorrow(borrow, r[i], a[i], 0);
            }
            if constexpr (SR > SA) {
                clear<SR-SA-1>(&r[SA+1]);
            }
        }
    }
    void sub_firstOperandIsBigger(Int r, ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        debugAssert(cmp(a, b, sa, sb) >= 0);
        auto borrow = subBorrow(0, r[0], a[0], b[0]);
        int i = 1;
        for (; i < sb; ++i) {
            borrow = subBorrow(borrow, r[i], a[i], b[i]);
        }
        for (; i < sa; ++i) {
            borrow = subBorrow(borrow, r[i], a[i], 0);
        }
    }


    // r = |a - b|
    // if a >= b returns false (result is not negative)
    // if a <  b returns true  (result is negative)
    template<int SR, int SA, int SB> bool sub(Int r, ConstInt a, ConstInt b) {
        auto isNegative = cmp<SA, SB>(a, b) < 0;
        if (isNegative) {
            std::swap(a, b);
        }
        sub_firstOperandIsBigger<SR, SA, SB>(r, a, b);
        return isNegative;
    }
    template<int S> bool sub(Int r, ConstInt a, ConstInt b) {
        return sub<S, S, S>(r, a, b);
    }
    bool sub(Int r, ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        auto isNegative = cmp(a, b, sa, sb) < 0;
        if (isNegative) {
            std::swap(a, b);
            std::swap(sa, sb);
        }
        sub_firstOperandIsBigger(r, a, b, sa, sb);
        return isNegative;
    }


    // r = a * b
    template<int SR, int SA> ALWAYS_INLINE void mulLimb(Int r, ConstInt a, Limb b) {
        mul128(r[1], r[0], a[0], b);
        for (int i = 1; i < min<SR - 1, SA>(); ++i) {
            Limb u0;
            mul128(r[i + 1], u0, a[i], b);
            r[i] += u0;
            r[i + 1] += (r[i] < u0);
        }
        if constexpr (SR == SA) {
            r[SR - 1] += a[SR - 1] * b;
        }
    }
    ALWAYS_INLINE void mulLimbNoLastLimb(Int r, ConstInt a, Limb b, int32_t s) {
        Limb u0, crec, c, p1, p0, r0;
        crec = 0;
        do {
            u0 = *a++;
            mul128(p1, p0, u0, b);
            r0 = p0 + crec;
            c = p0 > r0;
            crec = p1 + c;
            *r++ = r0;
        } while (--s != 0);
    }
    ALWAYS_INLINE void mulLimb(Int r, ConstInt a, Limb b, int32_t s) {
        Limb u0, crec, c, p1, p0, r0;
        crec = 0;
        do {
            u0 = *a++;
            mul128(p1, p0, u0, b);
            r0 = p0 + crec;
            c = p0 > r0;
            crec = p1 + c;
            *r++ = r0;
        } while (--s != 0);
        *r = crec;
    }


    // r += a * b
    // returns carry Limb
    template<int SA> Limb mulAddLimb(Int r, ConstInt a, Limb b) {
        Limb crec = 0;
        for (int i = 0; i < SA; ++i) {
            Limb u0, c, p1, p0, r0;
            u0 = *a++;
            mul128(p1, p0, u0, b);
            r0 = *r;
            p0 = r0 + p0;
            c = r0 > p0;
            p1 = p1 + c;
            r0 = p0 + crec;
            c = p0 > r0;
            crec = p1 + c;
            *r++ = r0;
        }
        return crec;
    }
    ALWAYS_INLINE Limb mulAddLimb(Int r, ConstInt a, Limb b, int32_t s) {
        Limb u0, crec, c, p1, p0, r0;
        crec = 0;
        do {
            u0 = *a++;
            mul128(p1, p0, u0, b);
            r0 = *r;
            p0 = r0 + p0;
            c = r0 > p0;
            p1 = p1 + c;
            r0 = p0 + crec;
            c = p0 > r0;
            crec = p1 + c;
            *r++ = r0;
        } while (--s != 0);

        return crec;
    }


    // r = a * b, where 'r' cannot share memory with either 'a' or 'b'
    // note: 'a' and 'b' can share memory, but in that case you should use 'sqr' method instead.
    template<int SR, int SA, int SB> void mul_SeperateMemoryLocation(Int r, ConstInt a, ConstInt b) {
        static_assert(SR > 1);
        static_assert(SA >= SB);
        debugAssert(r != a && r != b);
        //mulLimb<SR, SA>(r, a, b[0]);
        if constexpr (SR == SA) {
            mulLimbNoLastLimb(r, a, b[0], SA);
        } else {
            mulLimb(r, a, b[0], SA);
        }
        for (int i = 1; i < SB; ++i) {
            r[SA + i] = mulAddLimb(&r[i], a, b[i], SA);
        }
    }
    void mul_SeperateMemoryLocation(Int r, ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        debugAssert(sa >= sb);
        debugAssert(r != a && r != b);
        mulLimb(r, a, b[0], sa);
        for (int i = 1; i < sb; ++i) {
            r[sa + i] = mulAddLimb(&r[i], a, b[i], sa);
        }
    }

    // r = a * b, where 'r' cannot share memory. Only the lower S+1 Limbs of the result are calculated and stored
    template<int S> void mul_SeperateMemoryLocation_onlyLowLimbs_plus1(Int r, ConstInt a, ConstInt b) {
        debugAssert(r != a && r != b);
        mulLimb(r, a, b[0], S);
        for (int i = 1; i < S; ++i) {
            mulAddLimb(&r[i], a, b[i], S-i+1);
        }
    }
    void mul_SeperateMemoryLocation_onlyLowLimbs_plus1(Int r, ConstInt a, ConstInt b, int32_t s) {
        debugAssert(r != a && r != b);
        mulLimb(r, a, b[0], s);
        for (int i = 1; i < s; ++i) {
            mulAddLimb(&r[i], a, b[i], s-i+1);
        }
    }

    // r = a * b, where 'r' cannot share memory. Only the lower S Limbs of the result are calculated and stored
    template<int S> void mul_SeperateMemoryLocation_onlyLowLimbs(Int r, ConstInt a, ConstInt b) {
        debugAssert(r != a && r != b);
        if constexpr (S == 1) {
            r[0] = a[0] * b[0];
        } else {
            mulLimbNoLastLimb(r, a, b[0], S);
            for (int i = 1; i < S; ++i) {
                mulAddLimb(&r[i], a, b[i], S - i);
            }
        }
    }
    void mul_SeperateMemoryLocation_onlyLowLimbs(Int r, ConstInt a, ConstInt b, int32_t s) {
        debugAssert(r != a && r != b);
        if (s == 1) {
            r[0] = a[0] * b[0];
        } else {
            mulLimb(r, a, b[0], s);
            for (int i = 1; i < s; ++i) {
                mulAddLimb(&r[i], a, b[i], s - i);
            }
        }
    }

    template<int S> void muladd_SeperateMemoryLocation(Int r, ConstInt a, ConstInt b) {
        debugAssert(r != a && r != b);
        Limb c = 0;
        for (int i = 0; i < S; ++i) {
            auto res = mulAddLimb(&r[i], a, b[i], S);
            c = addCarry((uint8_t)c, r[S + i], r[S + i], res);
        }
    }

    // r = a * b
    template<int SR, int SA, int SB> void mul(Int r, ConstInt a, ConstInt b) {
        static_assert(SR >= SA && SR >= SB, "result needs to be at least as big as the operands");
        static_assert((SR == SA && SB == 1) || SR == SA + SB, "for now implementation only allows some size combinations");
        if constexpr (SR == 1) {
            r[0] = a[0] * b[0];
        } else if constexpr (SB > SA) {
            mul<SR, SB, SA>(r, b, a);
        } else {
            Limb aCopy[SA];
            Limb bCopy[SB];
            if (r == a) {
                copy<SA>(aCopy, a);
                a = aCopy;
            }
            if (r == b) {
                copy<SB>(bCopy, b);
                b = bCopy;
            }
            mul_SeperateMemoryLocation<SR, SA, SB>(r, a, b);
        }
    }
    template<int SA, int SB> void mul(Int r, ConstInt a, ConstInt b) {
        mul<SA+SB, SA, SB>(r, a, b);
    }
    template<int S> void mul(Int r, ConstInt a, ConstInt b) {
        mul<2*S, S, S>(r, a, b);
    }
    void mul(Int r, ConstInt a, ConstInt b, int32_t sa, int32_t sb) {
        debugAssert(r != a && r != b);
        if (sb > sa) {
            std::swap(a, b);
            std::swap(sa, sb);
        }
        mul_SeperateMemoryLocation(r, a, b, sa, sb);
    }


    // r = a >> c (right binary shift by 'c' bits, aka division by 2^c)
    template<int SR, int SA> void shr(Int r, ConstInt a, uint32_t c) {
        int partSkipCount = c / 64;
        int shiftSize = c - partSkipCount * 64;
        if (shiftSize == 0) {
            for (int i = 0; i < SA - partSkipCount; ++i) {
                r[i] = a[i + partSkipCount];
            }
        } else {
            for (int i = 0; i < SA - partSkipCount - 1; ++i) {
                r[i] = (a[i + partSkipCount] >> shiftSize) | (a[i + 1 + partSkipCount] << (64 - shiftSize));
            }
            if (SA - 1 - partSkipCount < SR) {
                r[SA - 1 - partSkipCount] = a[SA - 1] >> shiftSize;
            }
        }
        for (int i = SA - partSkipCount; i < SR; ++i) {
            r[i] = 0;
        }
    }

    // r = a << c (left binary shift by 'c' bits, aka multiplication by 2^c)
    template<int SR, int SA> void shl(Int r, ConstInt a, uint32_t c) {
        if constexpr (SR == 1) {
            r[0] = a[0] << c;
            return;
        }
        int partSkipCount = c / 64;
        int shiftSize = c - partSkipCount * 64;
        if (shiftSize == 0) {
            for (int i = std::min(SA - 1 + partSkipCount, SR - 1); i >= partSkipCount; --i) {
                r[i] = a[i - partSkipCount];
            }
        } else {
            int start = std::min(SA - 1 + partSkipCount, SR - 1);
            if (SR-1 > start) {
                r[start + 1] = a[start - partSkipCount] >> (64 - shiftSize);
            }
            for (int i = start; i > partSkipCount; --i) {
                r[i] = (a[i - partSkipCount] << shiftSize) | (a[i - 1 - partSkipCount] >> (64 - shiftSize));
            }
            r[partSkipCount] = a[0] << shiftSize;
        }
        for (int i = 0; i < partSkipCount; ++i) {
            r[i] = 0;
        }
    }
    int32_t shl(Int r, ConstInt a, uint32_t c, int32_t sa) {
        int partSkipCount = c / 64;
        int shiftSize = c - partSkipCount * 64;
        int32_t resultMaxIndex = sa + partSkipCount;
        if (shiftSize == 0) {
            r[resultMaxIndex] = 0;
            for (int i = resultMaxIndex-1; i >= partSkipCount; --i) {
                r[i] = a[i - partSkipCount];
            }
        } else {
            r[resultMaxIndex] = a[resultMaxIndex - 1 - partSkipCount] >> (64 - shiftSize);
            for (int i = resultMaxIndex-1; i > partSkipCount; --i) {
                r[i] = (a[i - partSkipCount] << shiftSize) | (a[i - 1 - partSkipCount] >> (64 - shiftSize));
            }
            r[partSkipCount] = a[0] << shiftSize;
        }
        for (int i = 0; i < partSkipCount; ++i) {
            r[i] = 0;
        }
        return resultMaxIndex + 1;
    }

    // r = a << 1 (left binary shift by 'c' bits, aka multiplication by 2)
    template<int S> void shl_1(Int r) {
        debugAssert(!(r[S - 1] & (1ull << 63)), "shifting left by 1 won't fit for given size");
        for (int i = S - 1; i > 0; --i) {
            r[i] = (r[i] << 1) | (r[i - 1] >> 63);
        }
        r[0] <<= 1;
    }
    void shl_1(Int r, int32_t sr) {
        debugAssert(!(r[sr - 1] & (1ull << 63)), "shifting left by 1 won't fit for given size");
        for (int i = sr - 1; i > 0; --i) {
            r[i] = (r[i] << 1) | (r[i - 1] >> 63);
        }
        r[0] <<= 1;
    }

    // r = a^2
    template<int SR, int SA> void sqr(Int r, ConstInt a) {
        static_assert(SR == 1 || SR == 2*SA, "for now implementation only allows some size combinations");
        if constexpr (SR == 1) {
            r[0] = a[0] * a[0];
        } else if constexpr (SA == 1) {
            mul128(r[1], r[0], a[0], a[0]);
        } else {
            Limb tp[2 * SA - 1];
            mulLimb(tp, &a[1], a[0], SA - 1);
            for (int i = 2; i < SA; ++i) {
                tp[SA + i - 2] = mulAddLimb(&tp[2 * i - 2], &a[i], a[i - 1], SA - i);
            }
            for (int i = 0; i < SA; ++i) {
                mul128(r[2 * i + 1], r[2 * i], a[i], a[i]);
            }
            tp[2 * SA - 2] = 0;
            shl_1<2 * SA - 1>(tp);
            add<2 * SA - 1>(r + 1, r + 1, tp);
        }
    }
    template<int SA> void sqr(Int r, ConstInt a) {
        sqr<2 * SA, SA>(r, a);
    }
    void sqr(Int r, ConstInt a, int32_t sa, Int tmpStorage) {
        if (sa == 1) {
            mul128(r[1], r[0], a[0], a[0]);
        } else {
            mulLimb(tmpStorage, &a[1], a[0], sa);
            for (int i = 2; i < sa; ++i) {
                tmpStorage[sa + i - 2] = mulAddLimb(&tmpStorage[2 * i - 2], &a[i], a[i - 1], sa - i);
            }
            for (int i = 0; i < sa; ++i) {
                mul128(r[2 * i + 1], r[2 * i], a[i], a[i]);
            }
            auto storageSize = 2 * sa - 1;
            tmpStorage[storageSize - 1] = 0;
            shl_1(tmpStorage, storageSize);
            add(r + 1, r + 1, tmpStorage, storageSize, storageSize);
        }
    }

    template<int SR, int SA, int SB, int SRem> void div(Int r, ConstInt a, ConstInt b, Int rem) {
        static_assert(SR >= SA - SB, "the division result might not fit in a given buffer");
        Limb aa[SA + 2];
        Limb bb[SB];
        Limb tmp[SB + 1];
        aa[SA + 1] = 0;
        aa[SA] = 0;

        auto realSizeB = sizeInBits<SB>(b);
        auto shiftSize = (64 - realSizeB % 64) % 64;
        shl<SA + 1, SA>(aa, a, shiftSize);
        shl<SB, SB>(bb, b, shiftSize);

        auto realSizeA = sizeInBits<SA + 1>(aa);
        auto realLimbSizeA = (realSizeA / 64 + 1) - (realSizeA % 64 == 0);
        auto realLimbSizeB = (realSizeB / 64 + 1) - (realSizeB % 64 == 0);

        clear<SR>(r);
        for (int i = realLimbSizeA - realLimbSizeB; i >= 0; --i) {
            auto d = div128(aa[i + realLimbSizeB], aa[i + realLimbSizeB - 1], bb[realLimbSizeB - 1]);
            mulLimb(tmp, bb, d, SB);
            if (sub<SB + 1>(&aa[i], &aa[i], tmp)) {
                d -= 1;
                while (sub<SB>(&aa[i], bb, &aa[i])) {
                    d -= 1;
                }
            }
            if (d != 0)
                r[i] = d;
        }
        shr<SB + 1, SB + 1>(aa, aa, shiftSize);
        copy<SRem>(rem, aa);
    }

    // r = a / b
    template<int SR, int SA, int SB> void div(Int r, ConstInt a, ConstInt b) {
        if constexpr (SR == 1 && SA == 1 && SB == 1) {
            r[0] = a[0] / b[0];
        } else {
            Limb rem[SB];
            div<SR, SA, SB, SB>(r, a, b, rem);
        }
    }

    // r = a % b
    template<int SR, int SA, int SB> void mod(Int r, ConstInt a, ConstInt b) {
        if constexpr (SR == 1 && SA == 1 && SB == 1) {
            r[0] = a[0] % b[0];
        } else {
            Limb divResult[SA];
            div<SA, SA, SB, SR>(divResult, a, b, r);
        }
    }

    // r = (a^-1) % mod
    int64_t gcdExtended(int64_t a, int64_t b, int64_t& x, int64_t& y) {
        auto a2 = a % b;
        auto c0 = a / b;
        if (a2 < 1) {
            x = 1;
            y = 0;
            return b;
        }
        auto a3 = b % a2;
        auto c1 = b / a2;
        if (a3 < 1) {
            x = -c0;
            y = 1;
            return a2;
        }
        auto g = gcdExtended(a2, a3, x, y);
        y -= c1 * x;
        x -= c0 * y;
        return g;
    }
    template<int S> void modInv(Int r, ConstInt a, ConstInt m) {
        Limb aCopy[S];
        Limb mCopy[S];
        Limb g[S];
        Limb s[S + 1];

        copy<S>(aCopy, a);
        copy<S>(mCopy, m);
        clear<S + 1>(s);

        mp_size_t sSize;
        mpn_gcdext(g, s, &sSize, aCopy, S, mCopy, S);
        if (sSize < 0) {
            sub_firstOperandIsBigger<S, S, S>(r, m, s);
        } else {
            copy<S>(r, s);
        }
    }
    void modInv(Int r, ConstInt a, ConstInt mod, int32_t s, Int tmpStorage) {
        copy(tmpStorage, a, s, s);
        mpn_sec_invert(r, tmpStorage, mod, s, 2 * 64 * s, tmpStorage + s);
    }

    // r = gcd(a, b)
    template<int S> void gcd(Int r, ConstInt a, ConstInt b) {
        Limb aCopy[S];
        Limb bCopy[S];
        copy<S>(aCopy, a);
        copy<S>(bCopy, b);
        clear<S>(r);
        mpn_gcd(r, aCopy, realSize<S>(a), bCopy, realSize<S>(b));
    }

    // modulo operation using Barret Reduction algorithm
    template<int S> void modBarret(Int r, ConstInt a, ConstInt mod, ConstInt R, uint32_t k) {
        if (cmp<2*S,S>(a, mod) < 0) {
            copy<S>(r, a);
            return;
        }
        Limb u0[3 * S];
        Limb u1[2 * S];
        mul<3 * S, 2 * S, S>(u0, a, R);
        shr<S, 3*S>(u0, u0, k);
        mul_SeperateMemoryLocation_onlyLowLimbs_plus1<S>(u1, u0, mod);
        sub<S, S+1, S+1>(r, a, u1);
        if (cmp<S>(r, mod) >= 0) {
            sub<S>(r, r, mod);
        }
    }

    template<int S> void montgomeryReductionCanUseInput(Int u, Int x, ConstInt k, ConstInt m, uint32_t b) {
        mul_SeperateMemoryLocation_onlyLowLimbs<S>(u, x, k);
        muladd_SeperateMemoryLocation<S>(x, u, m);
        copy<S>(u, x+S);
        if (cmp<S>(u, m) >= 0) {
            sub_firstOperandIsBigger<S,S,S>(u, u, m);
        }
    }

    template<int S> void montgomeryReduction(Int u, ConstInt x, ConstInt k, ConstInt m, uint32_t b) {
        Limb t[2*S];
        copy<2 * S>(t, x);
        mul_SeperateMemoryLocation_onlyLowLimbs<S>(u, t, k);
        muladd_SeperateMemoryLocation<S>(t, u, m);
        copy<S>(u, t + S);
        if (cmp<S>(u, m) >= 0) {
            sub_firstOperandIsBigger<S,S,S>(u, u, m);
        }
    }
    template<int S> void createMontgomeryReductionMod(Int k, ConstInt m, uint32_t& b) {
        Limb rInv[S+1];
        Limb r[S + 1];
        r[S] = 1;
        clear<S>(r);
        b = S;
        if (S == 1 && !(m[0] & (1ull << 63))) {
            Limb a[S];
            int64_t s, x;
            mod<S, S+1, S>(a, r, m);
            gcdExtended(m[0], a[0], s, x);
            if (s < 0) {
                rInv[0] = m[0] + s;
            } else {
                rInv[0] = s;
            }
            rInv[1] = 0;
        } else {
            Limb mCopy[S + 1];
            copy<S>(mCopy, m);
            mCopy[S] = 0;
            modInv<S + 1>(rInv, r, mCopy);
        }
        Limb r_rInv[2 * S + 1];
        clear<S>(r_rInv);
        copy<S+1>(r_rInv+S, rInv);
        Limb one = 1;
        sub<2 * S + 1, 2 * S + 1, 1>(r_rInv, r_rInv, &one);
        div<S, 2 * S, S>(k, r_rInv, m);
    }
    template<int S> void convertToMontgomeryForm(Int r, ConstInt a, ConstInt n, uint32_t b) {
        Limb mulRes[2 * S];
        
        clear<S>(mulRes);
        copy<S>(mulRes + S, a);

        mod<S, 2 * S, S>(r, mulRes, n);
    }

    template<int S> void montgomeryMult(Int r, ConstInt A, ConstInt B, ConstInt k, ConstInt m, uint32_t b) {
        uint64_t t[2 * S];
        if constexpr (S == 1) {
            mul128(t[1], t[0], A[0], B[0]);
            r[0] = t[0] * k[0];
            auto res = mulAddLimb(t, r, m[0], 1);
            addCarry(0, t[1], t[1], res);
        } else {
            mulLimb(t, A, B[0], S);
            for (int i = 1; i < S; ++i) {
                t[S + i] = mulAddLimb(&t[i], A, B[i], S);
            }
            mulLimbNoLastLimb(r, t, k[0], S);
            for (int i = 1; i < S; ++i) {
                mulAddLimb(&r[i], t, k[i], S - i);
            }
            Limb c = 0;
            for (int i = 0; i < S; ++i) {
                auto res = mulAddLimb(&t[i], r, m[i], S);
                c = addCarry((uint8_t)c, t[S + i], t[S + i], res);
            }
        }
        copy<S>(r, t + S);
        if (cmp<S>(r, m) >= 0) {
            sub_firstOperandIsBigger<S, S, S>(r, r, m);
        }
    }
    template<int S> void montgomerySqr(Int r, ConstInt A, ConstInt k, ConstInt m, uint32_t b) {
        uint64_t t[2 * S];
        sqr<S>(t, A);
        if constexpr (S == 1) {
            r[0] = t[0] * k[0];
            auto res = mulAddLimb(t, r, m[0], 1);
            addCarry(0, t[1], t[1], res);
        } else {
            mulLimbNoLastLimb(r, t, k[0], S);
            for (int i = 1; i < S; ++i) {
                mulAddLimb(&r[i], t, k[i], S - i);
            }
            Limb c = 0;
            for (int i = 0; i < S; ++i) {
                auto res = mulAddLimb(&t[i], r, m[i], S);
                c = addCarry((uint8_t)c, t[S + i], t[S + i], res);
            }
        }
        copy<S>(r, t + S);
        if (cmp<S>(r, m) >= 0) {
            sub_firstOperandIsBigger<S, S, S>(r, r, m);
        }
    }

    uint64_t div(uint64_t a, uint64_t m, uint32_t lm1) {
        uint64_t t1, t0;
        mul128(t1, t0, m, a);
        uint64_t q = ((t1 + ((a - t1) >> 1)) >> lm1);
        return q;
    }

    uint64_t mod(uint64_t a, uint64_t m, uint32_t d, uint32_t lm1) {
        return a - div(a, m, lm1) * d;
    }

    template<int S, int PowerModSize> uint32_t mod(ConstInt a, uint32_t m, uint64_t inverse, const uint32_t* constPowerModTable) {
        auto lm1 = (uint32_t)(::sizeInBits(m) - 1);
        
        if constexpr (S == 1) {
            return (uint32_t)mod(a[0], inverse, m, lm1);
        }

        uint64_t sum = mod(a[0], inverse, m, lm1);
        int i = 1;
        for (; i < min<S, PowerModSize>(); ++i) {
            sum += mod(a[i], inverse, m, lm1) * constPowerModTable[i - 1];
        }
        if constexpr (S > PowerModSize) {
            uint64_t powerModTable[S];
            if constexpr (PowerModSize >= 2) {
                auto powerMod1 = constPowerModTable[0];
                powerModTable[i - 1] = constPowerModTable[PowerModSize - 2];
                for (int j = i; j < S; ++j) {
                    powerModTable[j] = mod(powerModTable[j - 1] * powerMod1, inverse, m, lm1);
                }
            } else {
                powerModTable[1] = mod(mod(MaxU64, inverse, m, lm1) + 1, inverse, m, lm1);
                for (int i = 2; i < S; ++i) {
                    powerModTable[i] = mod(powerModTable[i - 1] * powerModTable[1], inverse, m, lm1);
                }
            }
            for (; i < S; ++i) {
                sum += mod(a[i], inverse, m, lm1) * powerModTable[i];
            }
        }

        return (uint32_t)mod(sum, inverse, m, lm1);
    }

    template<int S> uint32_t mod(ConstInt a, uint32_t m, int tableEntryId) {
        auto inverse = MultiplicativeInverses[tableEntryId];
        auto constPowerModTable = &PowerModTable(tableEntryId, 0);
        return mod<S, 8>(a, m, inverse, constPowerModTable);
    }

    template<int S> uint32_t mod(ConstInt a, uint32_t m, uint64_t inverse) {
        if constexpr (S == 1) {
            return a[0] % m;
        }
        return mod<S, 0>(a, m, inverse, nullptr);
    }

    template<int S> uint32_t mod(ConstInt a, uint32_t m) {
        if constexpr (S == 1) {
            return a[0] % m;
        }
        auto inverse = getInverse(m);
        return mod<S, 0>(a, m, inverse, nullptr);
    }
}
