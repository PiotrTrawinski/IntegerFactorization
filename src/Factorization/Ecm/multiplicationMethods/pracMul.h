#pragma once

#include "../common.h"
#include "../curves/common.h"
#include "../curves/montgomery.h"
#include <cstdint>
#include <algorithm>

constexpr auto MontgomeryAddCost = 6.0 /* number of multiplications in an addition */;
constexpr auto MontgomeryDblCost = 5.0 /* number of multiplications in a duplicate */;

double lucasCost(uint64_t n, double v) {
    uint64_t d, e, r;
    double c; /* cost */
    d = n;
    r = (uint64_t)((double)d * v + 0.5);
    if (r >= n)
        return (MontgomeryAddCost * (double)n);
    d = n - r;
    e = 2 * r - n;
    c = MontgomeryDblCost + MontgomeryAddCost; /* initial duplicate and final addition */
    while (d != e) {
        if (d < e) {
            r = d;
            d = e;
            e = r;
        }
        if (d - e <= e / 4 && ((d + e) % 3) == 0) { /* condition 1 */
            d = (2 * d - e) / 3;
            e = (e - d) / 2;
            c += 3.0 * MontgomeryAddCost; /* 3 additions */
        } else if (d - e <= e / 4 && (d - e) % 6 == 0) { /* condition 2 */
            d = (d - e) / 2;
            c += MontgomeryAddCost + MontgomeryDblCost; /* one addition, one duplicate */
        } else if ((d + 3) / 4 <= e) { /* condition 3 */
            d -= e;
            c += MontgomeryAddCost; /* one addition */
        } else if ((d + e) % 2 == 0) { /* condition 4 */
            d = (d - e) / 2;
            c += MontgomeryAddCost + MontgomeryDblCost; /* one addition, one duplicate */
        }
        /* now d+e is odd */
        else if (d % 2 == 0) { /* condition 5 */
            d /= 2;
            c += MontgomeryAddCost + MontgomeryDblCost; /* one addition, one duplicate */
        }
        /* now d is odd and e is even */
        else if (d % 3 == 0) { /* condition 6 */
            d = d / 3 - e;
            c += 3.0 * MontgomeryAddCost + MontgomeryDblCost; /* three additions, one duplicate */
        } else if ((d + e) % 3 == 0) { /* condition 7 */
            d = (d - 2 * e) / 3;
            c += 3.0 * MontgomeryAddCost + MontgomeryDblCost; /* three additions, one duplicate */
        } else if ((d - e) % 3 == 0) { /* condition 8 */
            d = (d - e) / 3;
            c += 3.0 * MontgomeryAddCost + MontgomeryDblCost; /* three additions, one duplicate */
        } else /* necessarily e is even: catches all cases */
        { /* condition 9 */
            e /= 2;
            c += MontgomeryAddCost + MontgomeryDblCost; /* one addition, one duplicate */
        }
    }
    return c;
}

template<template<typename, typename> typename CurveType, typename Type, typename ModType> void prac(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t k) {
    constexpr auto NV = 10;
    auto& tmp = curve.tmp[0];

    /* 1/val[0] = the golden ratio (1+sqrt(5))/2, and 1/val[i] for i>0
       is the real number whose continued fraction expansion is all 1s
       except for a 2 in i+1-st place */
    static double val[NV] = { 
        0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
        0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
        0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
        0.61807966846989581 
    };

    /* for small n, it makes no sense to try 10 different Lucas chains */
    auto nv = std::min<uint64_t>(NV, sizeInLimbs(curve.mod));

    uint64_t i = 0;
    if (nv > 1) {
        /* chooses the best value of v */
        double cmin = MontgomeryAddCost * (double)k;
        for (uint64_t d = 0; d < nv; d++) {
            double c = lucasCost(k, val[d]);
            if (c < cmin) {
                cmin = c;
                i = d;
            }
        }
    }

    auto d = k;
    auto r = (uint64_t)((double)d * val[i] + 0.5);

    /* first iteration always begins by Condition 3, then a swap */
    d = k - r;
    uint64_t e = 2 * r - k;
    auto& A = p;
    auto B = p;
    auto C = p;
    auto T = p;
    auto U = p;

    dbl(curve, p); /* A = 2*A */
    while (d != e) {
        if (d < e) {
            r = d;
            d = e;
            e = r;
            swap(A.x, B.x);
            swap(A.z, B.z);
        }
        /* do the first line of Table 4 whose condition qualifies */
        if (d - e <= e / 4 && ((d + e) % 3) == 0) { /* condition 1 */
            d = (2 * d - e) / 3;
            e = (e - d) / 2;
            diffAdd(curve, T, A, B, C);
            diffAdd(curve, U, T, A, B);
            diffAdd(curve, B, B, T, A);
            swap(A.x, U.x);
            swap(A.z, U.z);
        } else if (d - e <= e / 4 && (d - e) % 6 == 0) { /* condition 2 */
            d = (d - e) / 2;
            diffAdd(curve, B, A, B, C);
            dbl(curve, A);
        } else if ((d + 3) / 4 <= e) { /* condition 3 */
            d -= e;
            diffAdd(curve, T, B, A, C);
            /* circular permutation (B,T,C) */
            tmp = B.x;
            B.x = T.x;
            T.x = C.x;
            C.x = tmp;
            tmp = B.z;
            B.z = T.z;
            T.z = C.z;
            C.z = tmp;
        } else if ((d + e) % 2 == 0) { /* condition 4 */
            d = (d - e) / 2;
            diffAdd(curve, B, B, A, C);
            dbl(curve, A);
        }
        /* now d+e is odd */
        else if (d % 2 == 0) { /* condition 5 */
            d /= 2;
            diffAdd(curve, C, C, A, B);
            dbl(curve, A);
        }
        /* now d is odd, e is even */
        else if (d % 3 == 0) { /* condition 6 */
            d = d / 3 - e;
            dbl(curve, T, A);
            diffAdd(curve, U, A, B, C);
            diffAdd(curve, A, T, A, A);
            diffAdd(curve, T, T, U, C);
            /* circular permutation (C,B,T) */
            tmp = C.x;
            C.x = B.x;
            B.x = T.x;
            T.x = tmp;
            tmp = C.z;
            C.z = B.z;
            B.z = T.z;
            T.z = tmp;
        } else if ((d + e) % 3 == 0) { /* condition 7 */
            d = (d - 2 * e) / 3;
            diffAdd(curve, T, A, B, C);
            diffAdd(curve, B, T, A, B);
            dbl(curve, T, A);
            diffAdd(curve, A, A, T, A);
        } else if ((d - e) % 3 == 0) { /* condition 8 */
            d = (d - e) / 3;
            diffAdd(curve, T, A, B, C);
            diffAdd(curve, C, C, A, B);
            swap(B.x, T.x);
            swap(B.z, T.z);
            dbl(curve, T, A);
            diffAdd(curve, A, A, T, A);
        } else /* necessarily e is even here */
        { /* condition 9 */
            e /= 2;
            diffAdd(curve, C, C, B, A);
            dbl(curve, B);
        }
    }
    diffAdd(curve, A, A, B, C);
    debugAssert(d == 1);
}
