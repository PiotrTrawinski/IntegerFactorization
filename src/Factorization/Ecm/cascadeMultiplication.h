#pragma once

#include "../../PrecomputedTables/primeTable.h"
#include "curves/common.h"
#include "curves/twistedEdwards.h"
#include "common.h"
#include "multiplicationMethods/doubleAndAddMul.h"
#include "multiplicationMethods/wnafMul.h"
#include "multiplicationMethods/pracMul.h"

template<template<typename, typename> typename CurveType, typename Type, typename ModType>
void cascadeMulDoMultiplication(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    switch (context.mulMethod) {
    case EcmMulMethod::DoubleAndAdd: doubleAndAddMul(context, curve, p, n); break;
    case EcmMulMethod::Naf:          nafMul(context, curve, p, n);          break;
    case EcmMulMethod::WNaf3:        wnafMul(3, context, curve, p, n);      break;
    case EcmMulMethod::WNaf4:        wnafMul(4, context, curve, p, n);      break;
    case EcmMulMethod::DNaf:         dnafMul(context, curve, p, n);         break;
    case EcmMulMethod::Prac:         prac(context, curve, p, n);            break;
    }
}
template<template<typename, typename> typename CurveType, typename Type, typename ModType, typename T>
void cascadeMulDoMultiplication(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, const T& n) {
    switch (context.mulMethod) {
    case EcmMulMethod::DoubleAndAdd: doubleAndAddMulX(context, curve, p, n); break;
    case EcmMulMethod::Naf:          debugAssert(false);      break;
    case EcmMulMethod::WNaf3:        debugAssert(false);      break;
    case EcmMulMethod::WNaf4:        debugAssert(false);      break;
    case EcmMulMethod::DNaf:         debugAssert(false);      break;
    case EcmMulMethod::Prac:         debugAssert(false);      break;
    }
}

template<template<typename, typename> typename CurveType, typename Type, typename ModType> 
int ecmStage1Mul(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& point) {
    int i = 0;
    if (context.mulMethod == EcmMulMethod::Prac) {
        debugAssert(context.mulCascadeMethod == EcmMulCascadeMethod::Seperate);

        // prac requires multiplicands > 2
        for (uint64_t r = 2; r <= context.B1; r *= 2) {
            dbl(curve, point);
            context.out_dblCount += 1;
        }
        i = 2;
    } 
    if (context.mulCascadeMethod == EcmMulCascadeMethod::Seperate) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            for (uint64_t p = Primes_1_000_000[i]; p <= context.B1; p *= Primes_1_000_000[i]) {
                cascadeMulDoMultiplication(context, curve, point, (uint64_t)Primes_1_000_000[i]);
            }
        }
    } 
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::Powers) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            uint64_t q;
            do {
                q = p;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
            cascadeMulDoMultiplication(context, curve, point, q);
        }
    }
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntilOverflow) {
        uint64_t x = 1;
        uint64_t y = 1;
        // TODO: something is wrong here, not the same result as Seperate and Powers method.
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            do {
                y *= Primes_1_000_000[i];
                if (x > y || y > (1ull << 63)) {
                    cascadeMulDoMultiplication(context, curve, point, x);
                    y = Primes_1_000_000[i];
                }
                x = y;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
        }
        cascadeMulDoMultiplication(context, curve, point, x);
    }
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntil256Overflow) {
        BigIntFixedSize<4> x = 1;
        BigIntFixedSize<4> y = 1;
        BigIntFixedSize<4> limit = MaxU64;
        shl(limit, limit, 190);
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            do {
                mul(y, y, Primes_1_000_000[i]);
                if (x > y || y > limit) {
                    cascadeMulDoMultiplication(context, curve, point, x);
                    y = Primes_1_000_000[i];
                }
                x = y;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
        }
    }
    return i;
}
