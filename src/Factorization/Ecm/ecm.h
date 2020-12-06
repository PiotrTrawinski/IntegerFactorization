#pragma once
#include "../../BigInt/include.h"
#include "curves/common.h"
#include "common.h"
#include "cascadeMultiplication.h"
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <cmath>
#include <vector>

template<template<typename,typename> typename CurveType, typename Type, typename ModType> void mul(CurveType<Type,ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    EcmContext context;
    doubleAndAddMul(context, curve, p, n);
}

template<typename CurveType> typename CurveType::ValueType ecm_(EcmContext& context, CurveType& curve) {
    using T = typename CurveType::ValueType;

    T factor = T{ 1 };
    T zero = T{ 0 };
    T one = T{ 1 };
    auto modValue = getModValue(curve.mod);

    if (context.initialCurveSeed == MaxU64) {
        context.initialCurveSeed = curve.defaultSeed();
    }
    CurvePoint<T> point = curve.initializeCurveAndPoint(context.initialCurveSeed);

    for (std::size_t j = 0; j < context.curveCount; ++j) {
        context.out_curveDoneCount += 1;

        // Stage 1
        int i = ecmStage1Mul(context, curve, point);
        if (point.z != zero) {
            gcd(factor, point.z, curve.mod);
        }
        if (factor != one) {
            return factor;
        }

        // Stage 2
        if (context.B2 > context.B1) {
            uint64_t firstPrime = Primes_1_000_000[i];
            uint64_t prevPrime = firstPrime;
            i += 1;
            std::vector<int> primeDifferencesBetweenB1B2;
            for (; Primes_1_000_000[i] <= context.B2; ++i) {
                primeDifferencesBetweenB1B2.emplace_back(Primes_1_000_000[i] - prevPrime);
                prevPrime = Primes_1_000_000[i];
            }
            int maxDiff = *std::max_element(primeDifferencesBetweenB1B2.begin(), primeDifferencesBetweenB1B2.end());
            std::vector<CurvePoint<T>> diffTable(maxDiff / 2);
            diffTable[0] = point;
            dbl(curve, diffTable[0]);
            diffTable[1] = diffTable[0];
            dbl(curve, diffTable[1]);
            for (int j = 2; j < maxDiff/2; ++j) {
                diffTable[j] = diffTable[j - 1];
                if (context.mulMethod == EcmMulMethod::Prac) {
                    diffAdd(curve, diffTable[j], diffTable[j], diffTable[0], diffTable[j-1]);
                } else {
                    add(curve, diffTable[j], diffTable[0]);
                }
            }

            cascadeMulDoMultiplication(context, curve, point, firstPrime);
            
            T runningMult = point.z; // TODO: I'm not sure if how I'm calculating this is ok (just multiplying all point.z and gcd at end)
            for (int diff : primeDifferencesBetweenB1B2) {
                if (context.mulMethod == EcmMulMethod::Prac) {
                    debugAssert(false);
                    //diffAdd(curve, point, diffTable[diff / 2 - 1], point, );
                } else {
                    add(curve, point, diffTable[diff / 2 - 1]);
                }
                modMul(runningMult, runningMult, point.z, curve.mod);
            }
            if (runningMult != zero)
                gcd(factor, runningMult, curve.mod);
            if (factor != one && factor != modValue) {
                return factor;
            }
        }

        curve.generateNewCurveAndPoint(point);
    }

    return factor;
}

template<template<typename,typename> typename CurveType, typename ModType> BigIntValueType<ModType> ecm(EcmContext& context, const ModType& mod) {
    CurveType<BigIntValueType<ModType>, ModType> curve;
    curve.mod = mod;
    return ecm_(context, curve);
}

template<template<typename, typename> typename CurveType> BigInt ecm(EcmContext& context, const BigInt& n) {
    return n.visit([&context](auto&& a) { return BigInt{ ecm<CurveType>(context, a) }; });
}
