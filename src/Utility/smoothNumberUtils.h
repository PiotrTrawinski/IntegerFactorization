#pragma once

#include "../Factorization/TrialDivision.h"

template<typename T> bool isB1Smooth(const T& value, uint64_t B1) {
    auto a = value;
    while (true) {
        auto res = trialDivision(a, B1);
        if (res == a) {
            return res <= B1;
        }
        if constexpr (std::is_same_v<T, uint64_t>) {
            a = a / res;
        } else {
            a = a / res[0];
        }
    }
}

template<int Prime, typename T> uint64_t fastFactorPrime(T& a) {
    int count = 0;
    while (a % Prime == 0) {
        count += 1;
        a = a / Prime;
    }
    return Prime * count;
}

template<typename T> uint64_t smallestB1Value(const T& value, uint64_t bound) {
    auto a = value;
    uint64_t prevPrime = 1;
    uint64_t curPrime = 1;
    uint64_t maxPrime = 0;

    maxPrime = fastFactorPrime<3>(a);
    maxPrime = std::max(maxPrime, fastFactorPrime<5>(a));
    maxPrime = std::max(maxPrime, fastFactorPrime<7>(a));

    [[maybe_unused]] uint64_t biggestFactorChecked = 7;
    [[maybe_unused]] std::size_t currentTableEntryId = 0;
    [[maybe_unused]] std::size_t wheelIndex = 0;
    while (true) {
        T res;
        uint64_t u64Res;
        if constexpr (std::is_same_v<T, uint64_t>) {
            res = trialDivision(a, biggestFactorChecked, currentTableEntryId, wheelIndex, bound);
            u64Res = res;
        } else {
            res = trialDivision(a, bound);
            u64Res = res[0];
        }
        if (u64Res == prevPrime) {
            curPrime *= prevPrime;
        } else {
            prevPrime = u64Res;
            curPrime = prevPrime;
        }
        if (curPrime > maxPrime) {
            maxPrime = curPrime;
        }
        if (res == a) {
            return std::max<uint64_t>(u64Res, maxPrime);
        }
        a = a / u64Res;
    }
}
template<typename T> bool isB1PowerSmooth(const T& value, uint64_t B1) {
    return smallestB1Value(value, B1) <= B1;
}
