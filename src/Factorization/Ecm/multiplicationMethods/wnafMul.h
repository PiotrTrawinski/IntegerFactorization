#pragma once

#include "../common.h"
#include "../curves/common.h"
#include "../../../Utility/bitManipulation.h"
#include "../curves/twistedEdwards.h"
#include "../curves/shortWeierstrass.h"
#include <array>
#include <tuple>
#include <algorithm>

template<typename T, int N> struct StackVector {
    int size() const { return size_; }
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }
    auto begin() { return data.begin(); }
    auto end() { return data.begin() + size(); }
    void erase(int index) { 
        std::rotate(begin() + index, begin() + index + 1, end()); 
        size_ -= 1;
    }
    template<typename... Args> void emplace_back(Args&&... args) {
        new(&data[size_++]) T(args...);
    }

    std::array<T, N> data;
    int size_ = 0;
};

// includes optimization which changes from (4P - P) to (2P + P) => it's not "true" NAF
template<int64_t w=2> StackVector<int8_t, 64> wnaf(int64_t e) {
    StackVector<int8_t, 64> z;
    while (e > 0) {
        if (e % 2 == 1) {
            auto zi = e % (1 << w);
            if (w > 1 && zi >= (1 << (w - 1))) {
                zi -= (1 << w);
            }
            if (zi == 1 && z.size() >= 2 && z[z.size() - 2] == -1) {
                z[z.size() - 1] = 1;
                z.erase(0);
                z.emplace_back(1);
            } else {
                z.emplace_back((int8_t)zi);
            }
            e -= zi;
        } else {
            z.emplace_back(0);
        }
        e /= 2;
    }
    return z;
}
int absoluteMaxNaf(const StackVector<int8_t, 64>& naf) {
    int max = 1;
    for (int i = 0; i < naf.size(); ++i) {
        max = std::max(max, std::abs(naf[i]));
    }
    return max;
}
std::pair<int, int> nafDblAddCounts(const StackVector<int8_t, 64>& naf) {
    int dblCount = (naf[naf.size() - 1] != 1) + naf.size() - 1;
    int addCount = (absoluteMaxNaf(naf) + 1) / 2 - 1;
    for (int i = 0; i < naf.size() - 1; ++i) {
        if (naf[i] != 0) {
            addCount += 1;
        }
    }
    return { dblCount, addCount };
}
int nafCost(const StackVector<int8_t, 64>& naf, int dblCost, int addCost, int intermediateDblCost, int intermediateAddCost) {
    auto [dblCount, addCount] = nafDblAddCounts(naf);
    return (intermediateDblCost * (dblCount - addCount) + dblCost * addCount) + (intermediateAddCost * (addCount - 1) + addCost);
}

template<template<typename, typename> typename CurveType, typename Type, typename ModType> void nafMul(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;

    auto nafForm = wnaf<2>(n);
    CurvePoint<Type> q = p;
    for (int i = nafForm.size() - 2; i >= 0; --i) {
        dbl(curve, p);
        context.out_dblCount += 1;
        if (nafForm[i] == 1) {
            add(curve, p, q);
            context.out_addCount += 1;
        } else if (nafForm[i] == -1) {
            sub(curve, p, q);
            context.out_addCount += 1;
        }
    }
}
template<template<typename, typename> typename CurveType, typename Type, typename ModType> void wnafMul(int w, EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;

    debugAssert(w == 3 || w == 4 || w == 5 || w == 6);
    auto nafForm =
        (w == 3) ? wnaf<3>(n) :
        (w == 4) ? wnaf<4>(n) :
        (w == 5) ? wnaf<5>(n) :
                   wnaf<6>(n);
    int tableSize = (absoluteMaxNaf(nafForm) + 1) / 2;
    std::array<CurvePoint<Type>, 1 << 4> q; // NOTE: will fail to work for wNAF with w > 6!
    q[0] = p;
    dbl(curve, p);
    context.out_dblCount += 1;
    for (int i = 1; i < tableSize; ++i) {
        q[i] = q[i - 1];
        add(curve, q[i], p);
        context.out_addCount += 1;
    }
    int start = nafForm.size() - 3;
    if (nafForm[nafForm.size() - 1] != 1) {
        start += 1;
        p = q[(nafForm[nafForm.size() - 1] - 1) / 2];
    }
    for (int i = start; i >= 0; --i) {
        dbl(curve, p);
        context.out_dblCount += 1;
        if (nafForm[i] > 0) {
            add(curve, p, q[(nafForm[i] - 1) / 2]);
            context.out_addCount += 1;
        } else if (nafForm[i] < 0) {
            sub(curve, p, q[(-nafForm[i] - 1) / 2]);
            context.out_addCount += 1;
        }
    }
}

template<typename CostFunction> std::pair<int, StackVector<int8_t, 64>> getBestWNaf(uint64_t n, CostFunction costFunction) {
    std::array<std::tuple<int, StackVector<int8_t, 64>, double>, 5> nafForms = {
        std::tuple<int, StackVector<int8_t, 64>, double>{ 2, wnaf<2>(n), 0 },
        std::tuple<int, StackVector<int8_t, 64>, double>{ 3, wnaf<3>(n), 0 },
        std::tuple<int, StackVector<int8_t, 64>, double>{ 4, wnaf<4>(n), 0 },
        std::tuple<int, StackVector<int8_t, 64>, double>{ 5, wnaf<5>(n), 0 },
        std::tuple<int, StackVector<int8_t, 64>, double>{ 6, wnaf<6>(n), 0 }
    };
    for (auto& nafForm : nafForms) {
        std::get<2>(nafForm) = costFunction(std::get<1>(nafForm));
    }
    std::stable_sort(nafForms.begin(), nafForms.end(), [](auto& a, auto& b) {
        return std::get<2>(a) < std::get<2>(b);
    });
    auto& bestNafForm = nafForms[0];
    int bestWNaf = std::get<0>(bestNafForm);
    return { bestWNaf, std::get<1>(bestNafForm) };
}
template<template<typename, typename> typename CurveType, typename Type, typename ModType> void dnafMul(EcmContext& context, CurveType<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;

    int bestWNaf = 0;
    if constexpr (std::is_same_v<CurveType<Type, ModType>, TwistedEdwardsExtended<Type, ModType>>) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 8, 8, 8, 8); });
        bestWNaf = bestW;
    } else if constexpr (std::is_same_v<CurveType<Type, ModType>, ShortWeierstrassProjective<Type, ModType>>) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 12, 14, 12, 14); });
        bestWNaf = bestW;
    } else {
        debugAssert(false, "Unknown curve type");
    }
    
    if (bestWNaf == 2) {
        nafMul(context, curve, p, n);
    } else {
        wnafMul(bestWNaf, context, curve, p, n);
    }
}
