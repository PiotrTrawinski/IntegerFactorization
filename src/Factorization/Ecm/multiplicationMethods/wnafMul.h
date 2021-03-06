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
    StackVector() {}
    StackVector(const std::vector<T>& vec) {
        for (auto& v : vec) {
            emplace_back(v);
        }
    }
    void resize(int size, const T& value) {
        for (int i = 0; i < size; ++i) {
            emplace_back(value);
        }
    }

    int size() const { return size_; }
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }
    auto begin() { return data.begin(); }
    auto end() { return data.begin() + size(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.begin() + size(); }
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
            int8_t zi = e % (1 << w);
            if (w > 1 && zi >= (1 << (w - 1))) {
                zi -= (1 << w);
            }
            z.emplace_back(zi);
            e -= zi;
        } else {
            z.emplace_back(0);
        }
        e /= 2;
    }
    if (z.size() >= 3 && z[z.size() - 3] == -1) {
        z[z.size() - 3] = 1;
        z[z.size() - 2] = 1;
        z.size_ -= 1;
    }
    return z;
}
template<int64_t w, int S> StackVector<int8_t, 64*S> wnaf(BigIntFixedSize<S> e) {
    StackVector<int8_t, 64*S> z;
    while (e > 0) {
        if (e[0] & 1) {
            int8_t zi = e[0] & ((1 << w) - 1);
            if (w > 1 && zi >= (1 << (w - 1))) {
                zi -= (1 << w);
            }
            z.emplace_back(zi);
            e -= zi;
        } else {
            z.emplace_back(0);
        }
        e = e >> 1;
    }
    if (z.size() >= 3 && z[z.size() - 3] == -1) {
        z[z.size() - 3] = 1;
        z[z.size() - 2] = 1;
        z.size_ -= 1;
    }
    return z;
}
template<int64_t w> std::vector<int8_t> wnaf(BigIntGmp e) {
    std::vector<int8_t> z;
    while (e > 0) {
        if (e[0] & 1) {
            int8_t zi = e[0] & ((1 << w) - 1);
            if (w > 1 && zi >= (1 << (w - 1))) {
                zi -= (1 << w);
            }
            z.emplace_back(zi);
            sub(e, e, zi);
        } else {
            z.emplace_back(0);
        }
        shr(e, e, 1);
    }
    if (z.size() >= 3 && z[z.size() - 3] == -1) {
        z[z.size() - 3] = 1;
        z[z.size() - 2] = 1;
        z.erase(z.end() - 1);
    }
    return z;
}
template<int Size> int absoluteMaxNaf(const StackVector<int8_t, Size>& naf) {
    int max = 1;
    for (int i = 0; i < naf.size(); ++i) {
        max = std::max(max, std::abs(naf[i]));
    }
    return max;
}
int absoluteMaxNaf(const std::vector<int8_t>& naf) {
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
std::pair<int, int> nafDblAddCounts(const std::vector<int8_t>& naf) {
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
int nafCost(const std::vector<int8_t>& naf, int dblCost, int addCost, int intermediateDblCost, int intermediateAddCost) {
    auto [dblCount, addCount] = nafDblAddCounts(naf);
    return (intermediateDblCost * (dblCount - addCount) + dblCost * addCount) + (intermediateAddCost * (addCount - 1) + addCost);
}

template<typename Type, typename ModType> void nafMul(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
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
template<typename Type, typename ModType> void wnafMul(int w, EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;

    debugAssert(w == 3 || w == 4 || w == 5 || w == 6);
    auto nafForm =
        (w == 3) ? wnaf<3>(n) :
        (w == 4) ? wnaf<4>(n) :
        (w == 5) ? wnaf<5>(n) :
                   wnaf<6>(n);
    if (nafForm.size() == 1) {
        nafMul(context, curve, p, n);
        return;
    }
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
template<typename T, typename CostFunction> std::pair<int, std::vector<int8_t>> getBestWNaf(const T& n, CostFunction costFunction) {
    std::array<std::tuple<int, std::vector<int8_t>, double>, 5> nafForms = {
        std::tuple<int, std::vector<int8_t>, double>{ 2, wnaf<2>(n), 0 },
        std::tuple<int, std::vector<int8_t>, double>{ 3, wnaf<3>(n), 0 },
        std::tuple<int, std::vector<int8_t>, double>{ 4, wnaf<4>(n), 0 },
        std::tuple<int, std::vector<int8_t>, double>{ 5, wnaf<5>(n), 0 },
        std::tuple<int, std::vector<int8_t>, double>{ 6, wnaf<6>(n), 0 }
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

template<typename Type, typename ModType> void dnafMul(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;

    int bestWNaf = 0;
    if (curve.form == EllipticCurveForm::TwistedEdwards) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 8, 8, 8, 8); });
        bestWNaf = bestW;
    } else if (curve.form == EllipticCurveForm::Montgomery) {
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


#include "../bytecode.h"
template<typename T> void nafMul(bytecode::Writer& bc, const T& n) {
    bc.dbChainSTART();
    auto nafForm = wnaf<2>(n);
    for (int i = nafForm.size() - 2; i >= 0; --i) {
        bc.dbChainDBL();
        if (nafForm[i] == 1) {
            bc.dbChainADD(0);
        } else if (nafForm[i] == -1) {
            bc.dbChainSUB(0);
        }
    }
    bc.dbChainEND();
}

template<typename T> void wNafMul(bytecode::Writer& bc, const T& n, int w) {
    debugAssert(w == 3 || w == 4 || w == 5 || w == 6);
    auto nafForm =
        (w == 3) ? wnaf<3>(n) :
        (w == 4) ? wnaf<4>(n) :
        (w == 5) ? wnaf<5>(n) :
        wnaf<6>(n);
    if (nafForm.size() == 1) {
        nafMul(bc, n);
        return;
    }
    int tableSize = (absoluteMaxNaf(nafForm) + 1) / 2;

    bc.dbChainSTART(tableSize-1, (nafForm[nafForm.size() - 1] - 1) / 2);
    int start = nafForm.size() - 2;
    if (nafForm[nafForm.size() - 1] == 1 && tableSize != 1) {
        start -= 1;
    }
    for (int i = start; i >= 0; --i) {
        bc.dbChainDBL();
        if (nafForm[i] > 0) {
            bc.dbChainADD((nafForm[i] - 1) / 2);
        } else if (nafForm[i] < 0) {
            bc.dbChainSUB((-nafForm[i] - 1) / 2);
        }
    }
    bc.dbChainEND();
}

void dNafMul(bytecode::Writer& bc, uint64_t n, EllipticCurveForm curveForm) {
    int bestWNaf = 0;
    if (curveForm == EllipticCurveForm::TwistedEdwards) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 8, 8, 8, 8); });
        bestWNaf = bestW;
    } else if (curveForm == EllipticCurveForm::ShortWeierstrass) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 12, 14, 12, 14); });
        bestWNaf = bestW;
    } else {
        debugAssert(false, "Unknown curve type");
    }
    
    if (bestWNaf == 2) {
        nafMul(bc, n);
    } else {
        wNafMul(bc, n, bestWNaf);
    }
}

template<typename T> void dNafMul(bytecode::Writer& bc, const T& n, EllipticCurveForm curveForm) {
    int bestWNaf = 0;
    if (curveForm == EllipticCurveForm::TwistedEdwards) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 8, 8, 8, 8); });
        bestWNaf = bestW;
    } else if (curveForm == EllipticCurveForm::ShortWeierstrass) {
        auto [bestW, nafForm] = getBestWNaf(n, [](auto& a) { return nafCost(a, 12, 14, 12, 14); });
        bestWNaf = bestW;
    } else {
        debugAssert(false, "Unknown curve type");
    }
    
    if (bestWNaf == 2) {
        nafMul(bc, n);
    } else {
        wNafMul(bc, n, bestWNaf);
    }
}
