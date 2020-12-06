#pragma once
#include "common.h"
#include "../../../BigInt/modularArithmetic.h"
#include "../../../Utility/debugAssert.h"
#include <tuple>
#include <array>
#include <iostream>

template<typename T, typename U=T> struct ShortWeierstrassProjective {
    using ValueType = T;
    
    U mod;
    T a;
    std::array<T, 5> tmp;

    uint64_t defaultSeed() {
        return 2;
    }
    CurvePoint<T> initializeCurveAndPoint(uint64_t seed) {
        a = getConstant(seed, mod);
        return CurvePoint<T> { getConstant(1, mod), getConstant(2, mod), getConstant(1, mod), 0 };
    }
    void generateNewCurveAndPoint(CurvePoint<T>& out_point) {
        add(a, a, getConstant(1, mod));
        out_point.z = getConstant(1, mod);
    }
    CurvePoint<T> zero() {
        return CurvePoint<T> { getConstant(0, mod), getConstant(1, mod), getConstant(0, mod), 0 };
    }
};

// 12M + 2S + 7D
template<typename Type, typename ModType> void _addsub(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, bool isAdd) {
    auto [Px, Py, Pz, Qx, Qy, Qz, u0, u1, u2, u3, u4] = modArithm::createContext(
        curve.mod, p.x, p.y, p.z, q.x, q.y, q.z, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4]
    );
    if (isAdd) {
        u4 = Qy;
    } else {
        u4 = -Qy;
    }
    u0 = Px * Qz;
    u1 = Py * Qz;
    u2 = Pz * Qz;
    u3 = u4 * Pz - u1;
    Px = Qx * Pz - u0;
    Pz = sqr(Px);
    Py = Pz * u0;
    u0 = Pz * Px;
    Pz = sqr(u3) * u2 - u0 - Py - Py;
    u1 = u1 * u0;
    Px = Px * Pz;
    Py = ((Py - Pz) * u3) - u1;
    Pz = u0 * u2;
}
template<typename Type, typename ModType> void add(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _addsub(curve, p, q, true);
}
template<typename Type, typename ModType> void sub(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _addsub(curve, p, q, false);
}
template<typename Type, typename ModType> void diffAdd(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p, const CurvePoint<Type>& q, const CurvePoint<Type>& pmq) {
    r = p;
    add(curve, r, q);
}

// 6M + 6S + 12D
template<typename Type, typename ModType> void dbl(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& p) {
    auto [a, Px, Py, Pz, u0, u1, u2, u3, u4] = modArithm::createContext(
        curve.mod, curve.a, p.x, p.y, p.z, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4]
    );
    u0 = sqr(Px);
    u1 = Py * Pz;
    u1 = dbl(u1);
    u2 = Py * u1;
    u3 = sqr(u2);
    u4 = sqr(Pz) * a + u0 + u0 + u0;
    Pz = sqr(Px + u2) - u0 - u3;
    u2 = sqr(u4) - Pz - Pz;
    Py = (Pz - u2) * u4 - u3 - u3;
    Px = u2 * u1;
    Pz = sqr(u1) * u1;
}
template<typename Type, typename ModType> void dbl(ShortWeierstrassProjective<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p) {
    r = p;
    dbl(curve, r);
}

template<typename Type, typename ModType> void mulShortWeierstrassProjective(const ModType& m, const Type& a, CurvePoint<Type>& p, uint64_t n) {
    static ShortWeierstrassProjective<Type, ModType> tmpCurve;
    tmpCurve.mod = m;
    tmpCurve.a = a;
    mul(tmpCurve, p, n);
}
