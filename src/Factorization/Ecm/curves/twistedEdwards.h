#pragma once
#include "common.h"
#include "../../../BigInt/modularArithmetic.h"
#include "../../../Utility/debugAssert.h"
#include "shortWeierstrass.h"
#include <tuple>
#include <array>
#include <iostream>
#include <unordered_map>

template<typename T, typename U=T> struct TwistedEdwardsExtended {
    using ValueType = T;
    
    U mod;
    std::array<T, 4> tmp;
    uint64_t k;

    uint64_t defaultSeed() {
        return 2;
    }
    CurvePoint<T> initializeCurveAndPoint(uint64_t seed) {
        k = seed;
        return generateCurvePoint(*this, k);
    }
    void generateNewCurveAndPoint(CurvePoint<T>& out_point) {
        k += 1;
        out_point = generateCurvePoint(*this, k);
    }
    CurvePoint<T> zero() {
        return CurvePoint<T> { getConstant(0, mod), getConstant(1, mod), getConstant(0, mod), getConstant(0, mod) }; // TODO check what should it be
    }
};

// Input  = Extended
// Output =
// Extended:     8M + 0S + 10D
// Projective:   7M + 0S + 10D
// MontgomeryXY: 4M + 0S + 12D
template<typename Type, typename ModType> void _addsub(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, CoordinateSystem outputSystem = CoordinateSystem::Extended, bool isAdd = true) {
    debugAssert(outputSystem == CoordinateSystem::Extended || outputSystem == CoordinateSystem::Projective || outputSystem == CoordinateSystem::MontgomeryXY, "Unsupported output coordinate system");
    auto [Px, Py, Pz, Pt, Qx, Qy, Qz, Qt, u0, u1, u2] = modArithm::createContext(
        curve.mod, p.x, p.y, p.z, p.t, q.x, q.y, q.z, q.t, curve.tmp[0], curve.tmp[1], curve.tmp[2]
    );

    if (isAdd) {
        u0 = Qy + Qx;
        u1 = Qy - Qx;
    } else {
        u1 = Qy + Qx;
        u0 = Qy - Qx;
    }
    u2 = Py - Px;
    u0 = u0 * u2;
    u2 = Py + Px;
    u1 = u1 * u2;
    Px = u1 - u0;
    Py = u1 + u0;

    u1 = Pz + Pz;
    u1 = u1 * Qt;
    u2 = Qz + Qz;
    u2 = u2 * Pt;
    if (isAdd) {
        u0 = u2 - u1;
        u1 = u1 + u2;
    } else {
        u0 = u2 + u1;
        u2 = u2 - u1;
        u1 = u2;
    }
    if (outputSystem == CoordinateSystem::Extended || outputSystem == CoordinateSystem::Projective) {
        Pz = Px * Py;
        Px = Px * u1;
        Py = Py * u0;
        if (outputSystem == CoordinateSystem::Extended) {
            Pt = u0 * u1;
        }
    } else {
        // safe version
        //Pz = Px - u0;
        //Pz = Pz * Py;
        //Px = Px + u0;
        //Px = Px * Py;

        // short version (might produce additional point at infinity (2 :: 0) additionally to (0 :: 0)
        Pz = Px - u0;
        Px = Px + u0;
    }
}
template<typename Type, typename ModType> void add(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
    _addsub(curve, p, q, outputSystem, true);
}
template<typename Type, typename ModType> void sub(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
    _addsub(curve, p, q, outputSystem, false);
}
template<typename Type, typename ModType> void diffAdd(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p, const CurvePoint<Type>& q, const CurvePoint<Type>& pmq) {
    r = p;
    add(curve, r, q);
}


// Input  = Extended || Projective
// Output =
// Extended:   4M + 4S + 6D
// Projective: 3M + 4S + 6D
template<typename Type, typename ModType> void dbl(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& p, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
    debugAssert(outputSystem == CoordinateSystem::Extended || outputSystem == CoordinateSystem::Projective, "Unsupported output coordinate system");
    auto [X, Y, Z, T, u, v, xxx] = modArithm::createContext(curve.mod, p.x, p.y, p.z, p.t, curve.tmp[0], curve.tmp[1], curve.tmp[2]);
    
    u = sqr(X);
    v = sqr(Y);
    X = sqr(X + Y);
    Y = u + v;
    u = u - v;
    xxx = Y - X;
    v = sqr(Z);
    v = dbl(v) + u; 
    if (outputSystem == CoordinateSystem::Extended)
        T = xxx * Y;
    X = xxx * v;
    Y = Y * u;
    Z = u * v;
}
template<typename Type, typename ModType> void dbl(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p) {
    r = p;
    dbl(curve, r);
}


// Input  = Extended || Projective
// Output =
// Extended:   11M + 3S + 10D
// Projective:  9M + 3S + 10D
template<typename Type, typename ModType> void tpl(TwistedEdwardsExtended<Type, ModType>& curve, CurvePoint<Type>& p, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
    debugAssert(outputSystem == CoordinateSystem::Extended || outputSystem == CoordinateSystem::Projective, "Unsupported output coordinate system");
    auto [X, Y, Z, T, a, b, c, d] = modArithm::createContext(curve.mod, p.x, p.y, p.z, p.t, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3]);

    a = sqr(X);
    a = -a;
    b = sqr(Y);
    c = b + a;
    d = b - a;
    T = sqr(Z);
    T = T + T;
    T = T - c;
    T = T + T;
    c = c * d;
    a = a * T;
    b = b * T;
    T = c - b;
    b = b + c;
    d = c + a;
    c = a - c;
    b = X * b;
    c = Y * c;
    a = Z * T;
    if (outputSystem == CoordinateSystem::Projective) {
        X = b * T;
        Y = c * d;
        Z = a * d;
    } else {
        d = Z * d;
        X = b * a;
        Y = c * d;
        Z = a * d;
        T = b * c;
    }
}

template<typename Type, typename ModType> CurvePoint<Type> generateCurvePoint(TwistedEdwardsExtended<Type, ModType>& curve, uint64_t k, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
    debugAssert(outputSystem == CoordinateSystem::Extended || outputSystem == CoordinateSystem::Projective || outputSystem == CoordinateSystem::MontgomeryXY, "Unsupported output coordinate system");
    debugAssert(k != 0);
    
    static std::array<Type, 12> t_; // TODO: what if multithreaded?
    CurvePoint<Type> result;
    CurvePoint<Type> T;
    auto[Px,Py,Pz,Pt,Tx,Ty,Tz,Tt,U,V,W,u0,u1,u2,u3,u4,u5,u6,u7,u8] = modArithm::createContext(curve.mod, 
        result.x, result.y, result.z, result.t, T.x, T.y, T.z, T.t, t_[0], t_[1], t_[2], t_[3], t_[4], t_[5], t_[6], t_[7], t_[8], t_[9], t_[10], t_[11]
    );

    std::array<Type, 8>* constants;
    static std::unordered_map<ModType, std::array<Type, 8>> precomputedConstantsPerMod; // TODO: what if multithreaded?
    auto mapConstants = precomputedConstantsPerMod.find(curve.mod);
    if (mapConstants != precomputedConstantsPerMod.end()) {
        constants = &mapConstants->second;
    } else {
        auto inserted = precomputedConstantsPerMod.emplace(
            curve.mod,
            std::array<Type, 8> {
                getConstant(9747, curve.mod),
                getConstant(15, curve.mod),
                getConstant(378, curve.mod),
                getConstant(1, curve.mod),
                getConstant(144, curve.mod),
                getConstant(2985984, curve.mod),
                getConstant(96, curve.mod),
                getConstant(5, curve.mod)
            }
        );
        constants = &(*inserted.first).second;
    }

    u0 = (*constants)[0]; //getConstant(9747, curve.mod);
    u0 = -u0;
    Tx = (*constants)[1]; //getConstant(15, curve.mod);
    Ty = (*constants)[2]; //getConstant(378, curve.mod);
    Tz = (*constants)[3]; //getConstant(1, curve.mod);
    mulShortWeierstrassProjective(curve.mod, u0.value(), T, k);
    u0 = (*constants)[4]; //getConstant(144, curve.mod);
    U = Tx + Tz;
    U = U + Tz;
    U = U + Tz;
    U = U * u0;
    V = Ty;
    W = (*constants)[5]; //getConstant(2985984, curve.mod);
    W = W * Tz;

    u0 = (*constants)[6]; //getConstant(96, curve.mod);
    u0 = u0 * U;
    u1 = W - u0;
    u2 = sqr(u1);
    u3 = sqr(u0);
    u4 = (*constants)[7]; //getConstant(5, curve.mod);
    u4 = u4 * u3;
    u8 = u2 - u4;
    u5 = sqr(u8);
    u5 = u5 * u8;
    u7 = u1 + u1;
    u7 = dbl(u7);
    if (outputSystem == CoordinateSystem::MontgomeryXY) {

    } else {
        Tx = u1 - u0;
        u6 = (*constants)[7]; //getConstant(5, curve.mod);
        u6 = u6 * u0;
        u6 = u6 + u1;
        Tx = Tx * u6;
        u6 = (*constants)[7]; //getConstant(5, curve.mod);
        u6 = u6 * u3;
        Ty = u2 - u6;
        u6 = u6 + u2;
        Tx = Tx * u6;
        u6 = sqr(Ty);
        Ty = Ty * u6;
        u6 = u7 * u0;
        Tt = sqr(u6);
        Tt = Tt * u6;
        u6 = sqr(U);
        u6 = u6 * Tx;
        u2 = u3 * u0;
        u1 = u1 * u2;
        u1 = u1 * V;
        u1 = u1 * W;
        u1 = dbl(u1);
        Pz = Ty + Tt;
        Tz = Pz;
        Pz = Pz * u6;
        Px = Tz * u1;
        Py = Ty - Tt;
        if (outputSystem == CoordinateSystem::Extended) {
            Tz = Py;
        }
        Py = Py * u6;
        if (outputSystem == CoordinateSystem::Extended) {
            Pt = Tz * u1;
        }
    }

    // if montgomery calculate b here...

    return result;
}
