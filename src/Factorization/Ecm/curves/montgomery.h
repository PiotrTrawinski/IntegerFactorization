#pragma once
#include "common.h"
#include "../../../BigInt/modularArithmetic.h"
#include "../../../Utility/debugAssert.h"
#include <tuple>
#include <array>

// 2M + 2S + 1CM + 4D
template<typename Type, typename ModType> void montgomeryDbl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p) {
    auto [a24, Rx, Rz, Px, Pz, u0, u1, u2, u3, u4] = modArithm::createContext(
        curve.mod, curve.a24, r.x, r.z, p.x, p.z, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4]
    );
    u0 = Px + Pz;
    u1 = sqr(u0);
    u2 = Px - Pz;
    u3 = sqr(u2);
    u4 = u1 - u3;
    Rx = u1 * u3;
    Rz = (u4 * a24 + u3) * u4;
}

// 4M + 2S + 6D
template<typename Type, typename ModType> void montgomeryDiffAdd(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p, const CurvePoint<Type>& q, const CurvePoint<Type>& pmq) {
    auto [Rx, Rz, Px, Pz, Qx, Qz, PmQx, PmQz, u0, u1, u2, u3, u4, u5] = modArithm::createContext(
        curve.mod, r.x, r.z, p.x, p.z, q.x, q.z, pmq.x, pmq.z, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4], curve.tmp[5]
    );
    u0 = Px + Pz;
    u1 = Px - Pz;
    u2 = Qx + Qz;
    u3 = Qx - Qz;
    u4 = u0 * u3;
    u5 = u1 * u2;
    if (&r == &pmq) {
        u0 = sqr(u4 + u5) * PmQz;
        Rz = sqr(u4 - u5) * PmQx;
        Rx = u0;
    } else {
        Rx = sqr(u4 + u5) * PmQz;
        Rz = sqr(u4 - u5) * PmQx;
    }
}

template<typename Type, typename ModType> CurvePoint<Type> MontgomeryGenerateCurvePoint(EllipticCurve<Type, ModType>& curve, const Type& sigma) {
    // Brent-Suyama parametrization
    CurvePoint<Type> result;
    auto[Px,Pz,a24,s,a,u,v,vmu,v3u,u3,v4u3,c5] = modArithm::createContext(curve.mod, result.x, result.z, curve.a24, sigma, curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4], curve.tmp[5], curve.tmp[6], curve.tmp[7]);
    
    c5 = getConstant(5, curve.mod);;
    u = sqr(s) - c5;                        // u = sigma^2 - 5
    v = dbl(dbl(s));                        // v = 4*sigma
    vmu = v - u;                            // v - u
    v3u = dbl(u) + u + v;                   // 3u + v
    u3 = sqr(u) * u;                        // u^3
    v4u3 = inv(dbl(dbl(u3)) * v);           // 1 / (4u^3 * v)
    a = ((sqr(vmu) * vmu + v3u) * v4u3);    // a = ((v-u)^3 * (3u+v)) / (4u^3 * v)
    a24 = a;
    a24 >>= 2;
    Px = u3;                                // x = u^3
    Pz = sqr(v) * v;                        // z = v^3

    return result;
}
