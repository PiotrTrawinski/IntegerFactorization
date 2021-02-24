#pragma once
#include "common.h"
#include "../../../BigInt/modularArithmetic.h"
#include "../../../Utility/debugAssert.h"
#include <tuple>
#include <array>

// 13M + 3D
template<typename Type, typename ModType> void _twistedHessianAddsub(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, bool isAdd = true) {
    auto [a, Px, Py, Pz, Qx, Qy, Qz, A, B, C, D, E, F, CD, FA, BE] = modArithm::createContext(
        curve.mod, curve.a, p.x, p.y, p.z, q.x, q.y, q.z, curve.tmp[0], curve.tmp[1], 
        curve.tmp[2], curve.tmp[3], curve.tmp[4], curve.tmp[5], curve.tmp[6], curve.tmp[7], curve.tmp[8]
    );
    A  = Px * Qz;
    B  = Pz * Qz;
    C  = Py * Qx;
    D  = Py * Qy;
    E  = Pz * Qy;
    F  = a * Px * Qx;
    CD = C * D;
    FA = F * A;
    BE = B * E;
    Px = A * B - CD;
    Py = D * E - FA;
    Pz = F * C - BE;
}
template<typename Type, typename ModType> void twistedHessianAdd(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _twistedHessianAddsub(curve, p, q, true);
}
template<typename Type, typename ModType> void twistedHessianSub(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _twistedHessianAddsub(curve, p, q, false);
}


// 8M + 1S + 7D
template<typename Type, typename ModType> void twistedHessianDbl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    auto [d, X, Y, Z, P, P2, S, A, C, D, E, Y2, dX2P] = modArithm::createContext(
        curve.mod, curve.d, p.x, p.y, p.z, curve.tmp[0], curve.tmp[1], curve.tmp[2],
        curve.tmp[3], curve.tmp[4], curve.tmp[5], curve.tmp[6], curve.tmp[7], curve.tmp[8]
    );
    P  = Y * Z;
    P2 = P * getConstant(2, curve.mod);
    S  = Y + Z;
    A  = sqr(S) - P;
    C  = (A - P2) * S;
    D  = (Z - Y) * A;
    dX2P = d * X * P2;
    E  = C * getConstant(3, curve.mod) - dX2P;
    X  = X * D * getConstant(-2, curve.mod);
    Y2 = Z * (D - E);
    Z  = Y * (D + E);
    Y  = Y2;
}


// 10M + 6S + 12D
template<typename Type, typename ModType> void twistedHessianTpl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    auto [a, dinv, X, Y, Z, U, V, W, A, B, C, D, E, T] = modArithm::createContext(
        curve.mod, curve.a, curve.dinv, p.x, p.y, p.z, curve.tmp[0], curve.tmp[1], curve.tmp[2],
        curve.tmp[3], curve.tmp[4], curve.tmp[5], curve.tmp[6], curve.tmp[7], curve.tmp[8]
    );
    U = sqr(X) * X * a;
    V = sqr(Y) * Y;
    W = sqr(Z) * Z;
    A = sqr(U - V);
    B = sqr(U - W);
    C = sqr(V - W);
    D = A + C;
    E = A + B;
    T = B + D;
    X = (U + V + W) * dinv * T;
    T = V * (C - E);
    Y = U * C * getConstant(2, curve.mod) - T;
    T = U * (B - D);
    Z = V * B * getConstant(2, curve.mod) - T;
}
