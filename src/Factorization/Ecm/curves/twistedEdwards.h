#pragma once
#include "common.h"
#include "../../../BigInt/modularArithmetic.h"
#include "../../../Utility/debugAssert.h"
#include "shortWeierstrass.h"
#include <tuple>
#include <array>
#include <unordered_map>


// Input  = Extended
// Output =
// Extended:     8M + 0S + 10D
// Projective:   7M + 0S + 10D
// MontgomeryXY: 4M + 0S + 12D
template<typename Type, typename ModType> void _twistedEdwardsAddsub(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q, CoordinateSystem outputSystem = CoordinateSystem::Extended, bool isAdd = true) {
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
template<typename Type, typename ModType> void twistedEdwardsAdd(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _twistedEdwardsAddsub(curve, p, q, CoordinateSystem::Extended, true);
}
template<typename Type, typename ModType> void twistedEdwardsSub(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    _twistedEdwardsAddsub(curve, p, q, CoordinateSystem::Extended, false);
}


// Input  = Extended || Projective
// Output =
// Extended:   4M + 4S + 6D
// Projective: 3M + 4S + 6D
template<typename Type, typename ModType> void twistedEdwardsDbl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    CoordinateSystem outputSystem = CoordinateSystem::Extended;
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


// Input  = Extended || Projective
// Output =
// Extended:   11M + 3S + 10D
// Projective:  9M + 3S + 10D
template<typename Type, typename ModType> void twistedEdwardsTpl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    CoordinateSystem outputSystem = CoordinateSystem::Extended;
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


template<typename Type, typename ModType> CurvePoint<Type> twistedEdwardsGenerateCurvePoint(EllipticCurve<Type, ModType>& curve, uint64_t k, CoordinateSystem outputSystem = CoordinateSystem::Extended) {
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


enum class TwistedEdwardsParametrization {
    Old,
    Z6,
    Z8,
    Z2x4,
    Z2x8,
    Z12,
};

template<typename T> auto div(int a, int b, T& u0, T& u1) {
    int absA = std::abs(a);
    int absB = std::abs(b);
    u0 = absB;
    if (b < 0)
        u0 = -u0;
    u0 = inv(u0);
    u1 = absA;
    if (a < 0)
        u1 = -u1;
    return u1 * u0;
}
template<typename T> auto set(int a, T& u0) {
    int absA = std::abs(a);
    u0 = absA;
    if (a < 0)
        u0 = -u0;
    return u0;
}
template<typename ModType> void twistedEdwardsInitializeParametrization(TwistedEdwardsParametrization parametrization, EllipticCurve<BigIntValueType<ModType>, ModType>& curve) {
    using T = BigIntValueType<ModType>;
    auto [a1, a2, a3, a4, s, t, sBase, tBase, u0, u1] = modArithm::createContext(curve.mod.mod, 
        curve.a1, curve.a2, curve.a3, curve.a4, curve.s, curve.t, curve.sBase, curve.tBase, curve.tmp[0], curve.tmp[1]
    );
    switch (parametrization) {
    case TwistedEdwardsParametrization::Z6:
        a1    = 0;
        a3    = 0;
        a2    = div(-1, 2304, u0, u1);
        a4    = div(-5, 221184, u0, u1);
        sBase = div(1, 192, u0, u1);
        tBase = div(1, 4608, u0, u1);
        s     = div(13, 2304, u0, u1);
        t     = div(-5, 18432, u0, u1);
        break;
    case TwistedEdwardsParametrization::Z8:
        a1    = 0;
        a3    = 0;
        a2    = 0;
        a4    = 48;
        sBase = 4;
        tBase = set(-16, u0);
        s     = 1;
        t     = 7;
        break;
    case TwistedEdwardsParametrization::Z2x4:
        a1    = 0;
        a3    = 0;
        a2    = 0;
        a4    = div(-11664, 25, u0, u1);
        sBase = 36;
        tBase = div(-864, 5, u0, u1);
        s     = div(2601, 100, u0, u1);
        t     = div(73899, 1000, u0, u1);
        break;
    case TwistedEdwardsParametrization::Z2x8:
        a1    = 0;
        a3    = 0;
        a2    = 0;
        a4    = set(-8, u0);
        sBase = 12;
        tBase = 40;
        s     = 12;
        t     = 40;
        break;
    case TwistedEdwardsParametrization::Z12:
        a1    = 0;
        a3    = 0;
        a2    = 0;
        a4    = set(-12, u0);
        sBase = set(-2, u0);
        tBase = set(-4, u0);
        s     = 4;
        t     = 4;
        break;
    }
}
template<typename Type, typename ModType> CurvePoint<Type> twistedEdwardsGenerateNextCurvePoint(TwistedEdwardsParametrization parametrization, EllipticCurve<Type, ModType>& curve) {
    CurvePoint<Type> result;
    auto [a1, a2, a3, a4, s, t, sBase, tBase, u, v, sig, r, alpha, beta, lambda, xnum, xden, ynum, yden, v_1, v_2, v_3, v_4, v_5, v_6] = modArithm::createContext(curve.mod.mod,
        curve.a1, curve.a2, curve.a3, curve.a4, curve.s, curve.t, curve.sBase, curve.tBase, 
        curve.tmp[0], curve.tmp[1], curve.tmp[2], curve.tmp[3], curve.tmp[4], curve.tmp[5], 
        curve.tmp[6], curve.tmp[7], curve.tmp[8], curve.tmp[9], curve.tmp[10], curve.tmp[11], 
        curve.tmp[12], curve.tmp[13], curve.tmp[14], curve.tmp[15], curve.tmp[16]
    );
#define DIV(a, b, u) ((u = inv(b)), (a) * u)
    v_1 = 1;
    v_2 = 2;
    v_3 = 3;
    v_4 = 4;
    v_5 = 5;
    v_6 = 6;
    switch (parametrization) {
    case TwistedEdwardsParametrization::Z6: {
        u = 96;
        sig = inv(u * s) - v_1;
        r = inv(sqr(s)) * t;
        alpha = sqr(sig) - v_5;
        beta = dbl(dbl(sig));
        xnum = dbl(r * sig);
        u = sig + v_5;
        v = sqr(sig) + v_5;
        xden = (sig - v_1) * u * v;
        u = sqr(alpha) * alpha;
        v = sqr(beta) * beta;
        ynum = u - v;
        yden = u + v;
        break;
    }
    case TwistedEdwardsParametrization::Z8:
        //u = divide(2 * s2, t2);
        //v = divide(2 * s2 * s2 * s2, t2 * t2) - 1;
        //x1num = 2 * u * u;
        //x1den = 1;
        //y1num = x1num * x1num - 1;
        //y1den = v;
        u = DIV(dbl(s), t, r);
        v = DIV(dbl(sqr(s) * s), sqr(t), r) - v_1;
        xnum = dbl(u * v);
        xden = 1;
        ynum = sqr(xnum) - v_1;
        yden = v;
        break;
    case TwistedEdwardsParametrization::Z2x4:
        //u = s2;
        //v = t2;
        //x1num = 1;
        //x1den = 3;
        //y1num = (30 * u - 25 * v + 1944) * (30 * u - 25 * v + 1944);
        //y1den = -625 * v * v + 77760 * v + 1250 * u * u * u + 24300 * u * u - 3779136;
        xnum = 1;
        xden = 3;
        alpha = 30;
        alpha = alpha * s;
        beta = 25;
        beta = beta * t;
        r = 1944;
        ynum = sqr(alpha - beta + r);

        u = 77760;
        u = u * t;

        v = 1250;
        v = v * s * s * s;

        alpha = 24300;
        alpha = alpha * s * s;

        beta = set(-3779136, r);

        yden = set(-625, r);
        yden = yden * t * t + u + v + alpha + beta;
        break;
    case TwistedEdwardsParametrization::Z2x8:
        //alpha = divide(1, divide(t2 + 25, s2 - 9) + 1);
        //beta = divide(2 * alpha * (4 * alpha + 1), 8 * alpha * alpha - 1);
        //x1num = (2 * beta - 1) * (4 * beta - 3);
        //x1den = (6 * beta - 5);
        //y1num = (2 * beta - 1) * (t2 * t2 + 50 * t2 - 2 * s2 * s2 * s2 + 27 * s2 * s2 - 104);
        //y1den = (t2 + 3 * s2 - 2) * (t2 + s2 + 16);
        u = 25;
        v = 9;
        beta = DIV(t + u, s - v, r) + v_1;
        alpha = inv(beta);
        u = dbl(dbl(alpha)) + v_1;
        beta = DIV(dbl(alpha) * u, dbl(dbl(dbl(sqr(alpha)))) - v_1, r);
        u = dbl(dbl(beta)) - v_3;
        xnum = (dbl(beta) - v_1) * u;
        xden = (v_6 * beta) - v_5;
        
        v = 50;
        v = v * t;
        sig = dbl(sqr(s) * s);
        u = 27;
        r = sqr(s) * u;
        alpha = 104;
        u = sqr(t) + v - sig + r - alpha;
        ynum = (dbl(beta) - v_1) * u;
        
        v = 16;
        u = t + s + v;
        yden = (v_3 * s + t - v_2) * u;
        break;
    case TwistedEdwardsParametrization::Z12:
        //u = divide(t2, 2 * s2);
        //v = divide(s2 * s2 - 4 * s2 - 12, s2 * s2 + 12 * s2 - 12);
        //x0 = divide(3 * v * v + 1, 4 * v);
        //y0 = divide(s2 * s2 + 12, 8 * s2 * v * (u * u + 3));
        //x1num = x0 * (v - 1) * 4 * u;
        //x1den = y0 * (v * v - 1) * (u * u + 3);
        //y1num = x0 - 1;
        //y1den = x0 + 1;
        u = DIV(t, dbl(s), r);
        r = 12;
        alpha = dbl(dbl(s));
        beta = r * s;
        v = DIV(sqr(s) - alpha - r, sqr(s) + beta - r, sig);
        alpha = DIV(sqr(v) * v_3 + v_1, dbl(dbl(v)), sig);
        lambda = sqr(u) + v_3;
        beta = DIV(sqr(s) + r, dbl(dbl(dbl(s))) * v * lambda, sig);
        xnum = (v - v_1) * alpha * v_4 * u;
        xden = (sqr(v) - v_1) * beta * lambda;
        ynum = alpha - v_1;
        yden = alpha + v_1;
        break;
    }

    u = xnum * yden; convertToMontgomeryForm(result.x, u.value(), curve.mod);
    u = xden * ynum; convertToMontgomeryForm(result.y, u.value(), curve.mod);
    u = xden * yden; convertToMontgomeryForm(result.z, u.value(), curve.mod);
    u = xnum * ynum; convertToMontgomeryForm(result.t, u.value(), curve.mod);

    if (s.value() == sBase.value()) {
        u = a1 * s;
        v = inv(dbl(t) + u + a3);
        alpha = sqr(s) * v_3;
        beta = a1 * t;
        lambda = (dbl(a2 * s) + alpha + a4 - beta) * v;
    } else {
        v = inv(s - sBase);
        lambda = (t - tBase) * v;
    }
    u = a1 * lambda - a2 - sBase - s;
	s = sqr(lambda) + u;
    u = a1 * s;
	t = (sBase - s) * lambda - tBase - u - a3;

    return result;
#undef DIV
}
