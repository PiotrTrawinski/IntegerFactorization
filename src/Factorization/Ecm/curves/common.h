#pragma once
#include "../../../BigInt/modularArithmetic.h"

template<typename T> struct CurvePoint {
    T x;
    T y;
    T z;
    T t;
};

template<typename T> std::ostream& operator<<(std::ostream& out, const CurvePoint<T>& point) {
    out << "x = " << point.x << '\n';
    out << "y = " << point.y << '\n';
    out << "z = " << point.z << '\n';
    out << "t = " << point.t << '\n';
    return out;
}

enum class CoordinateSystem {
    Extended,
    Projective,
    MontgomeryXY
};

template<typename Type> CurvePoint<Type> toAffine(const CurvePoint<Type>& point, const Type& mod, CoordinateSystem coordinateSystem) {
    CurvePoint<Type> result;
    Type zInv;
    auto [X, Y, Z, T, Rx, Ry, Rz, Rt, zInvM] = modArithm::createContext(mod, point.x, point.y, point.z, point.t, result.x, result.y, result.z, result.t, zInv);
    if (coordinateSystem == CoordinateSystem::Extended) {
        
    } else if (coordinateSystem == CoordinateSystem::Projective) {
        zInvM = inv(Z);
        Rx = X * zInvM;
        Ry = Y * zInvM;
    }
    return result;
}
