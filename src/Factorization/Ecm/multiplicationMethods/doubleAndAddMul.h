#pragma once

#include "../common.h"
#include "../curves/common.h"
#include "../../../Utility/bitManipulation.h"

template<typename Type, typename ModType> void doubleAndAddMul(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    if (n == 0) { p = curve.zero(); return; }
    if (n == 1) return;
    
    CurvePoint<Type> q = p;
    for (auto i = mostSignificantBit(n) >> 1; i > 0; i >>= 1) {
        dbl(curve, p);
        context.out_dblCount += 1;
        if (n & i) {
            add(curve, p, q);
            context.out_addCount += 1;
        }
    }
}
template<typename Type, typename ModType, typename T> void doubleAndAddMulX(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const T& n) {
    if (isZero(n)) { p = curve.zero(); return; }

    CurvePoint<Type> q = p;
    for (auto i = n.sizeInBits()-1; i > 0; --i) {
        dbl(curve, p);
        context.out_dblCount += 1;
        if (n.bit(i)) {
            add(curve, p, q);
            context.out_addCount += 1;
        }
    }
}

#include "../bytecode.h"
void doubleAndAddMul(bytecode::Writer& bc, uint64_t n) {
    bc.nafSTART();
    for (auto i = mostSignificantBit(n) >> 1; i > 0; i >>= 1) {
        bc.nafDBL();
        if (n & i) {
            bc.nafADD(0);
        }
    }
    bc.nafEND();
}
