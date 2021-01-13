#pragma once
#include "common.h"
#include "shortWeierstrass.h"
#include "twistedEdwards.h"
#include "montgomery.h"

enum class EllipticCurveForm {
    ShortWeierstrass,
    TwistedEdwards,
    Montgomery
};

template<typename ValType, typename ModType = ValType> struct EllipticCurve {
    using ValueType = ValType;

    EllipticCurveForm form;
    ModType mod;
    std::array<ValType, 8> tmp;

    // ShortWeierstrass
    ValType a;

    // TwistedEdwards
    uint64_t k;

    // MontgomeryXY
    ValType a24;
    ValType sigma;

    EllipticCurve(EllipticCurveForm form) : form(form) {}

    uint64_t defaultSeed() {
        switch (form) {
        case EllipticCurveForm::ShortWeierstrass: return 2;
        case EllipticCurveForm::TwistedEdwards:   return 2;
        case EllipticCurveForm::Montgomery:     return 6;
        }
    }
    CurvePoint<ValType> initializeCurveAndPoint(uint64_t seed) {
        switch (form) {
        case EllipticCurveForm::ShortWeierstrass:
            a = getConstant(seed, mod);
            return CurvePoint<ValType> { getConstant(1, mod), getConstant(2, mod), getConstant(1, mod), 0 };
        case EllipticCurveForm::TwistedEdwards:
            k = seed;
            return twistedEdwardsGenerateCurvePoint(*this, k);
        case EllipticCurveForm::Montgomery:
            sigma = getConstant(seed, mod);
            return MontgomeryGenerateCurvePoint(*this, sigma);
        }
    }
    void generateNewCurveAndPoint(CurvePoint<ValType>& out_point) {
        switch (form) {
        case EllipticCurveForm::ShortWeierstrass:
            add(a, a, getConstant(1, mod));
            out_point.z = getConstant(1, mod);
            break;
        case EllipticCurveForm::TwistedEdwards:
            k += 1;
            out_point = twistedEdwardsGenerateCurvePoint(*this, k);
            break;
        case EllipticCurveForm::Montgomery:
            add(sigma, sigma, getConstant(1, mod));
            out_point = MontgomeryGenerateCurvePoint(*this, sigma);
            break;
        }
    }
    CurvePoint<ValType> zero() {
        switch (form) {
        case EllipticCurveForm::ShortWeierstrass: return CurvePoint<ValType> { getConstant(0, mod), getConstant(1, mod), getConstant(0, mod), 0 };
        case EllipticCurveForm::TwistedEdwards:   return CurvePoint<ValType> { getConstant(0, mod), getConstant(1, mod), getConstant(0, mod), getConstant(0, mod) };
        case EllipticCurveForm::Montgomery:       return CurvePoint<ValType> { getConstant(0, mod), 0, getConstant(0, mod), 0 }; // TODO: check if ok, pretty sure it is
        }
    }
};

template<typename Type, typename ModType> void add(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: shortWeierstrassAdd(curve, r, p); break;
    case EllipticCurveForm::TwistedEdwards:   twistedEdwardsAdd(curve, r, p); break;
    case EllipticCurveForm::Montgomery:       debugAssert(false, "cannot add using Montgomery form"); break;
    }
}
template<typename Type, typename ModType> void sub(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: shortWeierstrassSub(curve, r, p); break;
    case EllipticCurveForm::TwistedEdwards:   twistedEdwardsSub(curve, r, p); break;
    case EllipticCurveForm::Montgomery:       debugAssert(false, "cannot sub using Montgomery form"); break;
    }
}
template<typename Type, typename ModType> void diffAdd(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& r, const CurvePoint<Type>& p, const CurvePoint<Type>& q, const CurvePoint<Type>& pmq) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: debugAssert(false, "cannot diffAdd using ShortWeierstrass form"); break;
    case EllipticCurveForm::TwistedEdwards:   debugAssert(false, "cannot diffAdd using TwistedEdwards form"); break;
    case EllipticCurveForm::Montgomery:       montgomeryDiffAdd(curve, r, p, q, pmq); break;
    }
}
template<typename Type, typename ModType> void dbl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: shortWeierstrassDbl(curve, p); break;
    case EllipticCurveForm::TwistedEdwards:   twistedEdwardsDbl(curve, p); break;
    case EllipticCurveForm::Montgomery:       montgomeryDbl(curve, p, p); break;
    }
}
template<typename Type, typename ModType> void dbl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const CurvePoint<Type>& q) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: debugAssert(false, "not implemented"); break;
    case EllipticCurveForm::TwistedEdwards:   debugAssert(false, "not implemented"); break;
    case EllipticCurveForm::Montgomery:       montgomeryDbl(curve, p, q); break;
    }
}
template<typename Type, typename ModType> void tpl(EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p) {
    switch (curve.form) {
    case EllipticCurveForm::ShortWeierstrass: debugAssert(false, "cannot tpl using ShortWeierstrass form"); break;
    case EllipticCurveForm::TwistedEdwards:   twistedEdwardsTpl(curve, p); break;
    case EllipticCurveForm::Montgomery:       debugAssert(false, "cannot tpl using Montgomery form"); break;
    }
}
