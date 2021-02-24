#pragma once

#include "../../PrecomputedTables/primeTable.h"
#include "curves/common.h"
#include "curves/twistedEdwards.h"
#include "common.h"
#include "multiplicationMethods/doubleAndAddMul.h"
#include "multiplicationMethods/wnafMul.h"
#include "multiplicationMethods/pracMul.h"
#include "bytecode.h"

template<typename Type, typename ModType>
void cascadeMulDoMultiplication(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, uint64_t n) {
    switch (context.mulMethod) {
    case EcmMulMethod::DoubleAndAdd: doubleAndAddMul(context, curve, p, n); break;
    case EcmMulMethod::Naf:          nafMul(context, curve, p, n);          break;
    case EcmMulMethod::WNaf3:        wnafMul(3, context, curve, p, n);      break;
    case EcmMulMethod::WNaf4:        wnafMul(4, context, curve, p, n);      break;
    case EcmMulMethod::DNaf:         dnafMul(context, curve, p, n);         break;
    case EcmMulMethod::Prac:         prac(context, curve, p, n);            break;
    }
}
template<typename Type, typename ModType, typename T>
void cascadeMulDoMultiplication(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, const T& n) {
    switch (context.mulMethod) {
    case EcmMulMethod::DoubleAndAdd: doubleAndAddMulX(context, curve, p, n); break;
    case EcmMulMethod::Naf:          debugAssert(false);      break;
    case EcmMulMethod::WNaf3:        debugAssert(false);      break;
    case EcmMulMethod::WNaf4:        debugAssert(false);      break;
    case EcmMulMethod::DNaf:         debugAssert(false);      break;
    case EcmMulMethod::Prac:         debugAssert(false);      break;
    }
}

template<typename Type, typename ModType> 
int ecmStage1Mul(EcmContext& context, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& point) {
    int i = 0;
    if (context.mulMethod == EcmMulMethod::Prac) {
        // prac requires multiplicands > 2
        for (uint64_t r = 2; r <= context.B1; r *= 2) {
            dbl(curve, point);
            context.out_dblCount += 1;
        }
        i = 2;
    } 
    if (context.mulCascadeMethod == EcmMulCascadeMethod::Seperate) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            for (uint64_t p = Primes_1_000_000[i]; p <= context.B1; p *= Primes_1_000_000[i]) {
                cascadeMulDoMultiplication(context, curve, point, (uint64_t)Primes_1_000_000[i]);
            }
        }
    } 
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::Powers) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            uint64_t q;
            do {
                q = p;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
            cascadeMulDoMultiplication(context, curve, point, q);
        }
    }
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntilOverflow) {
        uint64_t x = 1;
        uint64_t y = 1;
        // TODO: something is wrong here, not the same result as Seperate and Powers method.
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            do {
                y *= Primes_1_000_000[i];
                if (x > y || y > (1ull << 63)) {
                    cascadeMulDoMultiplication(context, curve, point, x);
                    y = Primes_1_000_000[i];
                }
                x = y;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
        }
        cascadeMulDoMultiplication(context, curve, point, x);
    }
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntil256Overflow) {
        BigIntFixedSize<4> x = 1;
        BigIntFixedSize<4> y = 1;
        BigIntFixedSize<4> limit = MaxU64;
        shl(limit, limit, 190);
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            do {
                mul(y, y, Primes_1_000_000[i]);
                if (x > y || y > limit) {
                    cascadeMulDoMultiplication(context, curve, point, x);
                    y = Primes_1_000_000[i];
                }
                x = y;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
        }
    }
    return i;
}


template<typename ModType> void addBytecode(bytecode::Writer& bc, const EcmContext& context, EllipticCurveForm curveForm, const ModType& modValue, uint64_t n) {
    switch (context.mulMethod) {
    case EcmMulMethod::DoubleAndAdd: doubleAndAddMul(bc, n);                break;
    case EcmMulMethod::Naf:          nafMul(bc, n);                         break;
    case EcmMulMethod::WNaf3:        wNafMul(bc, n, 3);                     break;
    case EcmMulMethod::WNaf4:        wNafMul(bc, n, 4);                     break;
    case EcmMulMethod::DNaf:         dNafMul(bc, n, curveForm);             break;
    case EcmMulMethod::Prac:         pracMul(bc, n, sizeInLimbs(modValue)); break;
    }
}
template<typename ModType> std::pair<std::vector<uint8_t>, int> createBytecode(const EcmContext& context, EllipticCurveForm curveForm, const ModType& modValue) {
    bytecode::Writer bc;
    int i = 0;

    if (context.mulMethod == EcmMulMethod::Prac) {
        // prac requires multiplicands > 2
        bc.dbChainSTART(0, 0);
        for (uint64_t r = 2; r <= context.B1; r *= 2) {
            bc.dbChainDBL();
        }
        bc.dbChainEND();
        i = 2;
    } 
    if (context.mulCascadeMethod == EcmMulCascadeMethod::Seperate) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            for (uint64_t p = Primes_1_000_000[i]; p <= context.B1; p *= Primes_1_000_000[i]) {
                addBytecode(bc, context, curveForm, modValue, (uint64_t)Primes_1_000_000[i]);
            }
        }
    } 
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::Powers) {
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            uint64_t q;
            do {
                q = p;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
            addBytecode(bc, context, curveForm, modValue, q);
        }
    }
    else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntilOverflow) {
        uint64_t x = 1;
        uint64_t y = 1;
        // TODO: something is wrong here, not the same result as Seperate and Powers method.
        for (; Primes_1_000_000[i] <= context.B1; ++i) {
            uint64_t p = Primes_1_000_000[i];
            do {
                y *= Primes_1_000_000[i];
                if (x > y || y > (1ull << 63)) {
                    addBytecode(bc, context, curveForm, modValue, x);
                    y = Primes_1_000_000[i];
                }
                x = y;
                p *= Primes_1_000_000[i];
            } while (p <= context.B1);
        }
        addBytecode(bc, context, curveForm, modValue, x);
    }
    //else if (context.mulCascadeMethod == EcmMulCascadeMethod::MaxUntil256Overflow) {
    //    BigIntFixedSize<4> x = 1;
    //    BigIntFixedSize<4> y = 1;
    //    BigIntFixedSize<4> limit = MaxU64;
    //    shl(limit, limit, 190);
    //    for (; Primes_1_000_000[i] <= context.B1; ++i) {
    //        uint64_t p = Primes_1_000_000[i];
    //        do {
    //            mul(y, y, Primes_1_000_000[i]);
    //            if (x > y || y > limit) {
    //                cascadeMulDoMultiplication(context, curve, point, x);
    //                y = Primes_1_000_000[i];
    //            }
    //            x = y;
    //            p *= Primes_1_000_000[i];
    //        } while (p <= context.B1);
    //    }
    //}
    bc.END();
    return { bc.buffer.data, i };
}

template<typename Type, typename ModType> void runBytecode(const std::vector<uint8_t>& bytes, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& point) {
    bytecode::Reader bc(bytes);
    std::array<CurvePoint<Type>, 256> points;
    while (true) {
        auto block = bc.peekNextBlockOpCode();
        switch (block) {
        case bytecode::Block::Naf:
            runNafBlock(bc, curve, point, points);
            break;
        case bytecode::Block::DbChain:
            runDbChainBlock(bc, curve, point, points);
            break;
        case bytecode::Block::Prac:
            runPracBlock(bc, curve, point, points);
            break;
        case bytecode::Block::End:
            return;
        }
    }
}

template<typename Type, typename ModType> void runNafBlock(bytecode::Reader& bc, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, std::array<CurvePoint<Type>, 256>& points) {
    auto initPointCount = bc.peekDataBits();
    points[0] = p;
    if (initPointCount > 1) {
        points.back() = p;
        dbl(curve, points.back());
        for (int i = 1; i < initPointCount; ++i) {
            points[i] = points[i - 1];
            add(curve, points[i], points.back());
        }
        bc.skipByte();
        p = points[bc.peekByte()];
    }
    while (true) {
        bc.skipByte();
        switch (bc.peekNafOpCode()) {
        case bytecode::NafOpCode::ADD:
        case bytecode::NafOpCode::ADDn: {
            auto pInd = bc.peekDataBits();
            add(curve, p, points[pInd]);
            break;
        }
        case bytecode::NafOpCode::SUB:
        case bytecode::NafOpCode::SUBn: {
            auto pInd = bc.peekDataBits();
            sub(curve, p, points[pInd]);
            break;
        }
        case bytecode::NafOpCode::DBL:
        case bytecode::NafOpCode::DBLn:
            dbl(curve, p);
            break;
        case bytecode::NafOpCode::END:
            bc.skipByte();
            return;
        case bytecode::NafOpCode::FullMask: // should never happen, but added this case to avoid compiler warning
            return;
        }
    }
}

template<typename Type, typename ModType> void runDbChainBlock(bytecode::Reader& bc, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, std::array<CurvePoint<Type>, 256>& points) {
    auto initPointCount = bc.peekDataBits();
    points[0] = p;
    if (initPointCount > 1) {
        dbl(curve, p);
        for (int i = 1; i < initPointCount; ++i) {
            points[i] = points[i - 1];
            add(curve, points[i], p);
        }
        bc.skipByte();
        auto startPointIndex = bc.peekByte();
        if (startPointIndex != 0) {
            p = points[startPointIndex];
        }
    }
    bc.skipByte();
    while (true) {
        auto inst = bc.nextInstruction();
        for (int i = 0; i < inst.dblCount; ++i) {
            dbl(curve, p);
        }
        for (int i = 0; i < inst.tplCount; ++i) {
            tpl(curve, p);
        }
        if (!inst.skipAdd) {
            if (inst.isSub) {
                sub(curve, p, points[inst.index]);
            } else {
                add(curve, p, points[inst.index]);
            }
        }
        if (inst.isFinal) {
            break;
        }
    }
}

template<typename T> void pracSwap(CurvePoint<T>& A, CurvePoint<T>& B) {
    swap(A.x, B.x);
    swap(A.z, B.z);
}
template<typename T> void pracSwap3(T& tmp, CurvePoint<T>& A, CurvePoint<T>& B, CurvePoint<T>& C) {
    tmp = A.x;
    A.x = B.x;
    B.x = C.x;
    C.x = tmp;
    tmp = A.z;
    A.z = B.z;
    B.z = C.z;
    C.z = tmp;
}
template<typename Type, typename ModType> void runPracBlock(bytecode::Reader& bc, EllipticCurve<Type, ModType>& curve, CurvePoint<Type>& p, std::array<CurvePoint<Type>, 256>& points) {
    auto& tmp = curve.tmp[0];
    points[0] = p;
    points[1] = p;
    points[2] = p;
    points[3] = p;

    auto& A = p;
    auto& B = points[0];
    auto& C = points[1];
    auto& T = points[2];
    auto& U = points[3];
    
    dbl(curve, p);
    while (true) {
        bc.skipByte();
        int repCount = bc.peekRepCount();
        auto isSwap = bc.peekIfPracSwap();
        auto opCode = bc.peekPracOpCode();
        for (int i = 0; i < repCount; ++i) {
            if (isSwap) {
                pracSwap(A, B);
            }
            switch (opCode) {
            case bytecode::PracOpCode::Rule1:
                diffAdd(curve, T, A, B, C);
                diffAdd(curve, U, T, A, B);
                diffAdd(curve, B, B, T, A);
                pracSwap(A, U);
                break;
            case bytecode::PracOpCode::Rule2:
                diffAdd(curve, B, A, B, C);
                dbl(curve, A);
                break;
            case bytecode::PracOpCode::Rule3:
                diffAdd(curve, T, B, A, C);
                pracSwap3(tmp, B, T, C);
                break;
            case bytecode::PracOpCode::Rule4:
                diffAdd(curve, B, B, A, C);
                dbl(curve, A);
                break;
            case bytecode::PracOpCode::Rule5:
                diffAdd(curve, C, C, A, B);
                dbl(curve, A);
                break;
            case bytecode::PracOpCode::Rule6:
                dbl(curve, T, A);
                diffAdd(curve, U, A, B, C);
                diffAdd(curve, A, T, A, A);
                diffAdd(curve, T, T, U, C);
                pracSwap3(tmp, C, B, T);
                break;
            case bytecode::PracOpCode::Rule7:
                diffAdd(curve, T, A, B, C);
                diffAdd(curve, B, T, A, B);
                dbl(curve, T, A);
                diffAdd(curve, A, A, T, A);
                break;
            case bytecode::PracOpCode::Rule8:
                diffAdd(curve, T, A, B, C);
                diffAdd(curve, C, C, A, B);
                pracSwap(B, T);
                dbl(curve, T, A);
                diffAdd(curve, A, A, T, A);
                break;
            case bytecode::PracOpCode::Rule9:
                diffAdd(curve, C, C, B, A);
                dbl(curve, B);
                break;
            case bytecode::PracOpCode::End:
                goto pracEnd;
            }
        }
    }
pracEnd:
    bc.skipByte();
    diffAdd(curve, A, A, B, C);
}
