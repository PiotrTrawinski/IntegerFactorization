#pragma once

enum class EcmMulMethod {
    DoubleAndAdd, // Double-And-Add
    Naf,          // 2NAF aka NAF
    WNaf3,        // 3NAF
    WNaf4,        // 4NAF
    DNaf,         // dynamic NAF
    Prac          // Montgomery PRAC algorithm
};

enum class EcmMulCascadeMethod {
    Seperate,            // All primes are multiplied seperately                  (example: B1=11 => nP = P*2*2*2*3*3*5*7*11)
    Powers,              // Multiplication done in batches of prime powers        (example: B1=11 => nP = P*8*9*5*7*11)
    MaxUntilOverflow,    // Multiplication done in batches as big as u64 can hold (example: B1=11 => nP = P*27720)
    MaxUntil256Overflow  // Multiplication done in batches as big as 256-bit can hold
};

struct EcmContext {
    EcmContext() {}
    EcmContext(EcmMulMethod mulMethod, EcmMulCascadeMethod mulCascadeMethod, uint64_t B1, uint64_t B2, uint64_t curveCount) :
        mulMethod(mulMethod), mulCascadeMethod(mulCascadeMethod), B1(B1), B2(B2), curveCount(curveCount)
    {}
    EcmContext(EcmMulMethod mulMethod, EcmMulCascadeMethod mulCascadeMethod, uint64_t B1, uint64_t B2, uint64_t curveCount, uint64_t initialCurveSeed) :
        EcmContext(mulMethod, mulCascadeMethod, B1, B2, curveCount)
    {
        this->initialCurveSeed = initialCurveSeed;
    }

    EcmMulMethod mulMethod;
    EcmMulCascadeMethod mulCascadeMethod;

    uint64_t B1;
    uint64_t B2;
    uint64_t curveCount;

    uint64_t initialCurveSeed = MaxU64;

    int out_dblCount = 0;
    int out_addCount = 0;
    int out_curveDoneCount = 0;
};
