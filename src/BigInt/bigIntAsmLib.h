#pragma once
#include <stdint.h>

#ifdef USE_ASM_LIB
extern "C" void montgomeryMult1(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult2(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult3(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult4(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult5(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult6(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult7(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomeryMult8(uint64_t* r, const uint64_t* A, const uint64_t* B, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr1(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr2(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr3(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr4(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr5(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr6(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr7(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
extern "C" void montgomerySqr8(uint64_t* r, const uint64_t* A, const uint64_t* k, const uint64_t* m, uint32_t b);
#endif
