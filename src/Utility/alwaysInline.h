#pragma once
#include "compilerMacros.h"

#if defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    #define ALWAYS_INLINE inline __attribute__((__always_inline__))
#else
    #define ALWAYS_INLINE __forceinline
#endif
