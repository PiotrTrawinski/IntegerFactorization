#pragma once
#include "../BigInt/include.h"

template<typename T> bool millerRabinTest(const T& n, T d) {
    T a = n.size() > 1 ? random(2, std::numeric_limits<int>::max()) : random<int>(2, n[0] - 2);

    T x;
    modPow(x, a, d, n);

    auto nm1 = n - T{ 1 };
    if (x == 1 || x == nm1) {
        return true;
    }

    while (d != nm1) {
        modSqr(x, x, n);
        shl(d, d, 1);
        if (x == 1)   return false;
        if (x == nm1) return true;
    }

    return false;
}

template<typename T> bool millerRabinTest(T n, int k=24) {
    if (n == 2 || n == 3) {
        return true;
    }
    if (n == 0 || n == 1 || n[0] % 2 == 0) {
        return false;
    }

    // n = 2^s * d + 1  =>  d = (n-1)/(2^s)
    auto d = n - T{ 1 };
    while (d[0] == 0) {
        shr(d, d, 64);
    }
    shr(d, d, trailingZeroBitCount(d[0]));

    for (int i = 0; i < k; ++i) {
        if (!millerRabinTest(n, d)) {
            return false;
        }
    }
    return true;
}


