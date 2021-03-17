#pragma once
#include "../BigInt/include.h"

template<typename T> bool millerRabinTest(const T& n, T d) {
    T a;
    if constexpr (std::is_integral_v<T>) {
        a = random<int>(2, n - 2);
    } else {
        a = n.size() > 1 ? random(2, std::numeric_limits<int>::max()) : random<int>(2, n[0] - 2);
    }

    T x;
    modPow(x, a, d, n);

    auto one = T{ 1 };
    auto nm1 = n - T{ 1 };
    if (x == one || x == nm1) {
        return true;
    }

    while (d != nm1) {
        modSqr(x, x, n);
        shl(d, d, 1);
        if (x == one) return false;
        if (x == nm1) return true;
    }

    return false;
}

template<typename T> bool millerRabinTest(const T& n, int k=24) {
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

bool millerRabinTest(uint64_t n, int k=24) {
    return millerRabinTest(BigIntFixedSize<1>{n}, k);
}

bool millerRabinTest(const BigInt& n, int k=24) {
    return n.visitNoInit([k](auto&& a) { return millerRabinTest(a.mod, k); });
}
