#pragma once
#include <vector>
#include <cstdint>

struct PollardRhoParams {
    PollardRhoParams(uint64_t batchIterSize, uint64_t maxIterCount) : batchIterSize(batchIterSize), maxIterCount(maxIterCount), out_iterCount(0) {}
    uint64_t batchIterSize;
    uint64_t maxIterCount;
    uint64_t out_iterCount;
};

template<typename ModType> BigIntValueType<ModType> pollardRhoBrent(PollardRhoParams& params, const ModType& mod) {
    using T = BigIntValueType<ModType>;
    const T one = T{ 1 };
    const T oneConst = getConstant(1, mod);
    auto n = getModValue(mod);
    T x = getConstant(2, mod);
    T d = oneConst;
    T y, xs, dtmp;

    params.out_iterCount = 0;
    bool isDynamic = params.batchIterSize == 0;
    if (isDynamic) params.batchIterSize = 1;

    for (;;) {
        uint64_t r = 1;
        for (;;) {
            y = x;
            for (uint64_t i = 0; i < r; ++i) {
                modSqr(x, x, mod); add(x, x, oneConst);
            }
            params.out_iterCount += r;
            uint64_t k = 0;
            while (k < r) {
                xs = x;
                uint64_t end = std::min<uint64_t>(params.batchIterSize, r - k);
                for (uint64_t i = 0; i < end; ++i) {
                    modSqr(x, x, mod); add(x, x, oneConst);
                    absSub(dtmp, x, y);
                    modMul(d, d, dtmp, mod);
                }
                params.out_iterCount += end;
                if (isDynamic) {
                    params.batchIterSize = 10;
                    if (params.out_iterCount >= 1000) {
                        params.batchIterSize = 500;
                    } else if (params.out_iterCount >= 100) {
                        params.batchIterSize = 100;
                    }
                }
                if (!isZero(d)) {
                    gcd(d, d, n);
                    if (d != one)
                        return d;
                } else {
                    // try backtracking
                    for (uint64_t i = 0; i < params.batchIterSize - 1; ++i) {
                        modSqr(xs, xs, mod); add(xs, xs, oneConst);
                        absSub(d, xs, y);
                        params.out_iterCount += 1;
                        if (isZero(d)) {
                            break;
                        }
                        gcd(d, d, n);
                        if (d != one)
                            return d;
                    }
                    if (params.out_iterCount >= params.maxIterCount) {
                        return one;
                    }
                    break; // fail. Start again with new 'x'
                }
                if (params.out_iterCount >= params.maxIterCount) {
                    return one;
                }
                k += params.batchIterSize;
            }
            if (isZero(d)) {
                break;
            }
            r *= 2;
        }
        add(x, x, one); // new starting x
    }
    return one;
}

template<typename ModType> BigIntValueType<ModType> pollardRho(PollardRhoParams& params, const ModType& mod) {
    using T = BigIntValueType<ModType>;
    const T one = T{ 1 };
    const T oneConst = getConstant(1, mod);
    auto n = getModValue(mod);
    T x = getConstant(2, mod);
    T d = oneConst;
    T y, ys, xs, dtmp;
    params.out_iterCount = 0;
    bool isDynamic = params.batchIterSize == 0;
    if (isDynamic) params.batchIterSize = 1;
    for (;;) {
        y = x;
        for (;;) {
            xs = x;
            ys = y;
            for (uint64_t i = 0; i < params.batchIterSize; ++i) {
                modSqr(x, x, mod); add(x, x, oneConst);
                modSqr(y, y, mod); add(y, y, oneConst);
                modSqr(y, y, mod); add(y, y, oneConst);
                absSub(dtmp, y, x);
                modMul(d, d, dtmp, mod);
            }
            params.out_iterCount += params.batchIterSize;
            if (isDynamic) {
                params.batchIterSize = 10;
                if (params.out_iterCount >= 1000) {
                    params.batchIterSize = 500;
                } else if (params.out_iterCount >= 100) {
                    params.batchIterSize = 100;
                }
            }
            if (!isZero(d)) {
                gcd(d, d, n);
            } else {
                // try backtracking
                for (uint64_t i = 0; i < params.batchIterSize - 1; ++i) {
                    modSqr(xs, xs, mod); add(xs, xs, oneConst);
                    modSqr(ys, ys, mod); add(ys, ys, oneConst);
                    modSqr(ys, ys, mod); add(ys, ys, oneConst);
                    absSub(d, ys, xs);
                    params.out_iterCount += 1;
                    if (isZero(d)) {
                        break;
                    }
                    gcd(d, d, n);
                    if (d != one)
                        return d;
                }
                if (params.out_iterCount >= params.maxIterCount) {
                    return one;
                }
                break; // fail. Start again with new 'x'
            }
            if (d != one)
                return d;

            if (params.out_iterCount >= params.maxIterCount) {
                return one;
            }
        }
        add(x, x, one); // new starting x
    }
    return one;
}

BigInt pollardRhoBrent(PollardRhoParams& params, const BigInt& n) {
	return n.visit([&params](auto&& a) { return BigInt{ pollardRhoBrent(params, a) }; });
}
BigInt pollardRho(PollardRhoParams& params, const BigInt& n) {
	return n.visit([&params](auto&& a) { return BigInt{ pollardRho(params, a) }; });
}
