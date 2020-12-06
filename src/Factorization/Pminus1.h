#pragma once
#include "../PrecomputedTables/primeTable.h"
#include "../BigInt/common.h"
#include <vector>
#include <cstdint>
#include <numeric>

template<typename Type, typename ModType> void SquareAndMultiply(Type& a, const ModType& m, uint64_t n) {
    auto b = a;
    for (auto i = mostSignificantBit(n) >> 1; i > 0; i >>= 1) {
		modSqr(a, a, m);
        if (n & i) {
            modMul(a, a, b, m);
        }
    }
}

template<typename ModType> BigIntValueType<ModType> pMinus1(const ModType& n, uint64_t B1, uint64_t B2) {
	using T = BigIntValueType<ModType>;
	auto one = T{ 1 };
	//auto modValue = getModValue(n);
    auto oneConst = getConstant(1, n);

	// Stage 1
	T x = getConstant(2, n);
    int i = 0;
    for (; Primes_1_000_000[i] <= B1; ++i) {
		uint64_t p = Primes_1_000_000[i];
		uint64_t q;
		do {
			q = p;
			p *= Primes_1_000_000[i];
		} while (p <= B1);
		SquareAndMultiply(x, n, q);
	}
    T xm1;
	sub(xm1, x, oneConst);
	if (isZero(xm1)) {
		return one;
	}
	T a;
	gcd(a, xm1, n);
	if (a != one || B1 >= B2) {
		return a;
	}

	// Stage 2
    uint64_t firstPrime = Primes_1_000_000[i];
    uint64_t prevPrime = firstPrime;
    i += 1;
    std::vector<int> primeDifferencesBetweenB1B2;
    for (; Primes_1_000_000[i] <= B2; ++i) {
        primeDifferencesBetweenB1B2.emplace_back(Primes_1_000_000[i] - prevPrime);
        prevPrime = Primes_1_000_000[i];
    }
    int maxDiff = *std::max_element(primeDifferencesBetweenB1B2.begin(), primeDifferencesBetweenB1B2.end());
    std::vector<T> diffTable(maxDiff / 2);
    diffTable[0] = x;
    modSqr(diffTable[0], diffTable[0], n);
    diffTable[1] = diffTable[0];
    modSqr(diffTable[1], diffTable[1], n);
    for (int j = 2; j < maxDiff / 2; ++j) {
        diffTable[j] = diffTable[j - 1];
        modMul(diffTable[j], diffTable[j], diffTable[0], n);
    }

    SquareAndMultiply(x, n, firstPrime);
    T runningMult = x;
    int gcdInterval = 100;
    int gcdCount = 0;
    for (int diff : primeDifferencesBetweenB1B2) {
        modMul(x, x, diffTable[diff / 2 - 1], n);
        modMul(runningMult, runningMult, x, n);
        if (gcdCount == gcdInterval) {
            sub(xm1, runningMult, oneConst);
            if (!isZero(xm1)) {
                gcd(a, xm1, n);
                if (a != one) {
                    return a;
                }
            }
            gcdCount = 0;
        } else {
            gcdCount += 1;
        }
    }
    sub(xm1, runningMult, oneConst);
    if (isZero(xm1)) {
        return one;
    }
    gcd(a, xm1, n);
	return a;
}

BigInt pMinus1(const BigInt& n, uint64_t B1, uint64_t B2) {
    return n.visit([B1,B2](auto&& a) { return BigInt{ pMinus1(a, B1, B2) }; });
}
