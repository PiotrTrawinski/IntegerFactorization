#pragma once
#include "../PrecomputedTables/primeTable.h"
#include "../Utility/generalUtils.h"
#include "../Utility/avx2.h"
#include "../Utility/compilerMacros.h"
#include "../BigInt/include.h"
#include <vector>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <intrin.h>
#include <cmath>
#include <array>

constexpr uint64_t MaxIntegerDouble = 1ull << 53;

__m256d my256RemDoubleZeroMask(const AvxType<double>& a, const AvxType<double>& b) {
	return a == b * avx::floor(a / b);
}
__m256d my256RemDouble(const AvxType<double>& a, const AvxType<double>& b) {
	return a - b * avx::floor(a / b);
}
__m256d myConvertu32Tof64(const __m128i& a) {
	auto hasHighBitMask = _mm_cmplt_epi32(a, _mm_set1_epi32(0));
	auto hasHighBitMask64 = _mm256_cvtepi32_epi64(hasHighBitMask);
	auto allButHighBit = _mm_and_si128(a, _mm_set1_epi32(0x7fffffff));
	auto doubleAllButHighBit = _mm256_cvtepi32_pd(allButHighBit);
	auto doubleAddHightBitMask = _mm256_and_pd(_mm256_castsi256_pd(hasHighBitMask64), _mm256_set1_pd(0x80000000));
	return _mm256_add_pd(doubleAllButHighBit, doubleAddHightBitMask);
}
AvxType<uint64_t> Rem256_64(const AvxType<uint64_t>& a, const AvxType<uint64_t>& b) {
	/*
		a = H*2^12 + L
		a % b = (H*2^12 + L) % b = ((H % b) * (2^12 % b) + L%b) % b 
	*/
	auto aLowHigh = _mm256_permutevar8x32_epi32(a, _mm256_setr_epi32(0,2,4,6,1,3,5,7));
	auto aLow = myConvertu32Tof64(_mm256_castsi256_si128(_mm256_srli_epi32(_mm256_slli_epi32(aLowHigh, 20), 20)));
	auto aHigh = myConvertu32Tof64(_mm256_extracti128_si256(aLowHigh, 1));
	auto bLoww = _mm256_permutevar8x32_epi32(b, _mm256_setr_epi32(0,2,4,6,4,5,6,7));
	auto bLow = myConvertu32Tof64(_mm256_castsi256_si128(bLoww));

	auto aHighPartOfLow = myConvertu32Tof64(_mm256_castsi256_si128(_mm256_srli_epi32(aLowHigh, 12)));
	aHigh = _mm256_mul_pd(aHigh, _mm256_set1_pd(1024*1024));
	aHigh = _mm256_add_pd(aHigh, aHighPartOfLow);

	auto highModB = my256RemDouble(aHigh, bLow);
	auto cstModB = my256RemDouble(_mm256_set1_pd(4096), bLow);
	auto lowModB = my256RemDouble(aLow, bLow);
	auto res = _mm256_add_pd(_mm256_mul_pd(highModB, cstModB), lowModB);
	return _mm256_castpd_si256(my256RemDoubleZeroMask(res, bLow));
}

class AvxDividesResult {
	uint32_t value;
public:
	AvxDividesResult() {}
	AvxDividesResult(uint32_t value) : value(value) {}
	operator bool() { return value; }
	template<int Index> bool get() {
		if constexpr (Index == 0) return value & 0x000000ff;
		if constexpr (Index == 1) return value & 0x0000ff00;
		if constexpr (Index == 2) return value & 0x00ff0000;
		if constexpr (Index == 3) return value & 0xff000000;
		return false;
	}
};
AvxDividesResult divides(const AvxType<double>& a, const AvxType<double>& b, const AvxType<double>& bInverse) {
	return AvxDividesResult{ avx::mask8bit(a == b * avx::floor(a * bInverse)) };
}
AvxDividesResult divides(const AvxType<double>& a, const AvxType<double>& b) {
	return AvxDividesResult{ avx::mask8bit(a == b * avx::floor(a / b)) };
}
AvxDividesResult divides(const AvxType<uint64_t>& a, const AvxType<uint64_t>& b) {
	return AvxDividesResult{ avx::mask8bit(Rem256_64(a, b)) };
}
bool divides(double n, uint64_t div) {
	return n == div * std::floor(n / div);
}
bool divides(uint64_t n, uint64_t div) {
	return n % div == 0;
}

namespace TrialDivisionOption {
	enum {
		UseAvx = 1,
		UseTable = 2,
		UseWheel = 4,
		UseBound = 8,
		UseAll = UseAvx | UseTable | UseWheel
	};
}

#define IsAvx ((Options & TrialDivisionOption::UseAvx) != 0)
#define IsTable ((Options & TrialDivisionOption::UseTable) != 0)
#define IsWheel ((Options & TrialDivisionOption::UseWheel) != 0)
#define HasBound ((Options & TrialDivisionOption::UseBound) != 0)

template<typename T, int Options = TrialDivisionOption::UseAll>
uint64_t trialDivisionSingle(T n, uint64_t& biggestFactorChecked, std::size_t& currentTableEntryId, std::size_t& wheelIndex, const T* primes, uint64_t bound) {
	const double* mulInverses = MultiplicativeInversesDouble.data();
	if (IsTable && currentTableEntryId > 999'996)
		return trialDivisionSingle<T, Options & ~TrialDivisionOption::UseTable>(n, biggestFactorChecked, currentTableEntryId, wheelIndex, primes, bound);
	if (biggestFactorChecked < 7) {
		if (divides(n, 2)) return 2;
		if (divides(n, 3)) return 3;
		if (divides(n, 5)) return 5;
		if (divides(n, 7)) return 7;
		biggestFactorChecked = 7;
	}
	T sqrtN = (T)std::ceil(std::sqrt(n));
	if constexpr (IsAvx) {
		auto nVec = avx::setElements<T>(n);
		[[maybe_unused]] auto wheelInc = avx::setElements<T>(12, 12, 16, 14);
		bool wheelIncState = 0;
		auto factorVec = IsTable
			? avx::load(&primes[currentTableEntryId])
			: IsWheel
				? avx::setElements<T>((T)biggestFactorChecked, (T)biggestFactorChecked+4, (T)biggestFactorChecked+6, (T)biggestFactorChecked+10)
				: avx::setElements<T>((T)biggestFactorChecked, (T)biggestFactorChecked+1, (T)biggestFactorChecked+2, (T)biggestFactorChecked+3);
		[[maybe_unused]] auto factorVecInverses = avx::load(&mulInverses[currentTableEntryId]);
		
		// NOTE: Possibly one could gain significant speed-up by unrolling this loop
		while ((IsTable && primes[currentTableEntryId] <= sqrtN) || avx::getElement<0>(factorVec) <= sqrtN) {
			AvxDividesResult divisors;
			if constexpr (IsTable && std::is_same_v<T, double>) {
				divisors = divides(nVec, factorVec, factorVecInverses);
			} else {
				divisors = divides(nVec, factorVec);
			}
			if (divisors) {
				biggestFactorChecked = (uint64_t)avx::getElement<0>(factorVec);
				if (divisors.template get<0>()) return (uint64_t)avx::getElement<0>(factorVec);
				if (divisors.template get<1>()) return (uint64_t)avx::getElement<1>(factorVec);
				if (divisors.template get<2>()) return (uint64_t)avx::getElement<2>(factorVec);
				if (divisors.template get<3>()) return (uint64_t)avx::getElement<3>(factorVec);
			} else {
				if constexpr (IsTable) {
					currentTableEntryId += 4;
					if (currentTableEntryId + 4 > 1'000'000) {
						biggestFactorChecked = 15'485'827;
						return trialDivisionSingle<T, Options>(n, biggestFactorChecked, currentTableEntryId, wheelIndex, primes, bound);
					}
					if (HasBound) biggestFactorChecked = static_cast<uint64_t>(primes[currentTableEntryId-1]);
					factorVec = avx::load(primes + currentTableEntryId);
					factorVecInverses = avx::load(mulInverses + currentTableEntryId);
				} else if constexpr (IsWheel) {
					factorVec = factorVec + wheelInc;
					if (wheelIncState) {
						wheelInc = avx::setElements<T>(12, 12, 16, 14);
						if (HasBound) biggestFactorChecked += 18;
					} else {
						wheelInc = avx::setElements<T>(18, 18, 14, 16);
						if (HasBound) biggestFactorChecked += 12;
					}
					wheelIncState = !wheelIncState;
				} else {
					factorVec = factorVec + avx::setElements<T>(4);
					if (HasBound) biggestFactorChecked += 4;
				}
				if (HasBound) {
					if (biggestFactorChecked >= bound) {
						biggestFactorChecked = (uint64_t)avx::getElement<0>(factorVec); // TODO: its not fully accurate
						return (uint64_t)n; // n is coprime to all numbers below given bound
					}
				}
			}
		}
	} else {
		biggestFactorChecked = IsTable ? primes[currentTableEntryId] : biggestFactorChecked;
		std::array<T, 8> wheelInc = { 4, 2, 4, 2, 4, 6, 2, 6 };
		while (biggestFactorChecked <= sqrtN) {
			auto div = divides(n, biggestFactorChecked);
			if (div)
				return biggestFactorChecked;
			if constexpr (IsTable) {
				currentTableEntryId += 1;
				if (currentTableEntryId >= 1'000'000) {
					biggestFactorChecked = 15'485'827;
					return trialDivisionSingle<T, Options & ~TrialDivisionOption::UseTable>(n, biggestFactorChecked, currentTableEntryId, wheelIndex, primes, bound);
				}
				biggestFactorChecked = primes[currentTableEntryId];
			} else if constexpr (IsWheel) {
				biggestFactorChecked += wheelInc[wheelIndex];
				wheelIndex = (wheelIndex < wheelInc.size() - 1) ? wheelIndex + 1 : 0;
			} else {
				biggestFactorChecked += 1;
			}
			if (biggestFactorChecked >= bound) {
				return (uint64_t)n; // n is coprime to all numbers below given bound
			}
		}
	}
	return (uint64_t)n; // n is prime
}

template<typename T> T trialDivision(const T& n, std::size_t& currentTableEntryId, uint64_t bound) {
	debugAssert(bound < Primes_1_000_000.back()); // TODO: make it precise (and probably in release also?) depending on number of limbs of n.
	if (n[0] % 2 == 0) return T{ 2 };
	while (Primes_1_000_000[currentTableEntryId] <= bound) {
		if (mod(n, Primes_1_000_000[currentTableEntryId], (int)currentTableEntryId) == 0) {
			return T{ Primes_1_000_000[currentTableEntryId] };
		}
		currentTableEntryId += 1;
	}
	return n;
}

template<int Options=TrialDivisionOption::UseAll> 
uint64_t trialDivisionSingle(uint64_t n, uint64_t& biggestFactorChecked, std::size_t& currentTableEntryId, std::size_t& wheelIndex, uint64_t bound=MaxU64) {
	if (n <= MaxIntegerDouble) {
		if (bound == MaxU64) {
			return trialDivisionSingle<double, Options & ~TrialDivisionOption::UseBound>((double)n, biggestFactorChecked, currentTableEntryId, wheelIndex, &PrimesDouble64_1_000_000[0], bound);
		} else {
			return trialDivisionSingle<double, Options | TrialDivisionOption::UseBound>((double)n, biggestFactorChecked, currentTableEntryId, wheelIndex, &PrimesDouble64_1_000_000[0], bound);
		}
	} else {
		if (currentTableEntryId < 1'000'000) {
			auto newBound = std::min<uint64_t>(bound, Primes_1_000_000[Primes_1_000_000.size() - 2]);
			auto res = trialDivision(BigIntFixedSize<1>{n}, currentTableEntryId, newBound);
			biggestFactorChecked = Primes_1_000_000[currentTableEntryId];
			if (res[0] != n || bound <= newBound) {
				return res[0];
			}
		}
		if (bound == MaxU64) {
			return trialDivisionSingle<uint64_t, Options & ~TrialDivisionOption::UseBound>(n, biggestFactorChecked, currentTableEntryId, wheelIndex, &Primes64_1_000_000[0], bound);
		} else {
			return trialDivisionSingle<uint64_t, Options | TrialDivisionOption::UseBound>(n, biggestFactorChecked, currentTableEntryId, wheelIndex, &Primes64_1_000_000[0], bound);
		}
	}
}

//template<int Options=TrialDivisionOption::UseAll> 
//std::vector<uint64_t> trialDivisionSingleAll(uint64_t n) {
//	std::vector<uint64_t> factors;
//	uint64_t biggestFactorChecked = 2;
//	std::size_t currentTableEntryId = 0;
//	std::size_t wheelIndex = 0;
//	while (true) {
//		factors.emplace_back(trialDivisionSingle<Options>(n, biggestFactorChecked, currentTableEntryId, wheelIndex));
//		if (n == factors.back())
//			return factors;
//		n /= factors.back();
//	}
//}


// for users
uint64_t trialDivision(uint64_t n, uint64_t& biggestFactorChecked, std::size_t& currentTableEntryId, std::size_t& wheelIndex, uint64_t bound=MaxU64) {
	return trialDivisionSingle(n, biggestFactorChecked, currentTableEntryId, wheelIndex, bound);
}
uint64_t trialDivision(uint64_t n, uint64_t bound=MaxU64) {
	uint64_t biggestFactorChecked = 2;
	std::size_t currentTableEntryId = 0;
	std::size_t wheelIndex = 0;
	return trialDivision(n, biggestFactorChecked, currentTableEntryId, wheelIndex, bound);
}
//std::vector<uint64_t> trialDivisionAll(uint64_t n) {
//	return trialDivisionSingleAll(n);
//}

template<typename T> T trialDivision(const T& n, uint64_t bound) {
	std::size_t currentTableEntryId = 0;
	return trialDivision(n, currentTableEntryId, bound);
}

BigInt trialDivision(const BigInt& n, uint64_t bound=MaxU64) {
	return n.visitNoInit([bound](auto&& a) {
		if constexpr (BigInt::IsType<decltype(a), BigIntFixedSize<1>>()) {
			return BigInt{ trialDivision(a.mod[0], bound) };
		} else {
			debugAssert(bound != MaxU64, "bound cannot stay unset for input numbers above 64 bits");
			return BigInt{ trialDivision(a.mod, bound) };
		}
	});
}
BigInt trialDivision(const BigInt& n, uint64_t& biggestFactorChecked, std::size_t& currentTableEntryId, std::size_t& wheelIndex, uint64_t bound=MaxU64) {
	return n.visitNoInit([&](auto&& a) { return BigInt{ trialDivision(a.mod[0], biggestFactorChecked, currentTableEntryId, wheelIndex, bound) }; });
}


#undef IsAvx
#undef IsTable
#undef IsConst
#undef IsWheel
