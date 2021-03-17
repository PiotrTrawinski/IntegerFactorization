#pragma once
#include "../BigInt/include.h"
#include "../PrimalityTesting/millerRabin.h"
#include "../Factorization/TrialDivision.h"

constexpr uint64_t PrimeTesting1LimbTrialDivisionThreshold = 1ull << 48;
constexpr uint64_t PrimeTestingMultiLimbTrialDivisionThreshold = 1ull << 14;

bool isProbablyPrime(const BigInt& n, BigInt& factor, bool writeDebug=false) {
	factor = 1;
	if (n.IsType<BigIntFixedSize<1>>()) {
		auto nVal = n.get<BigIntFixedSize<1>>()[0];
		if (nVal < PrimeTesting1LimbTrialDivisionThreshold) {
			if (writeDebug) writeln("Running Trial Division with bound=", PrimeTesting1LimbTrialDivisionThreshold, "...");
			auto trialResult = trialDivision(nVal);
			if (trialResult != nVal) {
				factor = trialResult;
			}
			return trialResult == nVal;
		} else {
			if (writeDebug) writeln("Running Miller-Rabin primality test...");
			return millerRabinTest(nVal);
		}
	} else {
		if (writeDebug) writeln("Running Trial Division with bound=", PrimeTestingMultiLimbTrialDivisionThreshold, "...");
		auto trialResult = trialDivision(n, PrimeTestingMultiLimbTrialDivisionThreshold);
		if (trialResult != n) {
			factor = trialResult;
			return false;
		}
		if (writeDebug) writeln("Running Miller-Rabin primality test...");
		return millerRabinTest(n);
	}
}
bool isProbablyPrime(const BigInt& n) {
	BigInt factor;
	return isProbablyPrime(n, factor);
}
