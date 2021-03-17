#pragma once

#include "Ecm/ecm.h"
#include "TrialDivision.h"
#include "PollardRho.h"
#include "Pminus1.h"
#include "../PrimalityTesting/isProbablyPrime.h"

std::vector<std::pair<uint64_t, uint64_t>> B1_Curve_Pairs = {
	{     1629,    10}, // 40
	{     4537,    10}, // 45
	{    12322,     9}, // 50
	{    21905,    21}, // 60
	{    32918,    66}, // 70
	{   183849,   219}, // 90
	{  3071166,   649}, // 120
	{  9267681,  2399}, // 140
	{ 35158748,  6076}, // 160
	{491130495, 29584}, // 200
};

std::vector<BigInt> factor(BigInt n, bool writeDebug=false) {
	if (writeDebug) writeln("Started factorization of ", n);
	std::vector<BigInt> factors;
	BigInt factor;
	while (true) {
		if (isProbablyPrime(n, factor, writeDebug)) {
			factors.emplace_back(n);
			if (writeDebug) writeln("Remaining number is probably prime\nFinal factorization = ", factors);
			return factors;
		}
		if (!factor.isOne()) {
			factors.emplace_back(factor);
			n /= factor;
			if (writeDebug) writeln("Found factor=", factor, ". Remaining number=", n);
			continue;
		}
		
		PollardRhoParams pollardRhoParams(0, 1'000'000);
		if (writeDebug) writeln("Running Pollard Rho with max iteration count=", pollardRhoParams.maxIterCount, "...");
		factor = pollardRhoBrent(pollardRhoParams, n);
		for (int i = 0; factor.isOne() && i < B1_Curve_Pairs.size(); ++i) {
			auto B1 = B1_Curve_Pairs[i].first;
			auto curveCount = B1_Curve_Pairs[i].second;
			EcmContext ecmContext(EcmMulMethod::Naf, EcmMulCascadeMethod::Seperate, B1, B1, curveCount);

			if (writeDebug) writeln("Running P-1 with B1=", ecmContext.B1, "; B2=", ecmContext.B2, "...");
			factor = pMinus1(n, ecmContext.B1, ecmContext.B2);
			if (!factor.isOne())
				break;

			if (writeDebug) writeln("Running ECM with B1=", ecmContext.B1, "; B2=", ecmContext.B2, "; L=", ecmContext.curveCount, "...");
			factor = ecm(ecmContext, EllipticCurveForm::TwistedEdwards, n);
		}
		factors.emplace_back(factor);
		n /= factor;
		if (writeDebug) writeln("Found factor=", factor, ". Remaining number=", n);
	}
}
