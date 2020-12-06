#pragma once

#include <cstdint>
#include <vector>
#include "../../../PrecomputedTables/primeTable.h"
#include "../multiplicationMethods/wnafMul.h"

std::vector<uint64_t> createAllB1SmoothPrimes(uint64_t B1) {
	std::vector<uint64_t> result;
	for (int i = 0; Primes_1_000_000[i] <= B1; ++i) {
		for (uint64_t p = Primes_1_000_000[i]; p <= B1; p *= Primes_1_000_000[i]) {
			result.emplace_back(Primes_1_000_000[i]);
		}
	}
	return result;
}

struct NumberWithCost {
	NumberWithCost(uint64_t number, uint64_t B1) : number(number) {
		auto factNum = number;
		while (factNum != 1) {
			auto factor = trialDivision(factNum, B1);
			factors.emplace_back(factor);
			factNum /= factor;
		}
		double avgFactorSize = 0;
		for (auto factor : factors) {
			avgFactorSize += factor;
		}
		avgFactorSize /= factors.size();
		avgFactorSize = std::pow(avgFactorSize, 1.0 / 8.0);
		//avgFactorSize = std::sqrt(avgFactorSize);

		auto [bestWNaf, bestNafForm] = getBestWNaf(number, [&](auto& naf) {
			auto [dblCount, addCount] = nafDblAddCounts(naf);
			auto cost = (double)addCount / dblCount;
			return cost;
		});
		nafForm = bestNafForm;
		auto [dblCountX, addCountX] = nafDblAddCounts(nafForm);
		dblCount = dblCountX;
		addCount = addCountX;
		cost_ = (1 / avgFactorSize)*(double)addCount / dblCount;
	}
	uint64_t number;
	int dblCount;
	int addCount;
	double cost_;
	StackVector<int8_t, 64> nafForm;
	std::string method;
	std::vector<uint64_t> factors;

	double cost() {
		return cost_;// (double)setBitCount / bitCount;
	}
};

std::vector<NumberWithCost> convertToNumbersWithCost(const std::vector<uint64_t>& numbers, uint64_t B1) {
	std::vector<NumberWithCost> result;
	for (auto number : numbers) {
		result.emplace_back(number, B1);
	}
	return result;
}

bool contains(const std::vector<uint64_t>& primes, const std::vector<uint64_t>& factors) {
	std::size_t i = 0;
	std::size_t j = 0;
	for (;;) {
		if (primes[i] == factors[j]) {
			++j;
			if (j == factors.size()) {
				break;
			}
		} else if (primes[i] > factors[j]) {
			return false;
		}
		++i;
		if (i == primes.size()) {
			return false;
		}
	}
	return true;
}
void subtract(std::vector<uint64_t>& primes, const std::vector<uint64_t>& factors) {
	for (std::size_t i = 0, j = 0; j < factors.size(); ) {
		if (primes[i] == factors[j]) {
			++j;
			primes.erase(primes.begin() + i);
		} else {
			++i;
		}
	}
}
