#pragma once
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include "../Factorization/TrialDivision.h"

void createArrayFirstNBitPrimesInFile(const std::string& fileName, int primeCount, int maxBitCount) {
	std::ofstream file(fileName, std::ios::binary);
	std::vector<uint64_t> primes(primeCount, 0);
	file.write((const char*)primes.data(), primes.size() * sizeof(uint64_t)); // bitCount = 0
	file.write((const char*)primes.data(), primes.size() * sizeof(uint64_t)); // bitCount = 1
	for (int bitCount = 2; bitCount <= maxBitCount; ++bitCount) {
		int i = 0;
		while (i < primeCount) {
			uint64_t end = (bitCount == 64) ? std::numeric_limits<uint64_t>::max() : (1ull << bitCount);
			for (uint64_t n = (1ull << (bitCount - 1)); i < primeCount && n < end; ++n) {
				if (trialDivision(n) == n) {
					primes[i++] = n;
					std::cout << "\rbitCount = " << bitCount << " i = " << i << "   ";
				}
			}
		}
		file.write((const char*)primes.data(), primes.size() * sizeof(uint64_t));
	}
	std::cout << "\rDone                               \n";
}
std::vector<std::vector<uint64_t>> loadArrayFirstNBitPrimesFromFile(const std::string& fileName, int primeCount) {
	std::vector<std::vector<uint64_t>> primes;
	std::vector<uint64_t> nBitPrimes(primeCount);
	
	std::ifstream file(fileName, std::ios::binary);
	while (file.read((char*)nBitPrimes.data(), nBitPrimes.size() * sizeof(uint64_t))) {
		primes.emplace_back(nBitPrimes);
	}

	return primes;
}

static const auto Primes_first_100_per_bit = loadArrayFirstNBitPrimesFromFile("100_first_primes_up_to_64_bits.dat", 100);

