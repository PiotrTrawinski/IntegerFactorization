#pragma once
#include <cstdint>
#include <vector>
#include <fstream>
#include "../Utility/generalUtils.h"
#include "../PrecomputedTables/primeTable.h"
#include "../BigInt/64bitIntrinsics.h"

std::vector<uint32_t> createPowerModTable(int maxValue, int sizePerPrime) {
	std::vector<uint32_t> result;
	result.emplace_back(0);
	result.emplace_back(1);
	for (int i = 2; i <= maxValue; ++i) {
		uint32_t val1 = (MaxU64 % i + 1) % i;
		result.emplace_back(val1);
		uint32_t val = val1;
		for (int i = 2; i < sizePerPrime; ++i) {
			val = ((uint64_t)val * val1) % i;
			result.emplace_back(val);
		}
	}
	return result;
}
void createPowerModTableInFile(const std::string& fileName, int sizePerPrime) {
	std::ofstream file(fileName, std::ios::binary);
	for (auto m : Primes_1_000_000) {
		uint32_t val1 = (MaxU64 % m + 1) % m;
		file.write((const char*)&val1, sizeof(uint32_t));
		uint32_t val = val1;
		for (int i = 2; i < sizePerPrime; ++i) {
			val = ((uint64_t)val * val1) % m;
			file.write((const char*)&val, sizeof(uint32_t));
		}
	}
}
Array2d<uint32_t> loadPowerModTableFromFile(const std::string& fileName, int sizePerPrime) {
	Array2d<uint32_t> result((int)Primes_1_000_000.size(), sizePerPrime-1);
	
	std::ifstream file(fileName, std::ios::binary);
	file.read((char*)result.data.data(), Primes_1_000_000.size()*(sizePerPrime - 1)*sizeof(uint32_t));

	return result;
}


uint64_t getInverse(uint64_t d) {
	uint64_t l = sizeInBits(d);
	if (bitSetCount(d) == 1) {
		l -= 1;
	}
	uint64_t x[2];
	x[1] = (1ull << l) - d;
	x[0] = 0;
	auto m = div128(x[1], x[0], d);
	return m + 1;
}
std::vector<uint64_t> createMultiplicativeInverses(int maxValue) {
	std::vector<uint64_t> result;
	result.emplace_back(0);
	result.emplace_back(1);
	for (int i = 2; i <= maxValue; ++i) {
		result.emplace_back(getInverse(i));
	}
	return result;
}
void createMultiplicativeInversesInFile(const std::string& fileName) {
	std::ofstream file(fileName, std::ios::binary);
	for (uint32_t d : Primes_1_000_000) {
		uint64_t m = getInverse(d);
		file.write((const char*)&m, sizeof(uint64_t));
	}
}
void createMultiplicativeInversesDoubleInFile(const std::string& fileName) {
	std::ofstream file(fileName, std::ios::binary);
	for (uint32_t d : Primes_1_000_000) {
		double m = 1.0 / d;
		if (d * m < 1.0) {
			m += m * std::numeric_limits<double>::epsilon();
		}
		file.write((const char*)&m, sizeof(double));
	}
}
std::vector<uint64_t> loadMultiplicativeInversesTableFromFile(const std::string& fileName) {
	std::ifstream file(fileName, std::ios::binary);
	std::vector<uint64_t> result;
	result.resize(1'000'000);
	file.read((char*)result.data(), result.size() * sizeof(uint64_t));
	return result;
}
std::vector<double> loadMultiplicativeInversesDoubleTableFromFile(const std::string& fileName) {
	std::ifstream file(fileName, std::ios::binary);
	std::vector<double> result;
	result.resize(1'000'000);
	file.read((char*)result.data(), result.size() * sizeof(double));
	return result;
}

static const auto PowerModTable = loadPowerModTableFromFile("Power_Mod_Table_8.dat", 8);
static const auto MultiplicativeInverses = loadMultiplicativeInversesTableFromFile("Multiplicative_Inverses.dat");
static const auto MultiplicativeInversesDouble = loadMultiplicativeInversesDoubleTableFromFile("Multiplicative_Inverses_Double.dat");

static const auto PowerModTable1000 = createPowerModTable(1'000, 8);
static const auto MultiplicativeInverses1000 = createMultiplicativeInverses(1'000);
