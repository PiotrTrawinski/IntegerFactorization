#pragma once
#include <cstdint>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include "../BigInt/64bitIntrinsics.h"
#include "../Utility/generalUtils.h"

int getWheelNumber(int i) {
	constexpr std::array<int, 8> wheelAcc = { 0, 4, 6, 10, 12, 16, 22, 24 };
	return 7 + (i / 8) * 30 + wheelAcc[static_cast<int>(i % 8)];
}
int getWheelIndex(int n) {
	constexpr std::array<int, 30> wheelInv = { 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7 };
	return ((n - 7) / 30) * 8 + wheelInv[static_cast<int>((n - 7) % 30)];
}

std::vector<uint32_t> makePrimeArray() {
	std::vector<uint32_t> primes(1'000'000);
	primes[0] = 2;
	primes[1] = 3;
	primes[2] = 5;
	std::vector<bool> isComposite(15'485'864 * 8 / 30, false); 
	constexpr std::array<int, 8> wheelFirstValues = { 7, 11, 13, 17, 19, 23, 29, 31 };
	constexpr std::array<int, 8> wheelInc = { 4, 2, 4, 2, 4, 6, 2, 6 };
	auto sqrtSize = std::sqrt(isComposite.size());
	for (int i = 0; i < sqrtSize; ++i) {
		if (!isComposite[i]) {
			auto wheelNumber = getWheelNumber(i);
			int p = wheelNumber;

			std::array<int, 8> indexInc;
			int prev = i;
			for (int k = 0; k < 8; ++k) {
				int curr = getWheelIndex(p * wheelFirstValues[k]);
				indexInc[k] = curr - prev;
				prev = curr;
			}

			std::size_t wheelIndex = 1;
			int idx = i + indexInc[0];
			while (idx < (int)isComposite.size()) {
				isComposite[idx] = true;
				idx += indexInc[wheelIndex];
				wheelIndex = (wheelIndex + 1) % indexInc.size();
			}
		}
	}
	int j = 2;
	int value = wheelFirstValues[0];
	for (int i = 0; i < (int)isComposite.size(); ++i) {
		if (!isComposite[i]) {
			primes[++j] = value;
		}
		value += wheelInc[i % 8];
	}
	return primes;
}
template<typename T> std::vector<T> makePrimeArray(const std::vector<uint32_t>& arr) {
	std::vector<T> result(arr.size());
	for (int i = 0; i < (int)result.size(); ++i) {
		result[i] = arr[i];
	}
	return result;
}
void savePrimeArrayToFile(const std::string& fileName, const std::vector<uint32_t>& arr) {
	std::ofstream file(fileName, std::ios::binary);
	file.write((const char*)arr.data(), arr.size() * sizeof(uint32_t));
}
std::vector<uint32_t> loadPrimeArrayFromFile(const std::string& fileName, int size) {
	std::ifstream file(fileName, std::ios::binary);
	std::vector<uint32_t> primes;
	primes.resize(size);
	file.read((char*)primes.data(), primes.size() * sizeof(uint32_t));
	return primes;
}

static const auto Primes_1_000_000 = loadPrimeArrayFromFile("Primes_1_000_000.dat", 1'000'000);
static const auto Primes64_1_000_000 = makePrimeArray<uint64_t>(Primes_1_000_000);
static const auto PrimesDouble64_1_000_000 = makePrimeArray<double>(Primes_1_000_000);
