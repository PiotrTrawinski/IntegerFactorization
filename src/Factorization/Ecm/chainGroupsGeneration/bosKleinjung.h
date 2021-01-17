#pragma once

#include "../../../Utility/generalUtils.h"
#include "../../../Utility/bitManipulation.h"
#include "../../../PrecomputedTables/primeTable.h"
#include "../../../BigInt/include.h"
#include "../multiplicationMethods/wnafMul.h"
#include "../../../Utility/smoothNumberUtils.h"
#include "common.h"
#include <iostream>
#include <string>
#include <iomanip>

std::string nafToString(const StackVector<int8_t, 64>& vec) {
    std::string str;
    for (int i = 0; i < vec.size(); ++i) {
        str += (vec[i] >= 0 ? ' ' : '-');
        str += abs(vec[i]) + '0';
        str += ' ';
    }
    return str;
}

void nextNumberSparse(std::string& n) {
    int i = 0;
    // skip first group of 0's
    while (n[i] == '0') {
        i += 1;
    }

    // skip first group of 1's (and save as oneCount)
    int firstOneIndex = i;
    while (n[i] == '1' && n[i+1] == '0') {
        i += 2;
    }
    int oneCount = (i - firstOneIndex) / 2;

    // set first '0' after group of '1's to '1'.
    n[--i] = '1';

    // clear all bits below the new '1'
    for (int j = 0; j < i; ++j) {
        n[j] = '0';
    }

    // set oneCount-1 lower bits to '1'
    for (int j = 0; j < oneCount - 1; ++j) {
        n[j*2] = '1';
    }
}
template<typename T, int Capacity> struct BufferedWriter {
    std::vector<T> data;
    int size;
    std::ofstream outStream;
    BufferedWriter(const std::string& fileName) : data(Capacity), size(0), outStream(fileName, std::ios::binary) {}
    ~BufferedWriter() {
        flush();
    }

    void flush() {
        outStream.write((const char*)data.data(), size * sizeof(T));
        size = 0;
    }
    void write(const T& value) {
        if (size >= Capacity) {
            flush();
        }
        data[size++] = value;
    }
};
template<typename T, int Capacity> struct BufferedReader {
    std::vector<T> data;
    int size;
    int index;
    uint64_t remainingSizeInFile;
    uint64_t fileSize;
    std::ifstream inStream;
    BufferedReader(const std::string& fileName) : data(Capacity), size(0), index(0), inStream(fileName, std::ios::binary) {
        inStream.read((char*)&fileSize, sizeof(uint64_t));
        remainingSizeInFile = fileSize;
    }

    bool read(T& value) {
        if (index >= size) {
            if (remainingSizeInFile == 0) {
                return false;
            }
            size = std::min<int>((int)remainingSizeInFile, Capacity);
            inStream.read((char*)data.data(), size * sizeof(uint64_t));
            remainingSizeInFile -= size;
            index = 0;
        }
        value = data[index++];
        return true;
    }
};
void generateNumbersWithSparseBinaryRepresentations(const std::string& fileName, int doubleCount, int maxAddCount) {
    BufferedWriter<uint64_t, 4096> outFile(fileName);
    outFile.write(0);
    int numberCount = 0;
    for (int addCount = 1; addCount <= maxAddCount; ++addCount) {
        std::string number(doubleCount + 3, '0');
        number[doubleCount - 1] = '1';
        for (int i = 0; i < addCount-1; ++i) {
            number[2*i] = '1';
        }
        while (number[doubleCount] != '1') {
            uint64_t intNum = std::stoull(number.substr(0, doubleCount), nullptr, 2);
            outFile.write(intNum);
            numberCount += 1;
            nextNumberSparse(number);
        }
    }
    outFile.flush();
    outFile.outStream.seekp(std::ios::beg);
    outFile.write(numberCount);
}

void generateSpraseB1SmoothNumbers(const std::string& inFileName, const std::string& outFileName, uint64_t B1) {
    BufferedReader<uint64_t, 4096> inFile(inFileName);
    BufferedWriter<uint64_t, 4096> outFile(outFileName);

    uint64_t number;
    std::array<int64_t, 64> diffs;
    int count = 0;
    while (inFile.read(number)) {
        if (count % 10'000 == 0) {
            writef("\rcount = %/% (%)        ", count, inFile.fileSize, (double)count/inFile.fileSize);
        }
        count += 1;
        int diffsCount = 0;
        uint64_t val = number;
        for (int j = 0; val != 1; ++j) {
            if (val & 1) {
                diffs[diffsCount++] = 1ll << j;
            }
            val >>= 1;
        }
        for (uint64_t a = 0; a < (1ull << diffsCount); ++a) {
            val = mostSignificantBit(number);
            for (int i = 0; i < diffsCount; ++i) {
                val += (2ull*bit(a, i)-1) * diffs[i];
            }
            if (val % 2 == 0) {
                continue;
            }
            if (isB1PowerSmooth(val, B1)) {
                outFile.write(val);
            }
        }
    }
    writeln();
}

std::vector<uint64_t> readB1SmoothNumbers(const std::string& inFileName) {
    std::ifstream inFile(inFileName, std::ios::binary);
    std::vector<uint64_t> result;
    uint64_t number;
    while (inFile.read((char*)&number, sizeof(uint64_t))) {
        result.emplace_back(number);
    }
    return result;
}

void setRatios(std::vector<double>& ratios, const std::vector<uint64_t>& factors) {
    for (auto& ratio : ratios) {
        ratio = 0;
    }
    for (auto& factor : factors) {
        ratios[sizeInBits(factor)-1] += 1;
    }
    for (auto& ratio : ratios) {
        ratio /= factors.size();
    }
}

struct NumberBosKleinjung {
    NumberBosKleinjung() {}
    NumberBosKleinjung(uint64_t number, double cost, const std::vector<uint64_t>& factors, int ratioSize) : number(number), cost(cost), factors(factors) {
        ratios.resize(ratioSize, 0);
        setRatios(ratios, factors);
    }
    void calculateScore(const std::vector<double>& primeRatios) {
        score = 0;
        for (int i = 0; i < primeRatios.size(); ++i) {
            if (primeRatios[i] != 0) {
                score += ratios[i] / primeRatios[i];
            }
        }
    }

    uint64_t number;
    double cost;
    double score;
    std::vector<double> ratios;
    std::vector<uint64_t> factors;
};

bool canCreateNumberBosKleinjung(uint64_t n, uint64_t B1) {
    auto factors = trialDivisionAll(n, B1);
    return factors.back() <= B1;
}
NumberBosKleinjung createNumberBosKleinjung(uint64_t n, uint64_t B1, std::function<double(uint64_t)> costFunction) {
    return NumberBosKleinjung(n, costFunction(n), trialDivisionAll(n, B1), sizeInBits(B1));
}
std::vector<NumberBosKleinjung> createumbersBosKleinjung(const std::vector<uint64_t>& nums, uint64_t B1, std::function<double(uint64_t)> costFunction) {
    std::vector<NumberBosKleinjung> result;
    for (auto n : nums) {
        if (canCreateNumberBosKleinjung(n, B1)) {
            result.push_back(createNumberBosKleinjung(n, B1, costFunction));
        }
    }
    return result;
}
NumberWithCost convertToNumberWithCost(const NumberBosKleinjung& number) {
    return NumberWithCost(number.number, number.cost, number.factors);
}

void printPartitionUsingCost(const std::vector<NumberWithCost>& numbers, std::function<double(uint64_t)> costFunction) {
    writeln("results:");
    double totalCost = 0;
    for (auto& n : numbers) {
        std::cout << std::setw(20) << n.number << " ";
        std::cout << std::setw(50) << toString(n.factors) << " ";
        std::cout << std::endl;
        totalCost += costFunction(n.number);
    }
    writeln("totalCost = ", totalCost);
}

void bosKleinjungUpdate(std::vector<NumberBosKleinjung>& numbersWithCost, std::vector<double>& primeRatios, const std::vector<uint64_t>& primes) {
    setRatios(primeRatios, primes);

    numbersWithCost.erase(std::remove_if(numbersWithCost.begin(), numbersWithCost.end(), [&primes](auto& a) {
        return !contains(primes, a.factors);
    }), numbersWithCost.end());
    for (auto& number : numbersWithCost) {
        number.calculateScore(primeRatios);
    }
    std::sort(numbersWithCost.begin(), numbersWithCost.end(), [](auto& a, auto& b) {
        return a.score < b.score;
    });
}

std::vector<NumberWithCost> bosKleinjung(std::function<double(uint64_t)> costFunction, std::function<double(uint64_t)> evalFunction) {
    //generateNumbersWithSparseBinaryRepresentations("TEST_FILE_60_7.dat", 60, 7);
    //generateSpraseB1SmoothNumbers("TEST_FILE_60_7.dat", "SMOOTH_TEST_60_7_256.dat", 256);
    auto numbers = readB1SmoothNumbers("SMOOTH_TEST_60_7_256.dat");
    int B1 = 31;
    std::sort(numbers.begin(), numbers.end());
    if (numbers[0] == 1) { // should always be true, but checking to be sure.
        numbers.erase(numbers.begin());
    }

    auto numbersWithCost = createumbersBosKleinjung(numbers, B1, costFunction);
    auto primes = createAllB1SmoothPrimes(B1);
    while (primes[0] == 2) {
        primes.erase(primes.begin());
    }
    std::vector<double> primeRatios(sizeInBits(B1), 0);
    
    bosKleinjungUpdate(numbersWithCost, primeRatios, primes);

    auto [minCostIter, maxCostIter] = std::minmax_element(numbersWithCost.begin(), numbersWithCost.end(), [](auto& a, auto& b) {
        return a.cost < b.cost;
    });
    auto minCost = minCostIter->cost;
    auto maxCost = maxCostIter->cost;
    auto costIncrement = (maxCost - minCost) / 100.0;

    double skipChance = 0.50;
    int tryCount = 100;
    double bestTotalCost = std::numeric_limits<double>::max();

    std::vector<NumberWithCost> result;
    for (int j = 0; j < tryCount; ++j) {
        auto primesCopy = primes;
        auto numbersWithCostCopy = numbersWithCost;
        std::vector<NumberWithCost> potentialResult;
        for (auto costThreshold = minCost; costThreshold < maxCost; costThreshold += costIncrement) {
            bool found = false;
            do {
                found = false;
                for (int i = 0; i < numbersWithCostCopy.size(); ++i) {
                    if (numbersWithCostCopy[i].cost <= costThreshold && random() > skipChance) {
                        subtract(primesCopy, numbersWithCostCopy[i].factors);
                        potentialResult.emplace_back(convertToNumberWithCost(numbersWithCostCopy[i]));
                        bosKleinjungUpdate(numbersWithCostCopy, primeRatios, primesCopy);
                        found = true;
                        break;
                    }
                }
            } while (found);
        }
        if (primesCopy.size() > 0) {
            uint64_t remainingPrimesProduct = 1;
            for (auto& prime : primesCopy) {
                remainingPrimesProduct *= prime;
            }
            potentialResult.push_back(createNumberWithCost(remainingPrimesProduct, B1, costFunction));
        }
        potentialResult.push_back(createNumberWithCost(1 << (sizeInBits(B1) - 1), B1, costFunction));

        if (totalCost(potentialResult, evalFunction) < bestTotalCost) {
            result = potentialResult;
            bestTotalCost = totalCost(potentialResult, evalFunction);
        }
    }

    return result;
}

std::vector<NumberWithCost> bosKleinjung(std::function<double(uint64_t)> costFunction) {
    return bosKleinjung(costFunction, costFunction);
}
