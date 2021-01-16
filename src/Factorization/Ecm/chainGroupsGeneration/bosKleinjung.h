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


std::vector<NumberWithCost> bosKleinjung(std::function<double(uint64_t)> costFunction) {
    //generateNumbersWithSparseBinaryRepresentations("TEST_FILE_60_7.dat", 60, 7);
    //generateSpraseB1SmoothNumbers("TEST_FILE_60_7.dat", "SMOOTH_TEST_60_7_256.dat", 256);
    auto numbers = readB1SmoothNumbers("SMOOTH_TEST_60_7_256.dat");
    int B1 = 256;
    std::sort(numbers.begin(), numbers.end());
    if (numbers[0] == 1) { // should always be true, but checking to be sure.
        numbers.erase(numbers.begin());
    }

    auto numbersWithCost = createNumbersWithCosts(numbers, B1, costFunction);
    std::sort(numbersWithCost.begin(), numbersWithCost.end(), [](auto& a, auto& b) {
        return a.cost < b.cost;
    });

    auto primes = createAllB1SmoothPrimes(B1);
    while (primes[0] == 2) {
        primes.erase(primes.begin());
    }
    std::vector<NumberWithCost> result;
    std::size_t i = 0;
    while (primes.size() > 0 && i < numbersWithCost.size()) {
        if (contains(primes, numbersWithCost[i].factors)) {
            subtract(primes, numbersWithCost[i].factors);
            result.emplace_back(numbersWithCost[i]);
        }
        ++i;
    }
    if (primes.size() > 0) {
        uint64_t remainingPrimesProduct = 1;
        for (auto& prime : primes) {
            remainingPrimesProduct *= prime;
        }
        result.push_back(createNumberWithCost(remainingPrimesProduct, B1, costFunction));
    }
    result.push_back(createNumberWithCost(1 << (sizeInBits(B1) - 1), B1, costFunction));

    return result;
}
