#pragma once

#include "common.h"
#include "../../../Utility/generalUtils.h"
#include <cstdint>

void dixonLenstraGenerateTuples(int tupleSize, int i, uint64_t tupleProduct, const std::vector<uint64_t>& primes, std::vector<uint64_t>& result) {
    if (tupleSize == 0) {
        result.emplace_back(tupleProduct);
    } else {
        uint64_t baseTupleProduct = tupleProduct;
        for (; i < (int)primes.size(); ++i) {
            tupleProduct = baseTupleProduct * primes[i];
            dixonLenstraGenerateTuples(tupleSize - 1, i + 1, tupleProduct, primes, result);
        }
    }
}

std::vector<NumberWithCost> dixonLenstra(uint64_t B1, int tupleSize, std::function<double(uint64_t)> costFunction) {
    // get multiset of all prime factors of B1 besides '2'
    auto primes = createAllB1SmoothPrimes(B1);
    while (primes[0] == 2) {
        primes.erase(primes.begin());
    }

    // add auxilary '1' "primes" to act as padding.
    // we want not only tuples with exactly "tupleSize" elements, but all tuples sizes in range [1, tupleSize].
    for (int i = 0; i < tupleSize - 1; ++i) {
        primes.emplace_back(1);
    }
    std::vector<uint64_t> numbers;
    dixonLenstraGenerateTuples(tupleSize, 0, 1, primes, numbers);

    // calculate cost of all tuples
    auto numbersWithCosts = createNumbersWithCosts(numbers, B1, costFunction);
    std::sort(numbersWithCosts.begin(), numbersWithCosts.end(), [](auto& a, auto& b) {
        return a.cost < b.cost;
    });

    // greedly choose best tuples
    std::vector<NumberWithCost> result;
    std::size_t i = 0;
    while (primes.size() > 0 && i < numbersWithCosts.size()) {
        if (contains(primes, numbersWithCosts[i].factors)) {
            subtract(primes, numbersWithCosts[i].factors);
            result.push_back(numbersWithCosts[i]);
        }
        ++i;
    }

    // add all powers of 2 (which we previously exluded)
    result.push_back(createNumberWithCost(1 << (sizeInBits(B1)-1), B1, costFunction)); 

    return result;
}
