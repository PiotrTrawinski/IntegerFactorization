#pragma once

#include "common.h"
#include <cstdint>

void dixonLenstraRecurs(uint64_t B1, int tupleSize, int i, uint64_t tupleProduct, const std::vector<uint64_t>& primes, std::vector<NumberWithCost>& result) {
    if (tupleSize == 0) {
        result.emplace_back(tupleProduct, B1);
    } else {
        uint64_t baseTupleProduct = tupleProduct;
        for (; i < (int)primes.size(); ++i) {
            tupleProduct = baseTupleProduct * primes[i];
            dixonLenstraRecurs(B1, tupleSize - 1, i + 1, tupleProduct, primes, result);
        }
    }
}

void dixonLenstra(uint64_t B1, int tupleSize) {
    auto primes = createAllB1SmoothPrimes(B1);
    while (primes[0] == 2) {
        primes.erase(primes.begin());
    }

    for (int i = 0; i < tupleSize - 1; ++i) {
        primes.emplace_back(1);
    }
    std::vector<NumberWithCost> numbersWithCost;
    dixonLenstraRecurs(B1, tupleSize, 0, 1, primes, numbersWithCost);
    std::sort(numbersWithCost.begin(), numbersWithCost.end(), [](auto& a, auto& b) {
        return a.cost() < b.cost();
    });

    std::vector<NumberWithCost> result;
    //for (auto prime : primes) {
    //    result.emplace_back(prime, B1);
    //}
    std::size_t i = 0;
    while (primes.size() > 0 && i < numbersWithCost.size()) {
        if (contains(primes, numbersWithCost[i].factors)) {
            subtract(primes, numbersWithCost[i].factors);
            result.emplace_back(numbersWithCost[i]);
        }
        ++i;
    }

    writeln("results:");
    for (auto& res : result) {
        std::cout << std::setw(20) << res.number << " ";
        std::cout << std::setw(50) << toString(res.factors) << " ";
        std::cout << nafToString(res.nafForm) << '\n';
    }
    
    int dblCount = 0;
    int addCount = 0;
    for (auto& res : result) {
        dblCount += res.dblCount;
        addCount += res.addCount;
    }
    writeln("dblCount = ", dblCount);
    writeln("addCount = ", addCount);
    writeln("totalCost = ", dblCount + addCount);
}
