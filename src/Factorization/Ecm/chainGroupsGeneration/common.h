#pragma once

#include <cstdint>
#include <vector>
#include <functional>
#include <iomanip>
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
    NumberWithCost() {}
    NumberWithCost(uint64_t number, double cost, const std::vector<uint64_t>& factors) : number(number), cost(cost), factors(factors) {}

    uint64_t number;
    double cost;
    std::vector<uint64_t> factors;
};

NumberWithCost createNumberWithCost(uint64_t n, uint64_t B1, std::function<double(uint64_t)> costFunction) {
    return NumberWithCost(n, costFunction(n), trialDivisionAll(n, B1));
}
std::vector<NumberWithCost> createNumbersWithCosts(const std::vector<uint64_t>& nums, uint64_t B1, std::function<double(uint64_t)> costFunction) {
    std::vector<NumberWithCost> result;
    for (auto n : nums) {
        result.push_back(createNumberWithCost(n, B1, costFunction));
    }
    return result;
}
double naf2AddCountCost(uint64_t n) {
    return nafCost(wnaf<2>(n), 0, 1, 0, 1);
}
std::function<double(uint64_t)> costDividedBySize(std::function<double(uint64_t)> costFunction) {
    return [costFunction](uint64_t n) {
        return costFunction(n) / sizeInBits(n);
    };
}
double totalCost(const std::vector<NumberWithCost>& numbers, std::function<double(uint64_t)> costFunction) {
    double result = 0;
    for (auto& n : numbers) {
        result += costFunction(n.number);
    }
    return result;
}
double totalCost(const std::vector<NumberWithCost>& numbers) {
    double result = 0;
    for (auto& n : numbers) {
        result += n.cost;
    }
    return result;
}

void printPartition(const std::vector<NumberWithCost>& numbers) {
    writeln("results:");
    for (auto& n : numbers) {
        std::cout << std::setw(20) << n.number << " ";
        std::cout << std::setw(50) << toString(n.factors) << " ";
        std::cout << std::endl;
    }
    writeln("totalCost = ", totalCost(numbers));
}

bool contains(const std::vector<uint64_t>& primes, const std::vector<uint64_t>& factors) {
    if (primes.size() == 0) {
        return false;
    }
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
