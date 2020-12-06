#pragma once
#include "common.h"
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <limits>

std::vector<uint64_t> optimalPrimePartition(const std::vector<int>& multiplicities, const std::vector<uint64_t>& primes, const std::unordered_map<uint64_t, int>& costs) {
    /*
        Based on Knuth's "THE ART OF COMPUTER PROGRAMMING"
        algorithm M (Multipartitions in dereasing lexiographic order)
    */

    int m = (int)multiplicities.size();
    int n = std::accumulate(multiplicities.begin(), multiplicities.end(), 0);
    int size = n * m + 1;
    std::vector<int> c(size, 0);
    std::vector<int> u(size, 0);
    std::vector<int> v(size, 0);
    std::vector<int> f(n + 1, 0);

    std::vector<uint64_t> bestPartition;
    int bestCost = std::numeric_limits<int>::max();

    // Step M1
    for (int j = 0; j < m; ++j) {
        c[j] = j;
        u[j] = multiplicities[j];
        v[j] = multiplicities[j];
    }
    f[0] = 0;
    f[1] = m;
    int a = 0;
    int b = m;
    int lPart = 0;

    while (true) {
        while (true) {
            // Step M2
            int j = a;
            int k = b;
            int x = false;
            while (j < b) {
                u[k] = u[j] - v[j];
                if (u[k] == 0) {
                    x = true;
                } else if (!x) {
                    c[k] = c[j];
                    v[k] = std::min(v[j], u[k]);
                    x = u[k] < v[j];
                    k += 1;
                } else {
                    c[k] = c[j];
                    v[k] = u[k];
                    k += 1;
                }
                j += 1;
            }

            // Step M3
            if (k > b) {
                a = b;
                b = k;
                lPart += 1;
                f[lPart + 1] = b;
            } else {
                break;
            }
        }

        // Step M4
        int cost = 0;
        for (int i = 0; i <= lPart; ++i) {
            uint64_t number = 1;
            for (int j = f[i]; j < f[i + 1]; ++j) {
                if (v[j] > 0) {
                    for (int k = 0; k < v[j]; ++k) {
                        number *= primes[c[j]];
                    }
                }
            }
            cost += costs.at(number);
        }
        if (cost < bestCost) {
            bestPartition.clear();
            for (int i = 0; i <= lPart; ++i) {
                uint64_t number = 1;
                for (int j = f[i]; j < f[i + 1]; ++j) {
                    if (v[j] > 0) {
                        for (int k = 0; k < v[j]; ++k) {
                            number *= primes[c[j]];
                        }
                    }
                }
                bestPartition.emplace_back(number);
            }
            bestCost = cost;
            std::cout << cost << "\t";
            std::cout << bestPartition[0];
            for (std::size_t i = 1; i < bestPartition.size(); ++i) {
                std::cout << " " << bestPartition[i];
            }
            std::cout << '\n';
        }

        // Step M5
        while (true) {
            int j = b - 1;
            while (v[j] == 0) {
                j -= 1;
            }
            if (j == a && v[j] == 1) {
                // Step M6
                if (lPart == 0) {
                    return bestPartition;
                }
                lPart -= 1;
                b = a;
                a = f[lPart];
            } else {
                v[j] -= 1;
                for (int k = j + 1; k < b; ++k) {
                    v[k] = u[k];
                }
                break;
            }
        }
    }
}

bool nextMultisetSubset(const std::vector<int>& multiplicities, std::vector<int>& ar, int n) {
    bool changed = false;
    for (int i = n - 1; i >= 0 && !changed; i--) {
        if (ar[i] < multiplicities[i]) {
            ar[i]++;
            changed = true;
        } else {
            ar[i] = 0;
        }
    }
    if (!changed) {
        for (int i = 0; i < n; i++) {
            ar[i] = 0;
        }
    }
    return changed;
}

void optimalChainGroup(uint64_t B1) {
    auto primes = createAllB1SmoothPrimes(B1);
    while (primes[0] == 2) {
        primes.erase(primes.begin());
    }

    std::vector<int> primeMultiplicities;
    std::vector<uint64_t> primesUnique;
    primesUnique.emplace_back(primes[0]);
    primeMultiplicities.emplace_back(1);
    for (std::size_t i = 1; i < primes.size(); ++i) {
        if (primes[i] == primesUnique.back()) {
            primeMultiplicities.back() += 1;
        } else {
            primesUnique.emplace_back(primes[i]);
            primeMultiplicities.emplace_back(1);
        }
    }

    std::vector<int> counts(primeMultiplicities.size(), 0);
    std::unordered_map<uint64_t, int> costs;
    counts.back() = 1;
    do {
        uint64_t number = 1;
        for (std::size_t i = 0; i < counts.size(); ++i) {
            for (int j = 0; j < counts[i]; ++j) {
                uint64_t prevNumber = number;
                number *= primesUnique[i];
                if (prevNumber > number) {
                    std::cout << "ERROR\n";
                }
            }
        }
        auto [bestWNaf, bestNafForm] = getBestWNaf(number, [&](auto& naf) {
            return nafCost(naf, 8, 8, 8, 8);
		});
        costs[number] = nafCost(bestNafForm, 8, 8, 8, 8);
    } while (nextMultisetSubset(primeMultiplicities, counts, (int)counts.size()));
    auto bestPartition = optimalPrimePartition(primeMultiplicities, primesUnique, costs);
}
