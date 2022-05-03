#ifndef TBLIB_MILLER_RABIN_HPP
#define TBLIB_MILLER_RABIN_HPP

#include "tbl/prng/randomInteger_mt19937.hpp"
#include "tbl/helper/euclidean.hpp"
#include <unordered_map>

namespace tbl {

namespace helper {

static uint64_t powmod(uint64_t a, uint64_t b, uint64_t n) {
    typedef unsigned int uint128_t __attribute__((mode(TI)));
    uint64_t r = 1;
    while (b != 0) {
        if (b & 1) {
            r = ((uint128_t) r * a) % n;
        }
        a = ((uint128_t) a * a) % n;
        b >>= 1;
    }
    return r;
}

// k, the number of rounds of testing to perform
static constexpr int DEFAULT_K = 5;

// Miller Rabin Algorithm
bool isPrime(uint64_t n, int iteration = DEFAULT_K);

// Return a random N-bit prime number.
uint64_t getPrime(size_t N);

uint64_t pollardRho(uint64_t n);

void factor(uint64_t n, std::unordered_map<uint64_t, int> &res);

}

}

#endif //TBLIB_MILLER_RABIN_HPP
