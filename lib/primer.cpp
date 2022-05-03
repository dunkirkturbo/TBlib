#include "tbl/helper/primer.hpp"

namespace tbl {

namespace helper {

//static uint64_t powmod(uint64_t a, uint64_t b, uint64_t n) {
//    typedef unsigned int uint128_t __attribute__((mode(TI)));
//    uint64_t r = 1;
//    while (b != 0) {
//        if (b & 1) {
//            r = ((uint128_t)r * a) % n;
//        }
//        a = ((uint128_t)a * a) % n;
//        b >>= 1;
//    }
//    return r;
//}

// Miller Rabin Algorithm
bool isPrime(uint64_t n, int iteration) {
    if (n < 2 || (n > 2 && (n & 1) == 0)) {
        return false;
    }
    uint64_t ord_multiple = n - 1, b;
    while ((ord_multiple & 1) == 0) {
        ord_multiple >>= 1;
    }
    for (int i = 0; i < iteration; i++) {
        uint64_t a = getRandomRange(1, n);
        // or hard-coding a \in [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
        b = ord_multiple;
        uint64_t r = powmod(a, b, n);
        while (b != n - 1 && r != 1 && r != n - 1) {
            r = powmod(r, 2, n);
            b <<= 1;
        }
        if (r != n - 1 && (b & 1) == 0) {
            return false;
        }
    }
    return true;
}

// Return a random N-bit prime number.
uint64_t getPrime(size_t N) {
    uint64_t n = getRandomNBitInteger(N);
    while (!isPrime(n)) {
        n = getRandomNBitInteger(N);
    }
    return n;
}

uint64_t pollardRho(uint64_t n) {
    if (n == 1) return n;
    if ((n & 1) == 0) return 2;
    uint64_t x = getRandomRange(2, n);
    uint64_t y = x;
    uint64_t c = getRandomRange(1, n);
    uint64_t d = 1;
    while (d == 1) {
        x = (powmod(x, 2, n) + c) % n;
        y = (powmod((powmod(y, 2, n) + c) % n, 2, n) + c) % n;
        uint64_t abs = x > y ? x - y : y - x;
        d = gcd<uint64_t>(n, abs);
        if (d == n) {
            return pollardRho(n);
        }
    }
    return d;
}

void factor(uint64_t n, std::unordered_map<uint64_t, int> &res) {
    if (n == 1) return;
    if (isPrime(n)) {
        res[n]++;
        return;
    }
    uint64_t _n = pollardRho(n);
    factor(_n, res);
    factor(n / _n, res);
}

}

}