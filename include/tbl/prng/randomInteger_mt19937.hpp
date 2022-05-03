#ifndef TBLIB_RANDOMINTEGER_MT19937_HPP
#define TBLIB_RANDOMINTEGER_MT19937_HPP

#include <random>
#include <ctime>
#include <cassert>

namespace tbl {

//static std::random_device rd;
//static std::mt19937_64 ctx{rd()};
static std::mt19937_64 ctx(time(0));

// Return a random number at most N bits long.
uint64_t getRandomInteger(size_t N) {
    assert(N > 0 && N <= 64);
    uint64_t mask = ~(uint64_t)0;
    if (N != 64) mask = ((uint64_t)1 << N) - 1;
    return ctx() & mask;
}

// Return a random number *n* so that *a <= n < b*.
uint64_t getRandomRange(uint64_t a, uint64_t b) {
    assert(a < b);
    return (ctx() % (b - a)) + a;
}

// Return a random number with exactly N-bits,
// i.e. a random number between 2**(N-1) and (2**N)-1.
uint64_t getRandomNBitInteger(size_t N) {
    assert(N > 0 && N <= 64);
    uint64_t mask = ~(uint64_t)0;
    if (N != 64) mask = ((uint64_t)1 << N) - 1;
    uint64_t r = ctx();
    if (r & (1 << (N - 1))) {
        return r & mask;
    } else {
        return (~r) & mask;
    }
}

}

#endif //TBLIB_RANDOMINTEGER_MT19937_HPP