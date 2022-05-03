#ifndef TBLIB_PRIMER_SELECT_HPP
#define TBLIB_PRIMER_SELECT_HPP

#include "tbl/params.hpp"
#include "tbl/helper/primer.hpp"

namespace tbl {

namespace helper {

static constexpr int DEFAULT_MAX_PN = 100;

template<class T>
int genPrime(T* primes, int MAX_PN = DEFAULT_MAX_PN) {
    T N = params<T>::polyDegreeN;
    T p = ((T)1 << tbl::params<T>::kModulusBitsize) - 2 * N + 1;
    int PN = 0;
    T lower_bound = (T)1 << (tbl::params<T>::kModulusBitsize - 1);
    while (PN < MAX_PN && p > lower_bound) {
        if (isPrime(p, DEFAULT_K)) {
            primes[PN] = p;
            PN++;
        }
        p -= 2 * N;
    }
    return PN;
}

}

}

#endif //TBLIB_PRIMER_SELECT_HPP
