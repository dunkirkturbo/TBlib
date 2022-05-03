#ifndef TBLIB_PRIMITIVE_ROOT_SELECT_HPP
#define TBLIB_PRIMITIVE_ROOT_SELECT_HPP

#include "tbl/params.hpp"
#include "tbl/helper/primer.hpp"
#include <unordered_map>

namespace tbl {

namespace helper {

template<class T>
void genPrimitiveRoot(T* primitive_roots) {
    int pn = tbl::params<T>::kModulusMaxN;
    std::unordered_map<uint64_t, int> res;
    for (int i = 0; i < pn; i++) {
        res.clear();
        T p = tbl::params<T>::primes[i];
        factor(p - 1, res);
        for (T g = 2; g < p; g++) {
            bool isPrimitive = true;
            for (auto factor: res) {
                if (powmod(g, (p - 1) / factor.first, p) == 1) {
                    isPrimitive = false;
                    break;
                }
            }
            if (isPrimitive) {
                primitive_roots[i] = g;
                break;
            }
        }
    }
}

// eta^{2n} = 1 mod p
template<class T>
void genEta(T* etas, T* inv_etas) {
    using signed_value_type = typename params<T>::signed_value_type;
    int pn = tbl::params<T>::kModulusMaxN;
    T primitive_roots[pn];
    genPrimitiveRoot<T>(primitive_roots);
    int n = tbl::params<T>::polyDegreeN;
    std::cout << n << std::endl;
    for (int i = 0; i < pn; i++) {
        T p = tbl::params<T>::primes[i];
        etas[i] = (T)(powmod(primitive_roots[i], (p - 1) / (2 * n), p));
        inv_etas[i] = inverse<signed_value_type>(etas[i], p);
    }
}

// omega^{n} = 1 mod p
template<class T>
void genOmega(T* etas, T* omegas, T* inv_omegas) {
    using signed_value_type = typename tbl::params<T>::signed_value_type;
    int pn = tbl::params<T>::kModulusMaxN;
    for (int i = 0; i < pn; i++) {
        omegas[i] = (T)(powmod(etas[i], 2, tbl::params<T>::primes[i]));
        inv_omegas[i] = inverse<signed_value_type>(omegas[i], tbl::params<T>::primes[i]);
    }
}

}

}

#endif //TBLIB_PRIMITIVE_ROOT_SELECT_HPP
