#ifndef TBLIB_OPS_HPP
#define TBLIB_OPS_HPP

#include "tbl/params.hpp"
#include "tbl/helper/euclidean.hpp"
#include <cstddef>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <tuple>

namespace tbl {

namespace ops {

template<class T>
struct reducer_Barrett {
public:
    using value_type = typename params<T>::value_type;
    using greater_value_type = typename params<T>::greater_value_type;

    // construct
    explicit reducer_Barrett(value_type modulus) {
        // k-bit modulus
        kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        kModulusBitsize = kModulusRepresentationBitsize;
        for (int i = 1; i < kModulusRepresentationBitsize; i++) {
            if (((value_type)1 << i) > modulus) {
                kModulusBitsize = i;
                break;
            }
        }
        // yield approximately 2k+? bits, but let 2k+4 <= 2*kModulusRepresentationBitsize (avoid overflow)
        if (kModulusBitsize > kModulusRepresentationBitsize - 2) {
            throw std::invalid_argument("reducer_Barrett: modulus too large");
        }
        // pre-computation
        m_modulus = modulus;
        upper_bound = (greater_value_type)1 << (2 * kModulusBitsize);
        m_mu = upper_bound / m_modulus;
    }

    value_type get_modulus() {
        return m_modulus;
    }

    value_type reduce(greater_value_type x) {
        // x should be less than 2 ** (2 * kModulusBitsize), but not checked inside for efficiency
//        assert(x < upper_bound);
        greater_value_type q = x >> (kModulusBitsize - 1);
        q *= m_mu;
        q >>= (kModulusBitsize + 1);
        value_type r = x - q * m_modulus;
        while (r >= m_modulus) {
            r -= m_modulus;
        }
        return r;
    }

    value_type mulmod(value_type x, value_type y) {
        return reduce((greater_value_type)x * y);
    }

    // simple mode
    value_type reduce_simple(greater_value_type x) {
        return x % m_modulus;
    }

     value_type mulmod_simple(value_type x, value_type y) {
        return reduce_simple((greater_value_type)x * y);
    }

private:
    value_type m_modulus, m_mu;
    greater_value_type upper_bound; // for input of `reduce` function
    unsigned int kModulusRepresentationBitsize, kModulusBitsize;
};

template<class T>
struct reducer_Montgomery {
public:
    using value_type = typename params<T>::value_type;
    using signed_value_type = typename params<T>::signed_value_type;
    using greater_value_type = typename params<T>::greater_value_type;

    // construct
    explicit reducer_Montgomery(value_type modulus) {
        // k-bit modulus
        kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        kModulusBitsize = kModulusRepresentationBitsize - 1;
        // [!] use (kModulusRepresentationBitsize - 1) here but not kModulusRepresentationBitsize
        // is for avoiding the overflow of exgcd
        R = (value_type)1 << kModulusBitsize; // define R but not use R ==> use left/right shift
        if (modulus & R) {
            throw std::invalid_argument("reducer_Montgomery: modulus too large");
        }
        upper_bound = (greater_value_type)modulus << kModulusBitsize;

        // pre-computation
        m_modulus = modulus;
        auto res = helper::exgcd<signed_value_type>(R - m_modulus, m_modulus, 1, -1);
        _R = std::get<1>(res);
        _M = std::get<2>(res);
        if (_R < 0) {
            _R += m_modulus;
            _M = R - _M;
        } else {
            _M = -_M;
        }
        // R * _R - m_modulus * _M = 1
        // 0 < _R < m_modulus
        // 0 < _M < R
//        std::cout << "reducer_Montgomery: _R = " << _R << ", _M = " << _M << std::endl;

//        assert(_R > 0 && _R < m_modulus && _M > 0 && _M < R \
//        && (greater_value_type)_R * R - (greater_value_type)_M * m_modulus == 1);
        square_mod = ((greater_value_type)R << kModulusBitsize) % m_modulus;
    }

    value_type get_modulus() {
        return m_modulus;
    }

    // return (x * _R) % m_modulus
    value_type reduce(greater_value_type x) {
//        assert(x < upper_bound);
        greater_value_type m = ((x & (R - 1)) * (value_type)_M) & (R - 1);
        m = (greater_value_type)R - m; // 0 <= m < R
        greater_value_type X = m * m_modulus;
        if (x >= X) {
            X = (x - X) >> kModulusBitsize;
            while (X >= m_modulus) {
                X -= m_modulus;
            }
        } else {
            X = (X - x) >> kModulusBitsize;
            while (X > m_modulus) {
                X -= m_modulus;
            }
            X = m_modulus - X;
        }
//        assert(((X << kModulusBitsize) - x) % m_modulus == 0);
        return X;
    }

    value_type transform(value_type x) {
        return reduce((greater_value_type)x * square_mod);
    }

    value_type mulmod(value_type x, value_type y) {
        return reduce((greater_value_type)transform(x) * y);
    }

    value_type powmod(value_type a, value_type b) {
        value_type r = R;
        a = transform(a);
        while (b != 0) {
            if (b & 1) {
                r = reduce((greater_value_type)r * a);
            }
            a = reduce((greater_value_type)a * a);
            b >>= 1;
        }
        return reduce(r);
    }

    // simple mode
    value_type reduce_simple(greater_value_type x) {
        return x % m_modulus;
    }

    value_type mulmod_simple(value_type x, value_type y) {
        return reduce_simple((greater_value_type)x * y);
    }

    value_type powmod_simple(value_type a, value_type b) {
        value_type r = 1;
        while (b != 0) {
            if (b & 1) {
                r = mulmod_simple(r, a);
            }
            a = mulmod_simple(a, a);
            b >>= 1;
        }
        return reduce_simple(r);
    }

private:
    value_type m_modulus, R, square_mod;
    signed_value_type _M, _R; // R * _R - m_modulus * _M = 1
    greater_value_type upper_bound; // for input of `reduce` function
    unsigned int kModulusRepresentationBitsize, kModulusBitsize;
};

/*
 * Formal Reducer
 */

template<class T>
inline T mulmod(T x, T y, size_t p_index) {
    using greater_value_type = typename params<T>::greater_value_type;
    T modulus = params<T>::primes[p_index];
    T square_mod = params<T>::mont_square_mod[p_index];
    T _M = params<T>::mont_M[p_index];
    unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;

    greater_value_type z = (greater_value_type)x * square_mod;
    greater_value_type m = (T)((T)z * _M);
    T t = (z + m * modulus) >> kModulusRepresentationBitsize;
    z = (greater_value_type)t * y;
    m = (T)((T)z * _M);
    t = (z + m * modulus) >> kModulusRepresentationBitsize;
    return t;
}

template<class T>
inline T mulmod_barr(T x, T y, size_t p_index) {
    using greater_value_type = typename params<T>::greater_value_type;
    T modulus = params<T>::primes[p_index];
    T r = params<T>::barr_r[p_index];
    unsigned int kModulusBitsize = params<T>::kModulusBitsize;
    greater_value_type z = (greater_value_type)x * y;
    greater_value_type q = z >> (kModulusBitsize - 1);
    q *= r;
    q >>= (kModulusBitsize + 1);
    T t = z - q * modulus;
    while (t >= modulus) t -= modulus;
    return t;
}

template<class T>
inline T mulmod_simple(T x, T y, size_t p_index) {
    using greater_value_type = typename params<T>::greater_value_type;
    T modulus = params<T>::primes[p_index];
    return ((greater_value_type)x * y) % modulus;
}

}

}

#endif //TBLIB_OPS_HPP
