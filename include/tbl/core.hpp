#ifndef TBLIB_CORE_HPP
#define TBLIB_CORE_HPP

#include "tbl/poly.hpp"
#include "tbl/helper/primer.hpp"
#include <algorithm>

namespace tbl {

    template<class T, size_t _degree, size_t _kModulusN>
    typename poly<T, _degree, _kModulusN>::core poly<T, _degree, _kModulusN>::base;

    template<class T, size_t _degree, size_t _kModulusN>
    poly<T, _degree, _kModulusN>::core::core() {
//        std::cout << "core: construction()" << std::endl;
        static_assert(kModulusN <= params<T>::kModulusMaxN, "core: not enough modulus for the poly (see params.hpp)");
        static_assert(degree <= params<T>::polyDegreeN, "core: degree is larger than polyDegreeN (see params.hpp)");
        degreeBinl = 0;
        size_t curDegree = 1;
        while (curDegree < degree) {
            degreeBinl++;
            curDegree <<= 1;
        }
        if (curDegree != degree) {
            throw std::runtime_error("core: degree is not in the form of 2^{...}");
        }
        /*
         * if there are no calls for `base` (e.g. there are no calls for poly::ntt())
         * class core will be regarded as dead_code by the compiler!!
         *
         * Thank GStalker for helping me solve this problem.
         */
        initialize();
    }

    /*
     * construction
     */
    template<class T, size_t _degree, size_t _kModulusN>
    poly<T, _degree, _kModulusN>::poly() {
        set({0});
    }

    template<class T, size_t _degree, size_t _kModulusN>
    poly<T, _degree, _kModulusN>::poly(std::initializer_list<value_type> values, bool reduce_coef) {
        set(values, reduce_coef);
    }

    template<class T, size_t _degree, size_t _kModulusN>
    poly<T, _degree, _kModulusN>::poly(iterator values, size_t len, bool reduce_coef) {
        set(values, len, reduce_coef);
    }

    template<class T, size_t _degree, size_t _kModulusN>
    void poly<T, _degree, _kModulusN>::set(std::initializer_list<value_type> values, bool reduce_coef) {
        set(values.begin(), values.end(), reduce_coef);
    }

    template<class T, size_t _degree, size_t _kModulusN>
    void poly<T, _degree, _kModulusN>::set(iterator values, size_t len, bool reduce_coef) {
        if (len > degree && len != N) {
            throw std::runtime_error("core: initializer of size above degree but not equal to N");
        }

        auto _end = values + len;
        for (int i = 0; i < kModulusN; i++) {
            value_type p = params<T>::primes[i];
            int j = 0;
            for (; j < degree && values < _end; j++, values++) {
                _coef[i * degree + j] = reduce_coef ? (*values) % p : (*values);
            }
            for (; j < degree; j++) {
                _coef[i * degree + j] = 0;
            }
        }
    }

    template<class T, size_t _degree, size_t _kModulusN>
    template<class It>
    void poly<T, _degree, _kModulusN>::set(It _begin, It _end, bool reduce_coef) {
        size_t size = std::distance(_begin, _end);
        if (size > degree && size != N) {
            throw std::runtime_error("core: initializer of size above degree but not equal to N");
        }

        auto iter = _begin;
        for (int i = 0; i < kModulusN; i++) {
            value_type p = params<T>::primes[i];
            int j = 0;
            for (; j < degree && iter < _end; j++, iter++) {
                _coef[i * degree + j] = reduce_coef ? (*iter) % p : (*iter);
            }
            for (; j < degree; j++) {
                _coef[i * degree + j] = 0;
            }
        }
    }

    template<class T, size_t _degree, size_t _kModulusN>
    poly<T, _degree, _kModulusN>::operator bool() const {
        return std::find_if(begin(), end(), [](value_type v) { return v != 0; }) != end();
    }

    /*
     * ostream
     */
    template<class T, size_t _degree, size_t _kModulusN>
    std::ostream& operator<<(std::ostream& outs, poly<T, _degree, _kModulusN> const& pl)
    {
        outs << "[";
        for (int i = 0; i < _kModulusN - 1; i++) {
            outs << "{";
            for (int j = 0; j < _degree - 1; j++) {
                outs << pl(i, j) << ", ";
            }
            outs << pl(i, _degree - 1) << "}, ";
        }
        outs << "{";
        for (int j = 0; j < _degree - 1; j++) {
            outs << pl(_kModulusN - 1, j) << ", ";
        }
        return outs << pl(_kModulusN - 1, _degree - 1) << "}]";
    }

    /*
     * pre_compute parameters for NTT
     */
    template<class T, size_t _degree, size_t _kModulusN>
    void poly<T, _degree, _kModulusN>::core::initialize() {
//        std::cout << "core: initialize()" << std::endl;
        unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        for (int currentModulus = 0; currentModulus < kModulusN; currentModulus++) {
            value_type p = params<T>::primes[currentModulus];
            value_type eta = helper::powmod(params<T>::primitive_roots[currentModulus], (p - 1) / (2 * degree), p);
            value_type omega = helper::powmod(eta, 2, p);
            value_type inv_eta = helper::inverse<signed_value_type>(eta, p);
            value_type inv_omega = helper::inverse<signed_value_type>(omega, p);
            if (degree) {
                inv_N[currentModulus] = helper::inverse<signed_value_type>(degree, p);
                inv_N[currentModulus] = ((greater_value_type)inv_N[currentModulus] << kModulusRepresentationBitsize) % p;
            }

            omegas[currentModulus][0] = ((greater_value_type)1 << kModulusRepresentationBitsize) % p;
            etas[currentModulus][0] = omegas[currentModulus][0];
            inv_omegas[currentModulus][0] = omegas[currentModulus][0];
            inv_etas[currentModulus][0] = omegas[currentModulus][0];
            for (int power = 1; power < degree; power++) {
                omegas[currentModulus][power] = ((greater_value_type)omegas[currentModulus][power - 1] * omega) % p;
                inv_omegas[currentModulus][power] = ((greater_value_type)inv_omegas[currentModulus][power - 1] * inv_omega) % p;
                etas[currentModulus][power] = ((greater_value_type)etas[currentModulus][power - 1] * eta) % p;
                inv_etas[currentModulus][power] = ((greater_value_type)inv_etas[currentModulus][power - 1] * inv_eta) % p;
            }

            int k = 0, curPos = degree - 2;
            for (int r = 1; r <= degreeBinl; r++) {
                for (int i = 0; i < (1 << (r - 1)); i++) {
                    int rev_i = reverse_bit(i, r - 1);
                    rev_i <<= (degreeBinl - r);
                    K[currentModulus][k] = omegas[currentModulus][rev_i];
                    inv_K[currentModulus][curPos - k] = inv_omegas[currentModulus][rev_i];
                    k++;
                }
            }
        }
    }

    template<class T, size_t _degree, size_t _kModulusN>
    bool poly<T, _degree, _kModulusN>::core::_ntt_avx2_16(value_type* x, const value_type* K_p, size_t p_index, size_t binLen) {
        // require degree >= 32 (although unchecked here)
        __m256i vec_modp, vec_mont_M;
        __m256i veca, vecb, veck, vecbk, vecbk_lo, vecbk_hi, vec_tmp, veca_t, vecb_t;

        __m256i z0 = _mm256_setzero_si256();
        __m256i z1 = _mm256_set1_epi16(1);
        vec_modp = _mm256_set1_epi16(params<T>::primes[p_index]);
        vec_mont_M = _mm256_set1_epi16(params<T>::mont_M[p_index]);
        value_type p2 = params<T>::primes2[p_index];

        int k = 0, NP = degree >> 1;
        for (int r = 1; r <= binLen - 4; r++) {
            int i = 0, j;
            for (; i < degree; i = j + NP) {
                veck = _mm256_set1_epi16(K_p[k]);
                for (j = i; j < i + NP; j += 16) {
                    veca = _mm256_load_si256((__m256i *)(&x[j]));
                    vecb = _mm256_load_si256((__m256i *)(&x[j + NP]));
                    // Montgomery for (vecb * veck) % vec_modp ==================
                    vecbk_lo = _mm256_mullo_epi16(vecb, veck);
                    vecbk_hi = _mm256_mulhi_epu16(vecb, veck);
                    vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
                    vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                    vecbk = _mm256_add_epi16(vec_tmp, vecbk_hi);
                    vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
                    vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                    vecbk = _mm256_add_epi16(vecbk, vec_tmp);
                    //===========================================================
                    vecb = _mm256_sub_epi16(veca, vecbk);
                    veca = _mm256_add_epi16(veca, vecbk);
                    _mm256_store_si256((__m256i*)(&x[j]), veca);
                    _mm256_store_si256((__m256i*)(&x[j + NP]), vecb);
                    for (int k1 = j, k2 = j + NP; k1 < j + 16; k1++, k2++) {
                        if (x[k1] >= p2) x[k1] -= p2;
                        if (x[k2] > p2) x[k2] += p2; // use overflow
                    }
                }
                k++;
            }
            NP >>= 1;
        }

        // 4th layer from the end (NP = 8)
        for (int i = 0; i < degree; i += 32) {
            veck = _mm256_set_epi16(K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1],
                                    K_p[k], K_p[k], K_p[k], K_p[k], K_p[k], K_p[k], K_p[k], K_p[k]);
            k += 2;
            veca = _mm256_load_si256((__m256i *)(&x[i]));
            vecb = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca = _mm256_permute4x64_epi64(veca, 0b11011000);
            vecb = _mm256_permute4x64_epi64(vecb, 0b11011000);
            veca_t = _mm256_unpacklo_epi64(veca, vecb);
            vecb_t = _mm256_unpackhi_epi64(veca, vecb);
            veca = _mm256_permute4x64_epi64(veca_t, 0b11011000);
            vecb = _mm256_permute4x64_epi64(vecb_t, 0b11011000);

            vecbk_lo = _mm256_mullo_epi16(vecb, veck);
            vecbk_hi = _mm256_mulhi_epu16(vecb, veck);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecbk = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecbk = _mm256_add_epi16(vecbk, vec_tmp);

            vecb = _mm256_sub_epi16(veca, vecbk);
            veca = _mm256_add_epi16(veca, vecbk);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j1 = i, j2 = i + 16; j1 < i + 16; j1++, j2++) {
                if (x[j1] >= p2) x[j1] -= p2;
                if (x[j2] > p2) x[j2] += p2; // use overflow
            }
        }
        // 3rd layer from the end (NP = 4)
        for (int i = 0; i < degree; i += 32) {
            veck = _mm256_set_epi16(K_p[k + 3], K_p[k + 3], K_p[k + 3], K_p[k + 3], K_p[k + 2], K_p[k + 2], K_p[k + 2], K_p[k + 2],
                                    K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k + 1], K_p[k], K_p[k], K_p[k], K_p[k]);
            k += 4;
            veca_t = _mm256_load_si256((__m256i *)(&x[i]));
            vecb_t = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca = _mm256_unpacklo_epi64(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi64(veca_t, vecb_t);

            vecbk_lo = _mm256_mullo_epi16(vecb, veck);
            vecbk_hi = _mm256_mulhi_epu16(vecb, veck);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecbk = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecbk = _mm256_add_epi16(vecbk, vec_tmp);

            vecb = _mm256_sub_epi16(veca, vecbk);
            veca = _mm256_add_epi16(veca, vecbk);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j1 = i, j2 = i + 16; j1 < i + 16; j1++, j2++) {
                if (x[j1] >= p2) x[j1] -= p2;
                if (x[j2] > p2) x[j2] += p2; // use overflow
            }
        }
        // 2nd layer from the end (NP = 2)
        for (int i = 0; i < degree; i += 32) {
            veck = _mm256_set_epi16(K_p[k + 7], K_p[k + 7], K_p[k + 6], K_p[k + 6], K_p[k + 5], K_p[k + 5], K_p[k + 4], K_p[k + 4],
                                    K_p[k + 3], K_p[k + 3], K_p[k + 2], K_p[k + 2], K_p[k + 1], K_p[k + 1], K_p[k], K_p[k]);
            k += 8;
            veca = _mm256_load_si256((__m256i *)(&x[i]));
            vecb = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca_t = _mm256_unpacklo_epi32(veca, vecb);
            vecb_t = _mm256_unpackhi_epi32(veca, vecb);
            veca = _mm256_unpacklo_epi64(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi64(veca_t, vecb_t);

            vecbk_lo = _mm256_mullo_epi16(vecb, veck);
            vecbk_hi = _mm256_mulhi_epu16(vecb, veck);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecbk = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecbk = _mm256_add_epi16(vecbk, vec_tmp);

            vecb = _mm256_sub_epi16(veca, vecbk);
            veca = _mm256_add_epi16(veca, vecbk);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j1 = i, j2 = i + 16; j1 < i + 16; j1++, j2++) {
                if (x[j1] >= p2) x[j1] -= p2;
                if (x[j2] > p2) x[j2] += p2; // use overflow
            }
        }
        // 1st layer from the end (NP = 1)
        for (int i = 0; i < degree; i += 32) {
            veck = _mm256_set_epi16(K_p[k + 15], K_p[k + 14], K_p[k + 13], K_p[k + 12], K_p[k + 11], K_p[k + 10], K_p[k + 9], K_p[k + 8],
                                    K_p[k + 7], K_p[k + 6], K_p[k + 5], K_p[k + 4], K_p[k + 3], K_p[k + 2], K_p[k + 1], K_p[k]);
            k += 16;
            veca_t = _mm256_load_si256((__m256i *)(&x[i]));
            vecb_t = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca = _mm256_unpacklo_epi16(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi16(veca_t, vecb_t);
            veca_t = _mm256_unpacklo_epi32(veca, vecb);
            vecb_t = _mm256_unpackhi_epi32(veca, vecb);
            veca = _mm256_unpacklo_epi32(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi32(veca_t, vecb_t);

            vecbk_lo = _mm256_mullo_epi16(vecb, veck);
            vecbk_hi = _mm256_mulhi_epu16(vecb, veck);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecbk = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecbk = _mm256_add_epi16(vecbk, vec_tmp);

            vecb = _mm256_sub_epi16(veca, vecbk);
            veca = _mm256_add_epi16(veca, vecbk);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j1 = i, j2 = i + 16; j1 < i + 16; j1++, j2++) {
                if (x[j1] >= p2) x[j1] -= p2;
                if (x[j2] > p2) x[j2] += p2; // use overflow
            }
        }
        return true;
    }

    template<class T, size_t _degree, size_t _kModulusN>
    bool poly<T, _degree, _kModulusN>::core::_ntt(value_type* x, const value_type* K_p, size_t p_index, size_t binLen) {
        value_type p = params<T>::primes[p_index];
        value_type mont_M = params<T>::mont_M[p_index];
        value_type p2 = params<T>::primes2[p_index];
        unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        int k = 0, NP = degree >> 1;
        for (int r = 1; r <= binLen; r++) {
            int i = 0, j;
            for (; i < degree; i = j + NP) {
                for (j = i; j < i + NP; j++) {
                    greater_value_type tmp = (greater_value_type)K_p[k] * x[j + NP];
                    greater_value_type m = (value_type)((value_type)tmp * mont_M);
                    value_type t = (tmp + m * p) >> kModulusRepresentationBitsize;
                    x[j + NP] = x[j] - t;
                    x[j] = x[j] + t;
                    if (x[j] >= p2) x[j] -= p2;
                    if (x[j + NP] > p2) x[j + NP] += p2;
                }
                k++;
            }
            NP >>= 1;
        }
        return true;
    }

    template<class T, size_t _degree, size_t _kModulusN>
    void poly<T, _degree, _kModulusN>::core::ntt(poly &pl) {
        unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        if (kModulusRepresentationBitsize == 16) {
            __m256i vec_modp, vec_mont_M, vec_etas, vec_coef, vec_lo, vec_hi, vec_tmp;
            __m256i z0 = _mm256_setzero_si256();
            __m256i z1 = _mm256_set1_epi16(1);
            for (int currentModulus = 0; currentModulus < kModulusN; currentModulus++) {
                value_type p = params<T>::primes[currentModulus];
                value_type mont_M = params<T>::mont_M[currentModulus];
                auto coef_cm = pl.coef() + currentModulus * degree;

                vec_modp = _mm256_set1_epi16(p);
                vec_mont_M = _mm256_set1_epi16(mont_M);
                for (int i = 0; i < degree; i += 16) {
                    vec_etas = _mm256_load_si256((__m256i *)(&etas[currentModulus][i]));
                    vec_coef = _mm256_load_si256((__m256i *)(&coef_cm[i]));
                    vec_lo = _mm256_mullo_epi16(vec_etas, vec_coef);
                    vec_hi = _mm256_mulhi_epu16(vec_etas, vec_coef);
                    vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);
                    vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                    vec_coef = _mm256_add_epi16(vec_tmp, vec_hi);
                    vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);
                    vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                    vec_coef = _mm256_add_epi16(vec_coef, vec_tmp);

                    _mm256_store_si256((__m256i *)(&coef_cm[i]), vec_coef);
//                    for (int j = i; j < i + 16; j++) {
//                        if (coef_cm[j] >= p) coef_cm[j] -= p;
//                    }
                }
                poly::core::_ntt_avx2_16(coef_cm, K[currentModulus], currentModulus, degreeBinl);
            }
        } else {
            for (int currentModulus = 0; currentModulus < kModulusN; currentModulus++) {
                value_type p = params<T>::primes[currentModulus];
                value_type mont_M = params<T>::mont_M[currentModulus];
                auto coef_cm = pl.coef() + currentModulus * degree;
                for (int i = 0; i < degree; i++) {
                    greater_value_type tmp = (greater_value_type)etas[currentModulus][i] * coef_cm[i];
                    greater_value_type m = (value_type)((value_type)tmp * mont_M);
                    coef_cm[i] = (tmp + m * p) >> kModulusRepresentationBitsize;
//                    value_type t = (tmp + m * p) >> kModulusRepresentationBitsize;
//                    coef_cm[i] = t >= p ? t - p : t;
                }
                poly::core::_ntt(coef_cm, K[currentModulus], currentModulus, degreeBinl);
                /*
                 * the below is for Barrett Reduce in `dot product`
                 */
                for (int i = 0; i < degree; i++) {
                    if (coef_cm[i] >= p) coef_cm[i] -= p;
                }
            }
        }
    }

    template<class T, size_t _degree, size_t _kModulusN>
    bool poly<T, _degree, _kModulusN>::core::_inv_ntt_avx2_16(value_type* x, const value_type* inv_K_p, size_t p_index, size_t binLen) {
        // require degree >= 32 (although unchecked here)
        __m256i vec_modp, vec_mont_M, vec_modp2;
        __m256i veca, vecb, vec_invk, vecbk_lo, vecbk_hi, vec_tmp, veca_t, vecb_t;

        __m256i z0 = _mm256_setzero_si256();
        __m256i z1 = _mm256_set1_epi16(1);
        vec_modp = _mm256_set1_epi16(params<T>::primes[p_index]);
        vec_mont_M = _mm256_set1_epi16(params<T>::mont_M[p_index]);
        value_type p2 = params<T>::primes2[p_index];
        vec_modp2 = _mm256_set1_epi16(p2);

        int k = 0;
        // 1st layer from the end (NP = 1)
        for (int i = degree - 32; i > -1; i -= 32) {
            vec_invk = _mm256_set_epi16(inv_K_p[k], inv_K_p[k + 1], inv_K_p[k + 2], inv_K_p[k + 3], inv_K_p[k + 4], inv_K_p[k + 5], inv_K_p[k + 6], inv_K_p[k + 7],
                                        inv_K_p[k + 8], inv_K_p[k + 9], inv_K_p[k + 10], inv_K_p[k + 11], inv_K_p[k + 12], inv_K_p[k + 13], inv_K_p[k + 14], inv_K_p[k + 15]);
            k += 16;
            veca_t = _mm256_load_si256((__m256i *)(&x[i]));
            vecb_t = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca = _mm256_add_epi16(veca_t, vecb_t);
            vecb = _mm256_sub_epi16(veca_t, vecb_t);
            vecb = _mm256_add_epi16(vecb, vec_modp2);

            vecbk_lo = _mm256_mullo_epi16(vecb, vec_invk);
            vecbk_hi = _mm256_mulhi_epu16(vecb, vec_invk);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecb = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecb = _mm256_add_epi16(vecb, vec_tmp);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j = i; j < i + 16; j++) {
                if (x[j] >= p2) x[j] -= p2;
            }
        }
        // 2nd layer from the end (NP = 2)
        for (int i = degree - 32; i > -1; i -= 32) {
            vec_invk = _mm256_set_epi16(inv_K_p[k], inv_K_p[k], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 2], inv_K_p[k + 2], inv_K_p[k + 3], inv_K_p[k + 3],
                                        inv_K_p[k + 4], inv_K_p[k + 4], inv_K_p[k + 5], inv_K_p[k + 5], inv_K_p[k + 6], inv_K_p[k + 6], inv_K_p[k + 7], inv_K_p[k + 7]);
            k += 8;

            veca = _mm256_load_si256((__m256i *)(&x[i]));
            vecb = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca_t = _mm256_unpacklo_epi16(veca, vecb);
            vecb_t = _mm256_unpackhi_epi16(veca, vecb);
            veca = _mm256_unpacklo_epi32(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi32(veca_t, vecb_t);
            veca_t = _mm256_unpacklo_epi32(veca, vecb);
            vecb_t = _mm256_unpackhi_epi32(veca, vecb);

            veca = _mm256_add_epi16(veca_t, vecb_t);
            vecb = _mm256_sub_epi16(veca_t, vecb_t);
            vecb = _mm256_add_epi16(vecb, vec_modp2);

            vecbk_lo = _mm256_mullo_epi16(vecb, vec_invk);
            vecbk_hi = _mm256_mulhi_epu16(vecb, vec_invk);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecb = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecb = _mm256_add_epi16(vecb, vec_tmp);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j = i; j < i + 16; j++) {
                if (x[j] >= p2) x[j] -= p2;
            }
        }
        // 3rd layer from the end (NP = 4)
        for (int i = degree - 32; i > -1; i -= 32) {
            vec_invk = _mm256_set_epi16(inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1],
                                        inv_K_p[k + 2], inv_K_p[k + 2], inv_K_p[k + 2], inv_K_p[k + 2], inv_K_p[k + 3], inv_K_p[k + 3], inv_K_p[k + 3], inv_K_p[k + 3]);
            k += 4;

            veca_t = _mm256_load_si256((__m256i *)(&x[i]));
            vecb_t = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca = _mm256_unpacklo_epi32(veca_t, vecb_t);
            vecb = _mm256_unpackhi_epi32(veca_t, vecb_t);
            veca_t = _mm256_unpacklo_epi64(veca, vecb);
            vecb_t = _mm256_unpackhi_epi64(veca, vecb);

            veca = _mm256_add_epi16(veca_t, vecb_t);
            vecb = _mm256_sub_epi16(veca_t, vecb_t);
            vecb = _mm256_add_epi16(vecb, vec_modp2);

            vecbk_lo = _mm256_mullo_epi16(vecb, vec_invk);
            vecbk_hi = _mm256_mulhi_epu16(vecb, vec_invk);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecb = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecb = _mm256_add_epi16(vecb, vec_tmp);
            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j = i; j < i + 16; j++) {
                if (x[j] >= p2) x[j] -= p2;
            }
        }
        // 4th layer from the end (NP = 8)
        for (int i = degree - 32; i > -1; i -= 32) {
            vec_invk = _mm256_set_epi16(inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k], inv_K_p[k],
                                        inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1], inv_K_p[k + 1]);
            k += 2;

            veca = _mm256_load_si256((__m256i *)(&x[i]));
            vecb = _mm256_load_si256((__m256i *)(&x[i + 16]));
            veca_t = _mm256_unpacklo_epi64(veca, vecb);
            vecb_t = _mm256_unpackhi_epi64(veca, vecb);

            veca = _mm256_add_epi16(veca_t, vecb_t);
            vecb = _mm256_sub_epi16(veca_t, vecb_t);
            vecb = _mm256_add_epi16(vecb, vec_modp2);

            vecbk_lo = _mm256_mullo_epi16(vecb, vec_invk);
            vecbk_hi = _mm256_mulhi_epu16(vecb, vec_invk);
            vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
            vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
            vecb = _mm256_add_epi16(vec_tmp, vecbk_hi);
            vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
            vec_tmp = _mm256_add_epi16(vec_tmp, z1);
            vecb = _mm256_add_epi16(vecb, vec_tmp);

            veca = _mm256_permute4x64_epi64(veca, 0b11011000);
            vecb = _mm256_permute4x64_epi64(vecb, 0b11011000);
            veca_t = _mm256_unpacklo_epi64(veca, vecb);
            vecb_t = _mm256_unpackhi_epi64(veca, vecb);
            veca = _mm256_permute4x64_epi64(veca_t, 0b11011000);
            vecb = _mm256_permute4x64_epi64(vecb_t, 0b11011000);

            _mm256_store_si256((__m256i*)(&x[i]), veca);
            _mm256_store_si256((__m256i*)(&x[i + 16]), vecb);
            for (int j = i; j < i + 8; j++) {
                if (x[j] >= p2) x[j] -= p2;
            }
            for (int j = i + 16; j < i + 24; j++) {
                if (x[j] >= p2) x[j] -= p2;
            }
        }

        int NP = 16;
        for (int r = 5; r <= binLen; r++) {
            for (int i = degree - 2 * NP; i > -1; i -= 2 * NP) {
                vec_invk = _mm256_set1_epi16(inv_K_p[k]);
                for (int j = i; j < i + NP; j += 16) {
                    veca_t = _mm256_load_si256((__m256i *)(&x[j]));
                    vecb_t = _mm256_load_si256((__m256i *)(&x[j + NP]));

                    veca = _mm256_add_epi16(veca_t, vecb_t);
                    vecb = _mm256_sub_epi16(veca_t, vecb_t);
                    vecb = _mm256_add_epi16(vecb, vec_modp2);

                    vecbk_lo = _mm256_mullo_epi16(vecb, vec_invk);
                    vecbk_hi = _mm256_mulhi_epu16(vecb, vec_invk);
                    vec_tmp = _mm256_mullo_epi16(vecbk_lo, vec_mont_M);
                    vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                    vecb = _mm256_add_epi16(vec_tmp, vecbk_hi);
                    vec_tmp = _mm256_cmpeq_epi16(vecbk_lo, z0);
                    vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                    vecb = _mm256_add_epi16(vecb, vec_tmp);
                    _mm256_store_si256((__m256i*)(&x[j]), veca);
                    _mm256_store_si256((__m256i*)(&x[j + NP]), vecb);
                    for (int k1 = j; k1 < j + 16; k1++) {
                        if (x[k1] >= p2) x[k1] -= p2;
                    }
                }
                k++;
            }
            NP <<= 1;
        }
        return true;
    }

    /*
     * inv_ntt.GS operations' processing order can be completely the reverse of ntt.CT operations
     */
    template<class T, size_t _degree, size_t _kModulusN>
    bool poly<T, _degree, _kModulusN>::core::_inv_ntt(value_type* x, const value_type* inv_K_p, size_t p_index, size_t binLen) {
        value_type p = params<T>::primes[p_index];
        value_type mont_M = params<T>::mont_M[p_index];
        value_type p2 = params<T>::primes2[p_index];
        unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;

        int k = 0, NP = 1;
        for (int r = 1; r <= binLen; r++) {
            int i = degree - 1, j;
            for (; i > -1; i = j - NP) {
                for (j = i; j > i - NP; j--) {
                    value_type tmp = x[j - NP];
                    x[j - NP] = x[j] + tmp;
                    if (x[j - NP] >= p2) x[j - NP] -= p2;
                    x[j] = tmp - x[j];
                    if (x[j] > p2) x[j] += p2;

                    greater_value_type tmp1 = (greater_value_type)inv_K_p[k] * x[j];
                    greater_value_type m = (value_type)((value_type)tmp1 * mont_M);
                    x[j] = (tmp1 + m * p) >> kModulusRepresentationBitsize;
                }
                k++;
            }
            NP <<= 1;
        }
        return true;
    }

    template<class T, size_t _degree, size_t _kModulusN>
    void poly<T, _degree, _kModulusN>::core::inv_ntt(poly &pl) {
        unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
        if (kModulusRepresentationBitsize == 16) {
            __m256i vec_modp, vec_mont_M, vec_inv_etas, vec_inv_N, vec_coef, vec_lo, vec_hi, vec_tmp;
            __m256i z0 = _mm256_setzero_si256();
            __m256i z1 = _mm256_set1_epi16(1);
            for (int currentModulus = 0; currentModulus < kModulusN; currentModulus++) {
                value_type p = params<T>::primes[currentModulus];
                value_type mont_M = params<T>::mont_M[currentModulus];
                auto coef_cm = pl.coef() + currentModulus * degree;

                vec_modp = _mm256_set1_epi16(p);
                vec_mont_M = _mm256_set1_epi16(mont_M);
                vec_inv_N = _mm256_set1_epi16(inv_N[currentModulus]);
                poly::core::_inv_ntt_avx2_16(coef_cm, inv_K[currentModulus], currentModulus, degreeBinl);
                for (int i = 0; i < degree; i += 16) {
                    vec_inv_etas = _mm256_load_si256((__m256i *)(&inv_etas[currentModulus][i]));
                    vec_coef = _mm256_load_si256((__m256i *)(&coef_cm[i]));
                    vec_lo = _mm256_mullo_epi16(vec_inv_etas, vec_coef);
                    vec_hi = _mm256_mulhi_epu16(vec_inv_etas, vec_coef);
                    vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);
                    vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                    vec_coef = _mm256_add_epi16(vec_tmp, vec_hi);
                    vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);
                    vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                    vec_coef = _mm256_add_epi16(vec_coef, vec_tmp);

                    vec_lo = _mm256_mullo_epi16(vec_inv_N, vec_coef);
                    vec_hi = _mm256_mulhi_epu16(vec_inv_N, vec_coef);
                    vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);
                    vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                    vec_coef = _mm256_add_epi16(vec_tmp, vec_hi);
                    vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);
                    vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                    vec_coef = _mm256_add_epi16(vec_coef, vec_tmp);

                    _mm256_store_si256((__m256i *)(&coef_cm[i]), vec_coef);
                    for (int j = i; j < i + 16; j++) {
                        if (coef_cm[j] >= p) coef_cm[j] -= p;
                    }
                }
            }
        } else {
            for (int currentModulus = 0; currentModulus < kModulusN; currentModulus++) {
                value_type p = params<T>::primes[currentModulus];
                value_type mont_M = params<T>::mont_M[currentModulus];
                auto coef_cm = pl.coef() + currentModulus * degree;
                poly::core::_inv_ntt(coef_cm, inv_K[currentModulus], currentModulus, degreeBinl);
                for (int i = 0; i < degree; i++) {
                    greater_value_type tmp = (greater_value_type)inv_etas[currentModulus][i] * coef_cm[i];
                    greater_value_type m = (value_type)((value_type)tmp * mont_M);
                    value_type t = (tmp + m * p) >> kModulusRepresentationBitsize;
                    tmp = (greater_value_type)inv_N[currentModulus] * t;
                    m = (value_type)((value_type)tmp * mont_M);
                    t = (tmp + m * p) >> kModulusRepresentationBitsize;
                    coef_cm[i] = t >= p ? t - p : t;
                }
            }
        }
    }

    template<class T, size_t _degree, size_t _kModulusN>
    template<size_t _sigma, size_t _lambda, size_t _tau>
    void poly<T, _degree, _kModulusN>::get_noise(sample<_sigma, _lambda, _tau> &g_prng) {
        for (int i = 0; i < kModulusN; i++) {
            value_type p = params<T>::primes[i];
            int j = 0;
            for (int j = 0; j < degree; j++) {
                signed_value_type noise = g_prng.rand();
                if (noise < 0) {
                    _coef[i * degree + j] = (signed_value_type) p + noise;
                } else {
                    _coef[i * degree + j] = noise;
                }
            }
        }
    }

}

#endif //TBLIB_CORE_HPP