#ifndef TBLIB_POLY_HPP
#define TBLIB_POLY_HPP

#include "tbl/params.hpp"
#include "tbl/tools.hpp"
#include "tbl/prng/knuth_yao.hpp"
#include <cstddef>
#include <iostream>
#include <initializer_list>
#include <immintrin.h>
#include <avx2intrin.h>

namespace tbl {

    template<class T, size_t _degree, size_t _kModulusN>
    class poly {

    public:
        using value_type = typename params<T>::value_type;
        using greater_value_type = typename params<T>::greater_value_type;
        using signed_value_type = typename params<T>::signed_value_type;

        using iterator = T*;
        using const_iterator = T const*;

        static constexpr size_t degree = _degree;
        static constexpr size_t kModulusN = _kModulusN;
        static constexpr size_t kbits = params<T>::kModulusBitsize;
        static constexpr size_t aggregated_modulus_bits = kbits * kModulusN;    // for modulus used actually

        /*
         * construction
         */
        poly();
        poly(std::initializer_list<value_type> values, bool reduce_coef = true);
        poly(iterator values, size_t len, bool reduce_coef = true);

        void set(std::initializer_list<value_type> values, bool reduce_coef = true);
        void set(iterator values, size_t len, bool reduce_coef = true);
        template<class It> void set(It _begin, It _end, bool reduce_coef = true);

        /*
         * public API
         */
        static constexpr value_type get_modulus(size_t k_index) { return params<T>::primes[k_index]; }
        void ntt() { base.ntt(*this); }
        void inv_ntt() { base.inv_ntt(*this); }

        template<size_t _sigma, size_t _lambda, size_t _tau>
        void get_noise(sample<_sigma, _lambda, _tau> &g_prng);

        /*
         * explicit conversion (poly => bool)
         */
        explicit operator bool() const;

        /*
         * iterator
         */
        iterator begin() { return std::begin(_coef); }
        iterator end() { return std::end(_coef); }
        const_iterator const_begin() { return std::begin(_coef); }
        const_iterator const_end() { return std::end(_coef); }
        value_type* coef() { return _coef; }

        /*
         * operators
         */
        poly& operator=(std::initializer_list<value_type> values) { set(values); return *this; }
        value_type const& operator()(size_t k_index, size_t i) const { return _coef[k_index * degree + i]; }
        poly& operator+(poly &pl) {
            auto new_pl = alloc_aligned<poly<T, _degree, _kModulusN>, 32>(1);
            value_type* new_coef = new_pl->coef();
            value_type* pl_coef = pl.coef();
            for (int i = 0; i < kModulusN; i++) {
                value_type p = params<T>::primes[i];
                for (int j = i * degree; j < (i + 1) * degree; j++) {
                    new_coef[j] = ((greater_value_type)this->_coef[j] + pl_coef[j]) % p;
                }
            }
            return *new_pl;
        }
        poly& operator-(poly &pl) {
            auto new_pl = alloc_aligned<poly<T, _degree, _kModulusN>, 32>(1);
            value_type* new_coef = new_pl->coef();
            value_type* pl_coef = pl.coef();
            for (int i = 0; i < kModulusN; i++) {
                value_type p = params<T>::primes[i];
                for (int j = i * degree; j < (i + 1) * degree; j++) {
                    new_coef[j] = this->_coef[j] >= pl_coef[j] ?
                                  this->_coef[j] - pl_coef[j] : p - pl_coef[j] + this->_coef[j];
                }
            }
            return *new_pl;
        }
        // component-wise multiplication
        poly& operator*(poly &pl) {
            auto new_pl = alloc_aligned<poly<T, _degree, _kModulusN>, 32>(1);
            value_type* new_coef = new_pl->coef();
            value_type* pl_coef = pl.coef();
            unsigned int kModulusRepresentationBitsize = params<T>::kModulusRepresentationBitsize;
            if (kModulusRepresentationBitsize == 16) {
                __m256i vec_modp, vec_mont_M, vec_mont_sm, veca, vecb, vecc, vec_lo, vec_hi, vec_tmp;
                __m256i z0 = _mm256_setzero_si256();
                __m256i z1 = _mm256_set1_epi16(1);
                for (int i = 0; i < kModulusN; i++) {
                    vec_modp = _mm256_set1_epi16(params<T>::primes[i]);
                    vec_mont_M = _mm256_set1_epi16(params<T>::mont_M[i]);
                    vec_mont_sm = _mm256_set1_epi16(params<T>::mont_square_mod[i]);
                    for (int j = i * degree; j < (i + 1) * degree; j += 16) {
                        veca = _mm256_load_si256((__m256i *)(&(this->_coef[j])));
                        vecb = _mm256_load_si256((__m256i *)(pl_coef + j));
                        vec_lo = _mm256_mullo_epi16(veca, vec_mont_sm);
                        vec_hi = _mm256_mulhi_epu16(veca, vec_mont_sm);
                        vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);
                        vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                        veca = _mm256_add_epi16(vec_tmp, vec_hi);
                        vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);
                        vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                        veca = _mm256_add_epi16(veca, vec_tmp);

                        vec_lo = _mm256_mullo_epi16(veca, vecb);
                        vec_hi = _mm256_mulhi_epu16(veca, vecb);
                        vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);
                        vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);
                        vecc = _mm256_add_epi16(vec_tmp, vec_hi);
                        vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);
                        vec_tmp = _mm256_add_epi16(vec_tmp, z1);
                        vecc = _mm256_add_epi16(vecc, vec_tmp);
                        _mm256_store_si256((__m256i *)(new_coef + j), vecc);
                    }
                }
            } else {
                for (int i = 0; i < kModulusN; i++) {
                    // Barrett Reduce instead
                    value_type p = params<T>::primes[i];
                    value_type r = params<T>::barr_r[i];
                    unsigned int kModulusBitsize = params<T>::kModulusBitsize;
                    for (int j = i * degree; j < (i + 1) * degree; j++) {
                        greater_value_type z = (greater_value_type)this->_coef[j] * pl_coef[j];
                        greater_value_type q = z >> (kModulusBitsize - 1);
                        q *= r;
                        q >>= (kModulusBitsize + 1);
                        new_coef[j] = z - q * p;
                        while (new_coef[j] >= p) new_coef[j] -= p;
                    }
                }
            }
            return *new_pl;
        }

        /*
         * serializer
         */
        void serialize(std::ostream& outputstream) {
            outputstream.write(reinterpret_cast<char*>(_coef), N * sizeof(T));
        }
        void unserialize(std::istream& inputstream) {
            inputstream.read(reinterpret_cast<char*>(_coef), N * sizeof(T));
        }

    protected:
        class core {
        public:
            core();
            static bool _ntt(value_type* x, const value_type* K_p, size_t p_index, size_t binLen);
            static bool _inv_ntt(value_type* x, const value_type* inv_K_p, size_t p_index, size_t binLen);
            static bool _ntt_avx2_16(value_type* x, const value_type* K_p, size_t p_index, size_t binLen);
            static bool _inv_ntt_avx2_16(value_type* x, const value_type* inv_K_p, size_t p_index, size_t binLen);
            void ntt(poly& pl);
            void inv_ntt(poly& pl);
            void initialize();

        private:
            size_t degreeBinl;
            value_type omegas[kModulusN][degree] __attribute__((aligned(32)));
            value_type etas[kModulusN][degree] __attribute__((aligned(32)));
            value_type inv_omegas[kModulusN][degree] __attribute__((aligned(32)));
            value_type inv_etas[kModulusN][degree] __attribute__((aligned(32)));

            value_type K[kModulusN][degree] __attribute__((aligned(32)));   // factor for the CT operation(蝶件)
            value_type inv_K[kModulusN][degree] __attribute__((aligned(32)));

            value_type inv_N[kModulusN] __attribute__((aligned(32)));

            // reverse the lower bitLen-bit of integer i
            static int reverse_bit(int i, int bitLen) {
                if (bitLen == 0) return 0;
                int curPower = 1 << (bitLen - 1), rev_i = 0;
                while (i) {
                    if (i & 1) {
                        rev_i += curPower;
                    }
                    i >>= 1;
                    curPower >>= 1;
                }
                return rev_i;
            }

        } __attribute__((aligned(64)));

        static core base;

    private:
        static constexpr size_t N = _degree * _kModulusN;
        T _coef[N] __attribute__((aligned(32))); // aligned(32) for AVX2
    };

}

#endif //TBLIB_POLY_HPP
