#ifndef TBLIB_KNUTH_YAO_HPP
#define TBLIB_KNUTH_YAO_HPP

//#define _DEBUG

#include "tbl/prng/sample.hpp"
#include "tbl/prng/randomInteger_mt19937.hpp"
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <immintrin.h>
#include <avx2intrin.h>

namespace tbl {

    static std::bernoulli_distribution uniformBin(0.5);

    template<size_t _sigma, size_t _lambda, size_t _tau>
    constexpr int sample<_sigma, _lambda, _tau>::lower_bound;

    template<size_t _sigma, size_t _lambda, size_t _tau>
    constexpr int sample<_sigma, _lambda, _tau>::upper_bound;

    template<size_t _sigma, size_t _lambda, size_t _tau>
    typename sample<_sigma, _lambda, _tau>::knuth_yao sample<_sigma, _lambda, _tau>::ky;

    template<size_t _sigma, size_t _lambda, size_t _tau>
    sample<_sigma, _lambda, _tau>::knuth_yao::knuth_yao() {
//        std::cout << "knuth_yao: construction()" << std::endl;
        initialize();
    }

    template<size_t _sigma, size_t _lambda, size_t _tau>
    sample<_sigma, _lambda, _tau>::sample() {
//        std::cout << "sample: construction()" << std::endl;
    }

    template<size_t _sigma, size_t _lambda, size_t _tau>
    void sample<_sigma, _lambda, _tau>::knuth_yao::initialize() {
//        std::cout << "knuth_yao: initialize()" << std::endl;
        NTL::RR::SetPrecision(lambda + 10); // 10 or other constants

        NTL::ZZ curZ = NTL::conv<NTL::ZZ>(lower_bound);
        NTL::ZZ upperZ = NTL::conv<NTL::ZZ>(upper_bound);

        NTL::RR subProbSum = NTL::RR(0.0);

        NTL::RR probR[upper_bound - lower_bound + 1];

        int i = 0;
        while (curZ <= upperZ) {
            NTL::RR tmp = NTL::RR(0.5) * NTL::power(NTL::conv<NTL::RR>(curZ) / NTL::conv<NTL::RR>((double)sigma), 2.0);
            NTL::RR subProb = NTL::exp(-tmp);
            subProbSum += subProb;
            probR[i] = subProb;
            i++;
            curZ += NTL::ZZ(1);
        }
        for (int i = 0; i < upper_bound - lower_bound + 1; i++) {
            probR[i] = probR[i] / subProbSum;
            // probR[i] = mm * 2^ee
            NTL::ZZ mm = probR[i].mantissa();
            long ee = probR[i].exponent();
            long lower_j = (-ee) > lambda ? -ee - lambda - 1 : -1;
            for (long j = -ee - 1; j > lower_j; j--) {
                if ((mm & (NTL::ZZ(1) << j)) != NTL::ZZ(0)) {
                    prob[i][-ee - 1 - j] = 1;
                } else {
                    prob[i][-ee - 1 - j] = 0;
                }
            }
            for (long j = -ee; j < lambda; j++) {
                prob[i][j] = 0;
            }
        }

        /*
         * AVX2 optimized pre-compute weight
         */
        for (int j = 0; j < lambda; j += 32) {
            __m256i x = _mm256_load_si256((__m256i*)(&prob[0][j]));
            for (int i = 1; i < upper_bound - lower_bound + 1; i++) {
                __m256i y = _mm256_load_si256((__m256i*)(&prob[i][j]));
                x = _mm256_add_epi8(x, y);
            }
            _mm256_store_si256((__m256i*)(&weight[j]), x);
        }

#ifdef _DEBUG
        /*
         * output prob matrix (debug)
         */
        std::cout << "prob matrix = [";
        for (int i = 0; i < upper_bound - lower_bound; i++) {
            std::cout << "0.";
            for (int j = 0; j < lambda - 1; j++) {
                std::cout << (int)prob[i][j];
            }
            std::cout << (int)prob[i][lambda - 1] << ", ";
        }
        std::cout << "0.";
        for (int j = 0; j < lambda -1; j++) {
            std::cout << (int)prob[upper_bound - lower_bound][j];
        }
        std::cout << (int)prob[upper_bound - lower_bound][lambda - 1] << "]" << std::endl;

        std::cout << "weight of each column = [";
        for (int i = 0; i < lambda - 1; i++) {
            std::cout << (int)weight[i] << ", ";
        }
        std::cout << (int)weight[lambda - 1] << "]" << std::endl;
#endif
    }

    template<size_t _sigma, size_t _lambda, size_t _tau>
    int sample<_sigma, _lambda, _tau>::knuth_yao::_rand() {
        long d = 0;
        int r = 0;
        // for possible `timing side-channel attack`
        for (int j = 0; j < lambda; j++) {
            d = 2 * d + uniformBin(ctx);
            if (d < weight[j]) {
                for (int i = upper_bound - lower_bound; i > -1; i--) {
                    d = d - prob[i][j];
                    if (d == -1) {
                        r = lower_bound + i;
                        for (int i2 = i - 1; i2 > -1; i2--) {
                            d = d - prob[i2][j];
                        }
                        break;
                    }
                }
                for (int j2 = j + 1; j2 < lambda; j2++) {
                    d = 2 * d + uniformBin(ctx);
                    d = d - weight[j2];
                }
                break;
            } else {
                d = d - weight[j];
            }
        }
        return r;
    }

}

#endif //TBLIB_KNUTH_YAO_HPP
