#include <iostream>
#include "tbl/params.hpp"
#include "tbl/ops.hpp"
#include "tbl/prng/randomInteger_mt19937.hpp"
#include "tbl/helper/primer.hpp"
#include "tbl/helper/primer_select.hpp"
#include "tbl/helper/primitive_root_select.hpp"
#include "utils/bench.h"
#include "tbl/core.hpp"
#include "tbl/layer3/rlwe_demo.hpp"
#include <NTL/ZZ_pX.h>

#define TEST_MULMOD
#define TEST_MULMOD_AVX2
#define TEST_POLYMUL_16
//#define TEST_POLYMUL_32
#define TEST_RLWE_CRYPT

template<class T>
void run() {
    using value_type = typename tbl::params<T>::value_type;
    using greater_value_type = typename tbl::params<T>::greater_value_type;

    value_type modulus, x, y;
    int benches = 100, rounds = 1000;

#ifdef TEST_MULMOD
    modulus = tbl::params<T>::primes[0];
    for (int i = 0; i < 0xDEAD; i++) {
        x = tbl::getRandomRange(0, modulus);
        y = tbl::getRandomRange(0, modulus);
        assert(tbl::ops::mulmod_barr(x, y, 0) == tbl::ops::mulmod_simple(x, y, 0));
    }
    std::cout << "[PASSED] correctness of reducer_Barrett.mulmod" << std::endl;
    BENCH_START("Time per Barr_common.mulmod", benches);
    BENCH_ITEM(tbl::ops::mulmod_barr(x, y, 0), rounds);
    BENCH_FINAL(benches, rounds);

    for (int i = 0; i < 0xDEAD; i++) {
        x = tbl::getRandomRange(0, modulus);
        y = tbl::getRandomRange(0, modulus);
        assert(tbl::ops::mulmod(x, y, 0) % modulus == tbl::ops::mulmod_simple(x, y, 0));
    }
    std::cout << "[PASSED] correctness of reducer_Montgomery.mulmod" << std::endl;
    BENCH_START("Time per Mont_common.mulmod", benches);
    BENCH_ITEM(tbl::ops::mulmod(x, y, 0), rounds);
    BENCH_FINAL(benches, rounds);

    BENCH_START("Time per simple.mulmod", benches);
    BENCH_ITEM(tbl::ops::mulmod_simple(x, y, 0), rounds);
    BENCH_FINAL(benches, rounds);
#endif
}

void pre_compute_primes();
void pre_compute_primitive_root();

int main() {
//    pre_compute_primes();
//    pre_compute_primitive_root();
    run<uint64_t>();
    std::cout << std::endl;

    int benches, rounds;

#ifdef TEST_MULMOD_AVX2
    uint16_t a[16] __attribute__((aligned(32))) = {2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333, 2333};
    uint16_t b = 7777;
    __m256i veca, vecb, vecab, vec_tmp, vec_lo, vec_hi;
    __m256i vec_mont_M = _mm256_set1_epi16(tbl::params<uint16_t>::mont_M[0]);
    __m256i vec_modp = _mm256_set1_epi16(tbl::params<uint16_t>::primes[0]);
    __m256i z0 = _mm256_setzero_si256();
    __m256i z1 = _mm256_set1_epi16(1);

    benches = 100;
    rounds = 1000;
    BENCH_START("Time per avx2.mulmod (including 16 * uint16_t)", benches);
    BENCH_ITEM(veca = _mm256_load_si256((__m256i *)a);vecb = _mm256_set1_epi16(b);vec_lo = _mm256_mullo_epi16(veca, vecb);vec_hi = _mm256_mulhi_epu16(veca, vecb);vec_tmp = _mm256_mullo_epi16(vec_lo, vec_mont_M);vec_tmp = _mm256_mulhi_epu16(vec_tmp, vec_modp);vecab = _mm256_add_epi16(vec_tmp, vec_hi);vec_tmp = _mm256_cmpeq_epi16(vec_lo, z0);vec_tmp = _mm256_add_epi16(vec_tmp, z1);vecab = _mm256_add_epi16(vecab, vec_tmp);_mm256_store_si256((__m256i *)a, vecab);, rounds);
    BENCH_FINAL(benches, rounds);
    std::cout << std::endl;
#endif

    /*
     * AVX2 16-bit
     */
#ifdef TEST_POLYMUL_16
    NTL::ZZ p(tbl::params<uint16_t>::primes[0]);
    NTL::ZZ_p::init(p);
    const int d = 32;
    NTL::ZZ_pX f;
    NTL::SetCoeff(f, d);
    NTL::SetCoeff(f, 0, 1);
    NTL::ZZ_pXModulus F(f);
    NTL::ZZ_pX g, h;
    NTL::random(g, d - 1);
    NTL::random(h, d - 1);
    NTL::ZZ_pX res;

    benches = 50;
    rounds = 200;
    BENCH_START("Time per NTL.MulMod", benches);
    BENCH_ITEM(NTL::MulMod(res, g, h, F), rounds);
//    BENCH_ITEM(NTL::mul(res, g, h);, rounds); // mul use FFT(not mod), but slow also
    BENCH_FINAL(benches, rounds);

    auto pl1 = tbl::alloc_aligned<tbl::poly<uint16_t, d, 1>, 32>(1);
    auto pl2 = tbl::alloc_aligned<tbl::poly<uint16_t, d, 1>, 32>(1);
    uint16_t coef[2 * d];
    for (int i = 0; i < 2 * d; i++) {
        coef[i] = tbl::getRandomNBitInteger(15);
    }
    pl1->set(coef, d);
    pl2->set(coef + d, d);

    BENCH_START("Time per TBL.MulMod", benches);
    BENCH_ITEM(pl1->ntt();pl2->ntt();auto pl = (*pl1) * (*pl2);pl.inv_ntt();, rounds);
    BENCH_FINAL(benches, rounds);
#endif

    /*
     * 32-bit
     */
#ifdef TEST_POLYMUL_32
    NTL::ZZ p(tbl::params<uint32_t>::primes[0]);
    NTL::ZZ_p::init(p);
    const int d = 32;
    NTL::ZZ_pX f;
    NTL::SetCoeff(f, d);
    NTL::SetCoeff(f, 0, 1);
    NTL::ZZ_pX g, h;
    NTL::random(g, d - 1);
    NTL::random(h, d - 1);
    NTL::ZZ_pX res;

    benches = 50;
    rounds = 200;
    BENCH_START("Time per NTL.MulMod", benches);
    BENCH_ITEM(NTL::MulMod(res, g, h, f);, rounds);
    BENCH_FINAL(benches, rounds);

    auto pl1 = tbl::alloc_aligned<tbl::poly<uint32_t, d, 1>, 32>(1);
    auto pl2 = tbl::alloc_aligned<tbl::poly<uint32_t, d, 1>, 32>(1);
    uint32_t coef[2 * d];
    for (int i = 0; i < 2 * d; i++) {
        coef[i] = tbl::getRandomNBitInteger(31);
    }
    pl1->set(coef, d);
    pl2->set(coef + d, d);
    BENCH_START("Time per TBL.MulMod", benches);
    BENCH_ITEM(pl1->ntt();pl2->ntt();auto pl = (*pl1) * (*pl2);pl.inv_ntt();, rounds);
    BENCH_FINAL(benches, rounds);
#endif

#ifdef TEST_RLWE_CRYPT
    tbl::rlwe_demo rl;
    rl.rlwe_keygen();

    char* pt = (char*)"0xDktb@rlwe_demo";
    char ct[1024];
    char pt_edit[16];

    rl.rlwe_encrypt(ct, pt, 16);
    uint8_t *ct_out = (uint8_t*)ct;
    for (int i = 0; i < 1024; i++) {
        printf("%02x", ct_out[i]);
    }
    printf("\n");
    rl.rlwe_decrypt(pt_edit, ct, 16);
    for (int i = 0; i < 16; i++) {
        printf("%c", pt_edit[i]);
    }
    printf("\n");
#endif

    return 0;
}

void pre_compute_primes() {
    std::cout << "=========================================================" << std::endl;
    std::cout << "|   Generation of modulus<uint16_t/uint32_t/uint64_t>   |" << std::endl;
    std::cout << "=========================================================" << std::endl;
    uint16_t primes_16[tbl::helper::DEFAULT_MAX_PN];
    uint32_t primes_32[tbl::helper::DEFAULT_MAX_PN];
    uint64_t primes_64[tbl::helper::DEFAULT_MAX_PN];
    int PN_16 = tbl::helper::genPrime<uint16_t>(primes_16);
    int PN_32 = tbl::helper::genPrime<uint32_t>(primes_32);
    int PN_64 = tbl::helper::genPrime<uint64_t>(primes_64);
    if (PN_16) {
        std::cout << "primers_16[" << PN_16 << "] = {";
        for (int i = 0; i < PN_16 - 1; i++) {
            std::cout << primes_16[i] << ", ";
        }
        std::cout << primes_16[PN_16 - 1] << "}" << std::endl;
    }
    if (PN_32) {
        std::cout << "primers_32[" << PN_32 << "] = {";
        for (int i = 0; i < PN_32 - 1; i++) {
            std::cout << primes_32[i] << ", ";
        }
        std::cout << primes_32[PN_32 - 1] << "}" << std::endl;
    }
    if (PN_64) {
        std::cout << "primers_64[" << PN_64 << "] = {";
        for (int i = 0; i < PN_64 - 1; i++) {
            std::cout << primes_64[i] << ", ";
        }
        std::cout << primes_64[PN_64 - 1] << "}" << std::endl;
    }
}

void pre_compute_primitive_root() {
    std::cout << "=================================================================" << std::endl;
    std::cout << "|   Generation of primitive roots<uint16_t/uint32_t/uint64_t>   |" << std::endl;
    std::cout << "=================================================================" << std::endl;
    auto PN_16 = tbl::params<uint16_t>::kModulusMaxN;
    auto PN_32 = tbl::params<uint32_t>::kModulusMaxN;
    auto PN_64 = tbl::params<uint64_t>::kModulusMaxN;

    uint16_t prim_16[PN_16];
    tbl::helper::genPrimitiveRoot<uint16_t>(prim_16);
    std::cout << "prim_16[" << PN_16 << "] = {";
    for (int i = 0; i < PN_16 - 1; i++) {
        std::cout << prim_16[i] << ", ";
    }
    std::cout << prim_16[PN_16 - 1] << "}" << std::endl;

    uint32_t prim_32[PN_32];
    tbl::helper::genPrimitiveRoot<uint32_t>(prim_32);
    std::cout << "prim_32[" << PN_32 << "] = {";
    for (int i = 0; i < PN_32 - 1; i++) {
        std::cout << prim_32[i] << ", ";
    }
    std::cout << prim_32[PN_32 - 1] << "}" << std::endl;

    uint64_t prim_64[PN_64];
    tbl::helper::genPrimitiveRoot<uint64_t>(prim_64);
    std::cout << "prim_64[" << PN_64 << "] = {";
    for (int i = 0; i < PN_64 - 1; i++) {
        std::cout << prim_64[i] << ", ";
    }
    std::cout << prim_64[PN_64 - 1] << "}" << std::endl;
}
