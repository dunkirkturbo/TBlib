#ifndef TBLIB_RLWE_DEMO_HPP
#define TBLIB_RLWE_DEMO_HPP

#include "tbl/core.hpp"
#include "tbl/tools.hpp"

namespace tbl {

class rlwe_demo {
private:
    static constexpr size_t rlwe_degree = 128;
    static constexpr uint32_t rlwe_modulus = params<uint32_t>::primes[0];
    static constexpr uint32_t rlwe_half_modulus = rlwe_modulus / 2;
    sample<1, 192, 16> rlwe_prng;

public:
    using rlwe_poly = tbl::poly<uint32_t, rlwe_degree, 1>;

    rlwe_poly *rlwe_key, *rlwe_noise, *rlwe_half_q;
    rlwe_poly *rlwe_pt, *rlwe_ct;

    rlwe_demo() {
        // rlwe_key = {pub_g, pub_t, pub_s}
        rlwe_key = tbl::alloc_aligned<rlwe_poly, 32>(3);
        rlwe_pt = tbl::alloc_aligned<rlwe_poly, 32>(1);
        rlwe_ct = tbl::alloc_aligned<rlwe_poly, 32>(2);
        rlwe_noise = tbl::alloc_aligned<rlwe_poly, 32>(3);
        rlwe_half_q = tbl::alloc_aligned<rlwe_poly, 32>(1);
        uint32_t half_q_array[rlwe_degree];
        for (int i = 0; i < rlwe_degree; i++) {
            half_q_array[i] = rlwe_half_modulus;
        }
        rlwe_half_q->set(half_q_array, rlwe_degree);
    }

    void rlwe_keygen() {
        rlwe_poly *g = rlwe_key, *t = rlwe_key + 1, *s = rlwe_key + 2;
        rlwe_poly *e = rlwe_noise;
        s->get_noise(rlwe_prng);
        e->get_noise(rlwe_prng);
        uint32_t g_array[rlwe_degree];
        for (int i = 0; i < rlwe_degree; i++) {
            g_array[i] = tbl::getRandomRange(0, rlwe_modulus);
        }
        g->set(g_array, rlwe_degree);
        auto tmp = tbl::alloc_aligned<rlwe_poly, 32>(2);
        rlwe_poly *g2 = tmp, *s2 = tmp + 1;
        (*g2) = (*g);
        g2->ntt();
        (*s2) = (*s);
        s2->ntt();
        (*t) = (*g2) * (*s2);
        t->inv_ntt();
        (*t) = (*t) + (*e);
        free_aligned(2, tmp);
    }

    void rlwe_internal_encrypt(rlwe_poly &pt) {
        rlwe_poly *g = rlwe_key, *t = rlwe_key + 1;
        rlwe_poly *u = rlwe_ct, *v = rlwe_ct + 1;
        rlwe_poly *r = rlwe_noise, *e1 = rlwe_noise + 1, *e2 = rlwe_noise + 2;
        r->get_noise(rlwe_prng);
        e1->get_noise(rlwe_prng);
        e2->get_noise(rlwe_prng);

        auto tmp = tbl::alloc_aligned<rlwe_poly, 32>(3);
        rlwe_poly *g2 = tmp;
        (*g2) = (*g);
        g2->ntt();
        rlwe_poly *r2 = tmp + 1;
        (*r2) = (*r);
        r2->ntt();
        (*u) = (*g2) * (*r2);
        u->inv_ntt();
        (*u) = (*u) + (*e1);

        rlwe_poly *t2 = tmp + 2;
        (*t2) = (*t);
        t2->ntt();
        r->ntt();
        (*v) = (*t2) * (*r);
        v->inv_ntt();
        (*v) = (*v) + (*e2);
        (*v) = (*v) + ((*rlwe_half_q) * pt);    // without free, unsafe maybe

        free_aligned(3, tmp);
    }

    void rlwe_internal_decrypt(rlwe_poly &ct_u, rlwe_poly &ct_v) {
        rlwe_poly *s = rlwe_key + 2;
        auto tmp = tbl::alloc_aligned<rlwe_poly, 32>(2);
        rlwe_poly *u2 = tmp;
        (*u2) = ct_u;
        u2->ntt();
        rlwe_poly *s2 = tmp + 1;
        (*s2) = (*s);
        s2->ntt();
        (*rlwe_pt) = (*u2) * (*s2);
        rlwe_pt->inv_ntt();
        (*rlwe_pt) = ct_v - (*rlwe_pt);

        uint32_t *pt_array = rlwe_pt->coef();
        uint32_t lb = rlwe_modulus / 4, ub = 3 * lb;
        for (int i = 0; i < rlwe_degree; i++) {
            if (pt_array[i] > lb && pt_array[i] < ub) {
                pt_array[i] = 1;
            } else {
                pt_array[i] = 0;
            }
        }

//        std::cout << (*rlwe_pt) << std::endl;

        free_aligned(2, tmp);
    }

    void rlwe_encrypt(char* ct, const char* pt, size_t len) {
        rlwe_poly *u = rlwe_ct, *v = rlwe_ct + 1;
        int step_len = rlwe_degree * sizeof(uint32_t);
        uint32_t pt_array[rlwe_degree];
        for (int i = 0; i < len; i += 16) {
            int k = 128;
            if (i > len - 16) {
                k = (len - i) * 8;
            }
            for (int j = 0, i2 = i; j < k; j += 8, i2++) {
                pt_array[j] = (pt[i2] & 0x80) >> 7;
                pt_array[j + 1] = (pt[i2] & 0x40) >> 6;
                pt_array[j + 2] = (pt[i2] & 0x20) >> 5;
                pt_array[j + 3] = (pt[i2] & 0x10) >> 4;
                pt_array[j + 4] = (pt[i2] & 0x08) >> 3;
                pt_array[j + 5] = (pt[i2] & 0x04) >> 2;
                pt_array[j + 6] = (pt[i2] & 0x02) >> 1;
                pt_array[j + 7] = pt[i2] & 0x01;
            }
            rlwe_pt->set(pt_array, k);
            rlwe_internal_encrypt(*rlwe_pt);
            memcpy(ct, reinterpret_cast<char*>(u->coef()), step_len);
            ct += step_len;
            memcpy(ct, reinterpret_cast<char*>(v->coef()), step_len);
            ct += step_len;
        }
    }

    void rlwe_decrypt(char* pt, const char* ct, size_t len) {
        rlwe_poly *u = rlwe_ct, *v = rlwe_ct + 1;
        uint32_t *u_coef = u->coef(), *v_coef = v->coef(), *pt_array = rlwe_pt->coef();
        int step_len = rlwe_degree * sizeof(uint32_t);
        for (int i = 0; i < len; i += 16) {
            memcpy(reinterpret_cast<char*>(u_coef), ct, step_len);
            ct += step_len;
            memcpy(reinterpret_cast<char*>(v_coef), ct, step_len);
            ct += step_len;
            rlwe_internal_decrypt(*u, *v);
            int k = rlwe_degree;
            if (i > len - 16) {
                k = (len - i) * 8;
            }
            for (int j = 0, i2 = i; j < k; j += 8, i2++) {
                pt[i2] = (pt_array[j] << 7) | (pt_array[j + 1] << 6) | (pt_array[j + 2] << 5) | (pt_array[j + 3] << 4)
                        | (pt_array[j + 4] << 3) | (pt_array[j + 5] << 2) | (pt_array[j + 6] << 1) | pt_array[j + 7];
            }
        }
    }

    };

}

#endif //TBLIB_RLWE_DEMO_HPP
