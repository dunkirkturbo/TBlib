#ifndef TBLIB_EUCLIDEAN_HPP
#define TBLIB_EUCLIDEAN_HPP

#include <tuple>

namespace tbl {

namespace helper {

template<class T>
std::tuple<T, T, T> exgcd(T n, T d, T u = 1, T v = 0, T u_ = 0, T v_ = 1) {
    /*
     * u  * N +  v * D = n
     * u_ * N + v_ * D = d
     * custom (u, v, u_, v_) for the case
     * N or D is unsigned variable, and the sign bit == 1, avoid transforming into `negative` signed variable
     */
    bool geq = true;    // n >= d
    if (n < d) {
        geq = false;
        T temp = n; n = d; d = temp;
        temp = u; u = v_; v_ = temp;
        temp = v; v = u_; u_ = temp;
    }
    while (d != 0) {
        T q = n / d;

        T temp = u_; u_ = u - q * u_; u = temp;
        temp = v_; v_ = v - q * v_; v = temp;
        temp = d; d = n - q * d; n = temp;
    }
    if (!geq) {
        T temp = u; u = v; v = temp;
    }
    return std::make_tuple(n, u, v);
}

template<class T>
T inverse(T n, T d) {
    T u = 1, v = 0, u_ = 0, v_ = 1, d_bak = d;
    bool geq = true;    // n >= d
    if (n < d) {
        geq = false;
        T temp = n; n = d; d = temp;
    }
    while (d != 0) {
        T q = n / d;

        T temp = u_; u_ = u - q * u_; u = temp;
        temp = v_; v_ = v - q * v_; v = temp;
        temp = d; d = n - q * d; n = temp;
    }
    if (!geq) {
        T temp = u; u = v; v = temp;
    }
    return u < 0 ? u + d_bak : u;
}

template<class T>
T gcd(T n, T d) {
    if (n < d) {
        T temp = n; n = d; d = temp;
    }
    while (d != 0) {
        T q = n / d;
        T temp = d; d = n - q * d; n = temp;
    }
    return n;
}

}

}
#endif //TBLIB_EUCLIDEAN_HPP
