#ifndef TBLIB_SAMPLE_HPP
#define TBLIB_SAMPLE_HPP

namespace tbl {

/*
 * Discrete Gaussian Distribution, x
 * -(_tau * _sigma) <= x <= (_tau * _sigma)
 * Statistical distance to an ideal discrete Gaussian distribution <= (1/2) ^ {_lambda}
 * _tau ≈ sqrt{_lambda * 2ln2}
 */
template<size_t _sigma, size_t _lambda, size_t _tau>
class sample {

public:
    static constexpr size_t sigma = _sigma;
    static constexpr size_t lambda = _lambda;
    static constexpr size_t tau = _tau;
    static constexpr int lower_bound = (-1) * (int)tau * (int)sigma;
    static constexpr int upper_bound = tau * sigma;

    sample();
    int rand() { return ky._rand(); }

protected:
    class knuth_yao {
    public:
        knuth_yao();
        int _rand();
        void initialize();
    private:
        uint8_t prob[upper_bound - lower_bound + 1][lambda] __attribute__((aligned(32)));
        uint8_t weight[lambda] __attribute__((aligned(32)));    // assume upper_bound < 128 for AVX2
    };
    static knuth_yao ky;

private:

};

}

#endif //TBLIB_SAMPLE_HPP
