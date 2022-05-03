#include "bench.h"
#include <chrono>

typedef struct _bench_st
{
    /** Stores the time measured before the execution of the benchmark. */
    std::chrono::steady_clock::time_point before;
    /** Stores the time measured after the execution of the benchmark. */
    std::chrono::steady_clock::time_point after;
    /** Stores the sum of timings for the current benchmark. */
    long double total, result;
} bench_ctx;

static bench_ctx g_bench;

void bench_reset(void)
{
    g_bench.total = 0;
}

void bench_before(void)
{
    g_bench.before = std::chrono::steady_clock::now();
}

void bench_after(void)
{
    g_bench.after = std::chrono::steady_clock::now();
    auto diff = g_bench.after - g_bench.before;
    g_bench.total += (long double)(std::chrono::duration_cast<std::chrono::microseconds>(diff).count()); // us
}

void bench_compute(uint32_t benches)
{
    g_bench.result = g_bench.total / benches;
}

void bench_print(void)
{
    std::cout << g_bench.result << " us";

    if (g_bench.result < 0)
    {
        std::cout << " (overflow or bad overhead estimation)" << std::endl;
    }
    else
    {
        std::cout << std::endl;
    }
}
