#ifndef TBLIB_BENCH_H
#define TBLIB_BENCH_H

#include <iostream>
#include <iomanip>

/**
 * Runs a new benchmark.
 *
 * @param[in] LABEL			- the label for this benchmark.
 */
#define BENCH_START(_LABEL, _BENCHS)                                        \
   bench_reset();                                                           \
   std::cout << "BENCH: " << std::setw(48) << std::left << _LABEL << " ";   \
   for (int _b = 0; _b < _BENCHS; _b++)                                     \
   {

/**
 * Prints the average timing of each execution in the chosen metric.
 */
#define BENCH_FINAL(_BENCHS, _ROUNDS)   \
   }                                    \
   bench_compute(_BENCHS * _ROUNDS);    \
   bench_print();

/**
 * Measures the time of one execution and adds it to the benchmark total.
 *
 * @param[in] FUNCTION		- the function executed.
 */
#define BENCH_ITEM(_FUNCTION, _ROUNDS)  \
   _FUNCTION;                           \
   bench_before();                      \
   for (int _r = 0; _r < _ROUNDS; _r++) \
   {                                    \
      _FUNCTION;                        \
   }                                    \
   bench_after();

/*============================================================================*/
/* Function definitions                                                       */
/*============================================================================*/

/**
 * Resets the benchmark data.
 *
 * @param[in] label			- the benchmark label.
 */
void bench_reset(void);

/**
 * Measures the time before a benchmark is executed.
 */
void bench_before(void);

/**
 * Measures the time after a benchmark was started and adds it to the total.
 */
void bench_after(void);

/**
 * Computes the mean elapsed time between the start and the end of a benchmark.
 *
 * @param benches			- the number of executed benchmarks.
 */
void bench_compute(uint32_t benches);

/**
 * Prints the last benchmark.
 */
void bench_print(void);

#endif //TBLIB_BENCH_H
