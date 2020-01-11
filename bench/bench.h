#ifndef PWHD_BENCH_H
#define PWHD_BENCH_H

#include "../PwHammDistAlgorithmFactory.h"
#include "params.h"

static const char *const BENCH_LOG_FILENAME = "pwhd_res.txt";

void benchmark(uint8_t* sequences, PwHammDistAlgorithm* algorithm, BenchmarkParams &bParams, ExperimentParams &xParams);

#endif //PWHD_BENCH_H
