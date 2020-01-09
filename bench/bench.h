#ifndef DBPIT_BENCH_H
#define DBPIT_BENCH_H

#include "../SolverFactory.h"
#include "params.h"

static const char *const BENCH_LOG_FILENAME = "dbpit_res.txt";

void benchmark(uint8_t* sequences, DBPISolver* solver, BenchmarkParams &bParams, ExperimentParams &xParams);

#endif //DBPIT_BENCH_H
