#ifndef PWHD_BENCH_PARAMS_H
#define PWHD_BENCH_PARAMS_H

struct BenchmarkParams {

    bool verbose = true;
    bool verification = false;
    int repeats = 1;

    bool avoidTestDistruption = true;
    double med2minGuard = 1.030;
    double max2minGuard = 1.150;

};

#endif //PWHD_BENCH_PARAMS_H
