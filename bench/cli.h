#ifndef PWHD_CLI_H
#define PWHD_CLI_H

#include "params.h"
#include "../xp-params.h"
#include "../PwHammDistAlgorithmFactory.h"

void parseArgs(int argc, char *argv[], BenchmarkParams &bParams, ExperimentParams &xParams);

#endif //PWHD_CLI_H
