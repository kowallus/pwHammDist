#ifndef DBPIT_CLI_H
#define DBPIT_CLI_H

#include "params.h"
#include "../xp-params.h"

void parseArgs(int argc, char *argv[], BenchmarkParams &bParams, ExperimentParams &xParams);

#endif //DBPIT_CLI_H
