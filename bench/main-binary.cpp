#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"
#include "../SolverFactory.h"
#include "bench.h"

using namespace std;

uint8_t* generateSequences(BenchmarkParams &bParams, ExperimentParams &xParams) {
    if (bParams.verbose) cout << "Generation of data..." << std::endl;
    xParams.bytesPerSequence = ceilDivisionBySmallInteger(xParams.m, xParams.bitsPerPacked)
            * ceilDivisionBySmallInteger(xParams.bitsPerPacked, 8);
    if (xParams.alignSequences)
        xParams.bytesPerSequence = ceilDivisionBySmallInteger(xParams.bytesPerSequence, 8) * 8;
    uint8_t* sequences = new uint8_t[xParams.d * xParams.bytesPerSequence]();
    getRandomValues(sequences, xParams);
    return sequences;
}

int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;
    xParams.enableBinaryMode();

    parseArgs(argc, argv, bParams, xParams);

    uint8_t* sequences = generateSequences(bParams, xParams);
    DBPISolver* solver = getSolverInstance(xParams);

    benchmark(sequences, solver, bParams, xParams);

    delete(solver);
    delete(sequences);

    return 0;
}

