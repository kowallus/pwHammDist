#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"
#include "../PwHammDistAlgorithmFactory.h"
#include "bench.h"

using namespace std;

uint8_t* generateSequences(BenchmarkParams &bParams, ExperimentParams &xParams) {
    if (bParams.verbose) cout << "Generation of data..." << std::endl;
    xParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(xParams.m, xParams.bitsPerPacked)
            * ceilDivisionBySmallInteger(xParams.bitsPerPacked, 8);
    if (xParams.alignSequences)
        xParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(xParams.bytesPerSequence, 16) * 16;
    uint8_t* sequences = new uint8_t[(size_t) xParams.d * xParams.bytesPerSequence]();
    getRandomValues(sequences, xParams);
    return sequences;
}

int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;
    xParams.enableBinaryMode();

    parseArgs(argc, argv, bParams, xParams);

    uint8_t* sequences = generateSequences(bParams, xParams);
    PwHammDistAlgorithm* algorithm = getPwHammDistAlgorithmInstance(xParams);

    benchmark(sequences, algorithm, bParams, xParams);

    delete(algorithm);
    delete(sequences);

    return 0;
}

