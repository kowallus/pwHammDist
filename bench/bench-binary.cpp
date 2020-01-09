#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"
#include "../SolverFactory.h"

void resultsToStream(ostream &outStream, DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams,
                     const vector<double> &times);

void logResults(DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams, vector<double> &times);

using namespace std;

void benchmark(DBPISolver* solver, BenchmarkParams &bParams, ExperimentParams &xParams) {
    if (bParams.verbose) cout << "Generation of data..." << std::endl;
    xParams.bytesPerSequence = ceilDivisionBySmallInteger(xParams.m, xParams.bitsPerPacked)
            * ceilDivisionBySmallInteger(xParams.bitsPerPacked, 8);
    if (xParams.alignSequences)
        xParams.bytesPerSequence = ceilDivisionBySmallInteger(xParams.bytesPerSequence, 8) * 8;
    uint8_t* sequences = new uint8_t[xParams.d * xParams.bytesPerSequence]();
    getRandomValues(sequences, xParams);

    if (bParams.verbose) cout << "Solving... " << endl;

    vector<double> times;
    int brute = 0;
    for(int i = 0; i < bParams.repeats; i++) {
        cleanCache();
        time_checkpoint();
        brute += solver->findSimilarSequences(sequences).size();
        times.push_back(time_millis());
    }
    logResults(solver, bParams, xParams, times);

    cout << std::endl << "check: " << (brute / bParams.repeats) << std::endl;
    if (bParams.verification) cout << "Veryfication uninmplemented :(";

    delete(sequences);
    if (bParams.verbose) cout << "The end..." << std::endl;
}

void logResults(DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams, vector<double> &times) {
    sort(times.begin(), times.end());
    ofstream fout("dbpit_res.txt", ios_base::out | ios_base::binary | ios_base::app);

    resultsToStream(fout, solver, bParams, xParams, times);
    if (bParams.verbose) {
        cout << endl << "total time [ms]; algID; m; d; k; ones [%]; bits_packed";
        if (bParams.repeats > 1)
            cout << "; repeats; max/min time [ms]";
        cout <<  endl;
    }
    resultsToStream(cout, solver, bParams, xParams, times);
}

void resultsToStream(ostream &outStream, DBPISolver* solver, const BenchmarkParams &bParams,
        const ExperimentParams &xParams, const vector<double> &times) {
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    outStream << medianTime << "\t" << solver->getName() << "\t" << xParams.m << "\t" << xParams.d << "\t" << xParams.k << "\t"
         << xParams.onesInPromiles << "\t" << (int) xParams.bitsPerPacked;
    if (bParams.repeats > 1)
        outStream << "\t" << bParams.repeats << "\t" << maxTime << "\t" << minTime << "\t";
    outStream << endl;
}


int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;
    xParams.enableBinaryMode();

    parseArgs(argc, argv, bParams, xParams);

    DBPISolver* solver = getSolverInstance(xParams);

    benchmark(solver, bParams, xParams);

    delete(solver);

    return 0;
}

