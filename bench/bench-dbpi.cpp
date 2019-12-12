#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"
#include "../SolverFactory.h"

void logResults(ostream &outStream, const BenchmarkParams &bParams, const ExperimentParams &xParams,
                const vector<double> &times);

using namespace std;

template<typename t_packed>
void benchmark(DBPISolver* solver, BenchmarkParams &bParams, ExperimentParams &xParams) {
    if (bParams.verbose) cout << "Generation of data..." << std::endl;
    uint16_t noOfPacked = ((xParams.m - 1) / xParams.bits_per_packed) + 1;
    t_packed* dataArray = new t_packed[xParams.d * noOfPacked];
    t_packed (&data)[noOfPacked] = reinterpret_cast<t_packed (&)[noOfPacked]>(*dataArray);
    getRandomValues((t_packed*) data, xParams.m, xParams.d, xParams.bits_per_packed, xParams.ones_in_promiles);

    if (bParams.verbose) {
        cout << "Solving... ";
    }

    double buildTime = 0;

    vector<double> times;
    int brute = 0;
    for(int i = 0; i < bParams.repeats; i++) {
        cleanCache();
        time_checkpoint();
        brute += solver->findSimilarSequences((t_packed*) data).size();
        times.push_back(time_millis());
    }
    std::sort(times.begin(), times.end());
    ofstream fout("dbpit_res.txt", ios::out | ios::binary | ios::app);

    logResults(fout, bParams, xParams, times);
    if (bParams.verbose) {
        cout << std::endl << "total time [ms]; algID; m; d; k; ones [%]; bits_packed";
        if (bParams.repeats > 1)
            cout << "; repeats; max/min time [ms]";
        cout <<  std::endl;
    }
    logResults(cout, bParams, xParams, times);

    cout << brute << std::endl;
    if (bParams.verification) cout << "Veryfication uninmplemented :(";

    delete(dataArray);
    if (bParams.verbose) cout << "The end..." << std::endl;
}

void logResults(ostream& outStream, const BenchmarkParams &bParams, const ExperimentParams &xParams,
                const vector<double> &times) {
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    outStream << medianTime << "\t" << xParams.solverID << "\t" << xParams.m << "\t" << xParams.d << "\t" << xParams.k << "\t"
         << xParams.ones_in_promiles << "\t" << (int) xParams.bits_per_packed;
    if (bParams.repeats > 1)
        outStream << "\t" << bParams.repeats << "\t" << maxTime << "\t" << minTime << "\t";
    outStream << endl;
}


int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;

    parseArgs(argc, argv, bParams, xParams);

    DBPISolver* solver = getSolverInstance(xParams);

    benchmark<uint16_t>(solver, bParams, xParams);

    delete(solver);

    return 0;
}

