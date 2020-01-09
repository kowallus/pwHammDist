#include "bench.h"

void resultsToStream(ostream &outStream, DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams,
                     const vector<double> &times);

void logResults(DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams, vector<double> &times);

void verifySolverResult(const uint8_t *sequences, DBPISolver *solver, ExperimentParams &xParams) {
    if (xParams.solverID != NAIVE_BRUTE_FORCE_ID) {
            cout << "Verification... " << endl;
            auto res = solver->findSimilarSequences(sequences);
            xParams.solverID == NAIVE_BRUTE_FORCE_ID;
            DBPISolver *const modelSolver = getSolverInstance(xParams);
            auto modelRes = modelSolver->findSimilarSequences(sequences);
            delete(modelSolver);
            sort(res.begin(), res.end());
            sort(modelRes.begin(), modelRes.end());
            if (res == modelRes) {
                cout << "OK!" << endl;
            } else {
                cout << "ERROR!" << endl;
                for(size_t i = 0; i < res.size(); i++) {
                    if (res[i] != modelRes[i]) {
                        if (res[i] < modelRes[i]) {
                            cout << "An example of false positive: " << res[i].first << " and " << res[i].second << endl;
                        } else {
                            cout << "An example of false negative: " << modelRes[i].first << " and " << modelRes[i].second << endl;
                        }
                    }
                }
            }
        } else {
            cout << "Model algorithm cannot be verified ";
        }
}

void benchmark(uint8_t* sequences, DBPISolver* solver, BenchmarkParams &bParams, ExperimentParams &xParams) {

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
    if (bParams.verification) {
        verifySolverResult(sequences, solver, xParams);
    }
    if (bParams.verbose) cout << "The end..." << std::endl;
}

void logResults(DBPISolver* solver, const BenchmarkParams &bParams, const ExperimentParams &xParams, vector<double> &times) {
    sort(times.begin(), times.end());
    ofstream fout(BENCH_LOG_FILENAME, ios_base::out | ios_base::binary | ios_base::app);

    resultsToStream(fout, solver, bParams, xParams, times);
    if (bParams.verbose) {
        cout << endl << "time[ms]\t  algID\t    m\t    d\t    k" <<
            (xParams.isOnesInPromilesEnabled()?"\tones[%]":"") <<
            (xParams.isBitsPerPackedEnabled()?"\tbits_packed":"");
        if (bParams.repeats > 1)
            cout << "\trepeats\tmax/min time [ms]";
        cout <<  endl;
    }
    resultsToStream(cout, solver, bParams, xParams, times);
}

void resultsToStream(ostream &outStream, DBPISolver* solver, const BenchmarkParams &bParams,
                     const ExperimentParams &xParams, const vector<double> &times) {
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    outStream << alignRight(toString(medianTime), 8) << "\t" << alignRight(solver->getName(), 7) <<
        "\t" << alignRight(toString(xParams.m), 5) << "\t" << alignRight(toString(xParams.d), 5) <<
        "\t" << alignRight(toString(xParams.k), 5);
    if (xParams.isOnesInPromilesEnabled())
        outStream << "\t" << alignRight(toString(xParams.onesInPromiles),7);
    if (xParams.isBitsPerPackedEnabled())
        outStream << "\t" << alignRight(toString(xParams.bitsPerPacked), 11);
    if (bParams.repeats > 1)
        outStream << "\t" << alignRight(toString(bParams.repeats), 7) << "\t" << maxTime << "\t" << minTime << "\t";
    outStream << endl;
}

