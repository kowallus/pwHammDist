#include "bench.h"

void resultsToStream(ostream &outStream, PwHammDistAlgorithm* algorithm, const BenchmarkParams &bParams, const ExperimentParams &xParams,
                     const vector<double> &times, bool isErrorFound);

void logResults(PwHammDistAlgorithm* algorithm, const BenchmarkParams &bParams, const ExperimentParams &xParams,
        vector<double> &times, bool isErrorFound);

bool isTestDisrupted(vector<double> times, const BenchmarkParams &bParams);

bool isAlgorithmResultValid(uint8_t *sequences, PwHammDistAlgorithm *algorithm, ExperimentParams &xParams) {
    if (xParams.verbose) cout << "Verification... " << endl;
    auto res = algorithm->findSimilarSequences(sequences);
    ExperimentParams mxParams = xParams;
    mxParams.setModelAlgorithm();
    PwHammDistAlgorithm *const modelAlgorithm = getPwHammDistAlgorithmInstance(mxParams);
    auto modelRes = modelAlgorithm->findSimilarSequences(sequences);
    delete(modelAlgorithm);
    sort(res.begin(), res.end());
    sort(modelRes.begin(), modelRes.end());
    bool verifyRes = res == modelRes;
    if (xParams.verbose) {
        if (verifyRes) {
            cout << "OK!" << endl;
        } else {
            cout << "ERROR!" << endl;
            for (size_t i = 0; i < res.size(); i++) {
                if (res[i] != modelRes[i]) {
                    if (res[i] < modelRes[i]) {
                        cout << "An example of false positive: " << res[i].first << " and " << res[i].second << endl;
                    } else {
                        cout << "An example of false negative: " << modelRes[i].first << " and " << modelRes[i].second
                             << endl;
                    }
                    break;
                }
            }
        }
    }
    return verifyRes;
}

void benchmark(uint8_t* sequences, PwHammDistAlgorithm* algorithm, BenchmarkParams &bParams, ExperimentParams &xParams) {

    if (bParams.verbose) cout << "Solving... " << endl;

    vector<double> times;
    uint64_t cummPairsMatched = 0, cummPreStageTimeInUsec;
    bool testAccepted = true;
    do {
        times.clear();
        cummPairsMatched = 0;
        cummPreStageTimeInUsec = 0;
        for (int i = 0; i < bParams.repeats; i++) {
            cleanCache();
            xParams.resetStats();
            time_checkpoint();
            xParams.pairsFound = algorithm->findSimilarSequences(sequences).size();
            times.push_back(time_micros());
            cummPairsMatched += xParams.pairsFound;
#ifdef XP_STATS
            algorithm->cummulateStats(xParams);
            cummPreStageTimeInUsec += xParams.preStageTimeInUsec;
#endif
        }
        if (bParams.avoidTestDistruption)
            testAccepted = !isTestDisrupted(times, bParams);
        if (bParams.verbose and !testAccepted) cout << "d";
    } while (testAccepted);
    bool errorFound = false;
    if (bParams.verification) {
        errorFound = !isAlgorithmResultValid(sequences, algorithm, xParams);
    }
    xParams.pairsFound = (cummPairsMatched / bParams.repeats);
    xParams.preStageTimeInUsec = (cummPreStageTimeInUsec / bParams.repeats);
    logResults(algorithm, bParams, xParams, times, errorFound);
    if (bParams.verbose) cout << std::endl << "pairs count (average per repeat): " << xParams.pairsFound << std::endl;
    if (bParams.verbose) cout << "The end..." << std::endl;
}

bool isTestDisrupted(vector<double> times, const BenchmarkParams &bParams) {
    sort(times.begin(), times.end());
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    return minTime * bParams.max2minGuard >= maxTime || minTime * bParams.med2minGuard >= medianTime;
}

void logResults(PwHammDistAlgorithm* algorithm, const BenchmarkParams &bParams, const ExperimentParams &xParams,
        vector<double> &times, bool isErrorFound) {
    sort(times.begin(), times.end());
    ofstream fout(BENCH_LOG_FILENAME, ios_base::out | ios_base::binary | ios_base::app);

    resultsToStream(fout, algorithm, bParams, xParams, times, isErrorFound);
    if (bParams.verbose) {
        cout << endl << "time[ms]\t                    algID\t    m\t    d\tsigma\t    k" <<
            (xParams.isOnesInPromilesEnabled()?"\tones[%]":"") <<
            (xParams.isBitsPerPackedEnabled()?"\tbits_packed":"");
#ifdef XP_STATS
        cout << "  matched pairs\t verified pairs\tfilter accepted\tfilter rejected\t  ignored pairs"
                "\t  pre-time [us]";
#endif
        if (bParams.repeats > 1)
            cout << "\trepeats\tmax/min time [us]";
        cout <<  endl;
    }
    resultsToStream(cout, algorithm, bParams, xParams, times, isErrorFound);
}

string microSecToMillis(double medianTime, uint8_t minDigits2show) {
    string millisTime = toString(medianTime);
    if (millisTime.length() <= 3) {
        while(millisTime.length() < 3)
            millisTime.insert(0, "0");
        millisTime.insert(0, "0.");
    } else {
        if (toString(medianTime + 1).length() > millisTime.length())
            millisTime = toString(++medianTime);
        if (minDigits2show + 3 <= millisTime.length()) {
            if (millisTime[millisTime.length() - 3] >= '5')
                medianTime = medianTime / 1000 + 1;
            else
                medianTime /= 1000;
            millisTime = toString(medianTime);
        } else {
            millisTime.insert(millisTime.length() - 3, ".");
            bool inc = false;
            while (millisTime.length() > minDigits2show + 1) {
                inc = millisTime[millisTime.length() - 1] >= '5';
                millisTime.pop_back();
            }
            for (int i = millisTime.length(); i-- > 0 && inc;) {
                if (millisTime[i] == '.')
                    continue;
                if (millisTime[i] == '9')
                    millisTime[i] = '0';
                else {
                    millisTime[i]++;
                    inc = false;
                }
            }
            if (inc)
                millisTime.insert(0, "1");
        }
    }
    return alignRight(millisTime, 8);
}

void resultsToStream(ostream &outStream, PwHammDistAlgorithm* algorithm, const BenchmarkParams &bParams,
                     const ExperimentParams &xParams, const vector<double> &times, bool isErrorFound) {
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    outStream << microSecToMillis(medianTime, 3) << "\t" << alignRight(algorithm->getName(), 25) <<
              "\t" << alignRight(toString(xParams.m), 5) << "\t" << alignRight(toString(xParams.d), 5) <<
              "\t" << alignRight(toString(xParams.alphabetSizeUpperBound), 5) << "\t" << alignRight(toString(xParams.k), 5);
    if (xParams.isOnesInPromilesEnabled())
        outStream << "\t" << alignRight(toString(xParams.onesInPromiles),7);
    if (xParams.isBitsPerPackedEnabled())
        outStream << "\t" << alignRight(toString(xParams.bitsPerPacked), 11);
#ifdef XP_STATS
    outStream << "\t" << alignRight(toString(xParams.pairsFound), 15)
        << "\t" << alignRight(toString(xParams.pairsVerified), 15)
        << "\t" << alignRight(toString(xParams.pairsFilterAccepted), 15)
        << "\t" << alignRight(toString(xParams.pairsFilterRejected), 15)
        << "\t" << alignRight(toString(xParams.pairsIgnored), 15)
        << "\t" << alignRight(toString(xParams.preStageTimeInUsec), 15);
#endif
    if (bParams.repeats > 1)
        outStream << "\t" << alignRight(toString(bParams.repeats), 7) << "\t" << (int) maxTime << "\t" << (int) minTime << "\t";
    if (isErrorFound)
        outStream << "ERROR";
    outStream << endl;
}
