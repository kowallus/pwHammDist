#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"

using namespace std;

fstream fout("dbpit_res.txt", ios::out | ios::binary | ios::app);

template<typename t_packed>
void benchmark(BenchmarkParams &bParams, ExperimentParams &xParams) {
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
//        brute += dbpi_trasposed(values, m, d, k, bits_packed);
        times.push_back(time_millis());
    }
    std::sort(times.begin(), times.end());
    double maxTime = times[bParams.repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    if (bParams.verbose) cout << std::endl << "total time [ms]; m; d; k; ones [%]; bits_packed; max/min time [ms]" << std::endl;
    cout << medianTime << "\t" << xParams.m << "\t" << xParams.d << "\t" << xParams.k << "\t"
            << xParams.ones_in_promiles << "\t" << (int) xParams.bits_per_packed
            << "\t" << maxTime << "\t" << minTime << "\t" << std::endl;
    fout << medianTime << "\t" << xParams.m << "\t" << xParams.d << "\t" << xParams.k << "\t"
            << xParams.ones_in_promiles << "\t" << (int) xParams.bits_per_packed
            << "\t" << maxTime << "\t" << minTime << "\t" << std::endl;
    cout << brute << std::endl;
    if (bParams.verification) cout << "Veryfication uninmplemented :(";

    delete(dataArray);
    if (bParams.verbose) cout << "The end..." << std::endl;
}


int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;

    parseArgs(argc, argv, bParams, xParams);

    benchmark<uint16_t>(bParams, xParams);

    return 0;
}

