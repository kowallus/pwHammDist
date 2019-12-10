#include "../utils/helper.h"
#include "../utils/testdata.h"
#include <unistd.h>

using namespace std;

bool verbose = true;
bool verification = false;
int repeats = 1;
fstream fout("dbpit_res.txt", ios::out | ios::binary | ios::app);

template<typename t_packed>
void benchmark(uint16_t m, uint16_t d, uint16_t k, uint8_t bits_packed, uint16_t ones_in_promiles) {
    if (verbose) cout << "Generation of data..." << std::endl;
    uint16_t noOfPacked = ((m - 1) / bits_packed) + 1;
    t_packed* dataArray = new t_packed[d * noOfPacked];
    t_packed (&data)[noOfPacked] = reinterpret_cast<t_packed (&)[noOfPacked]>(*dataArray);
    getRandomValues((t_packed*) data, m, d, bits_packed, ones_in_promiles);

    if (verbose) {
        cout << "Solving... ";
    }

    double buildTime = 0;

    vector<double> times;
    int brute = 0;
    for(int i = 0; i < repeats; i++) {
        cleanCache();
        time_checkpoint();
//        brute += dbpi_trasposed(values, m, d, k, bits_packed);
        times.push_back(time_millis());
    }
    std::sort(times.begin(), times.end());
    double maxTime = times[repeats - 1];
    double medianTime = times[times.size()/2];
    double minTime = times[0];
    if (verbose) cout << std::endl << "total time [ms]; m; d; k; ones [%]; bits_packed; max/min time [ms]" << std::endl;
    cout << medianTime << "\t" << m << "\t" << d << "\t" << k << "\t" << ones_in_promiles << "\t" << (int) bits_packed
         << "\t" << maxTime << "\t" << minTime << "\t" << std::endl;
    fout << medianTime << "\t" << m << "\t" << d << "\t" << k << "\t" << ones_in_promiles << "\t" << (int) bits_packed
         << "\t" << maxTime << "\t" << minTime << "\t" << std::endl;
    cout << brute << std::endl;
    if (verification) cout << "Veryfication uninmplemented :(";

    delete(dataArray);
    if (verbose) cout << "The end..." << std::endl;
}


int main(int argc, char *argv[]) {

    int opt; // current option

    while ((opt = getopt(argc, argv, "r:vq?")) != -1) {
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'v':
                verification = true;
                break;
            case 'r':
                repeats = atoi(optarg);
                if (repeats <= 0) {
                    fprintf(stderr, "%s: Expected number of repeats >=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-r noOfRepeats] [-v] [-q] ones_in_promiles m d k\n\n",
                        argv[0]);
                fprintf(stderr, "\n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind != (argc - 4)) {
        fprintf(stderr, "%s: Expected 4 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    uint16_t ones_in_promiles = atoi(argv[optind++]);
    uint16_t m = atoi(argv[optind++]);
    uint16_t d = atoi(argv[optind++]);
    uint16_t k = atoi(argv[optind++]);

    uint8_t bits = 12;

    benchmark<uint16_t>(m, d, k, bits, ones_in_promiles);

    return 0;
}

