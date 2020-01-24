#include "cli.h"

#include <unistd.h>
#include <iostream>

void parseArgs(int argc, char *argv[], BenchmarkParams &bParams, ExperimentParams &xParams) {

    int opt; // current option

    std::string paramsStr = std::string(xParams.isBitsPerPackedEnabled()?"b:":"") + std::string("r:a:Avq?");

    while ((opt = getopt(argc, argv, paramsStr.c_str())) != -1) {
        switch (opt) {
            case 'b':
                if (atoi(optarg) > 0)
                    xParams.bitsPerPacked = atoi(optarg);
                else {
                    xParams.bitsPerPacked = 1;
                    xParams.binaryMode = false;
                }
                break;
            case 'a':
                xParams.algorithmID = std::string(optarg);
                break;
            case 'A':
                xParams.alignSequencesTo256bits = false;
                break;
            case 'q':
                bParams.verbose = false;
                break;
            case 'v':
                bParams.verification = true;
                break;
            case 'r':
                bParams.repeats = atoi(optarg);
                if (bParams.repeats <= 0) {
                    fprintf(stderr, "%s: Expected number of repeats >=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-a algorithmID] [-r noOfRepeats] [-v] [-q] [-A] ", argv[0]);
                if (xParams.isBitsPerPackedEnabled())
                    fprintf(stderr, "[-b bitsPerPacked] ");
                if (xParams.isOnesInPromilesEnabled())
                    fprintf(stderr, "onesInPromiles ");
                if (xParams.isInDatasetMode())
                    fprintf(stderr, "datasetFileName ");
                else
                    fprintf(stderr, "m d ");
                fprintf(stderr, "k\n\n");
                fprintf(stderr, "algorithm ID : algorithm name\n");
                for(pair<string, PwHammDistAlgorithmFactory*> p: pwHammDistAlgorithmTypesMap) {
                    fprintf(stderr, "%s : %s\n", p.first.c_str(), p.second->getAlgorithmName().c_str());
                }
                fprintf(stderr, "\n-A ignore aligning sequences to 128-bit");
                fprintf(stderr, "\n-v verify results \n-q quiet output (only parameters)\n\n");
                if (xParams.isBitsPerPackedEnabled())
                    fprintf(stderr, "To disable binary mode apply: -b 0\n\n");
                exit(EXIT_FAILURE);
        }
    }

    int fixedArgsCount = (xParams.isInDatasetMode()?1:2) + 1 + (xParams.isOnesInPromilesEnabled()?1:0);

    if (optind != (argc - fixedArgsCount)) {
        fprintf(stderr, "%s: Expected %d arguments after options (found %d)\n", argv[0], fixedArgsCount, argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    if (xParams.isOnesInPromilesEnabled())
        xParams.onesInPromiles = atoi(argv[optind++]);
    if (xParams.isInDatasetMode())
        xParams.datasetFileName = argv[optind++];
    else {
        xParams.m = atoi(argv[optind++]);
        xParams.d = atoi(argv[optind++]);
    }
    xParams.k = atoi(argv[optind++]);
}