#include "cli.h"

#include <unistd.h>
#include <iostream>

void parseArgs(int argc, char *argv[], BenchmarkParams &bParams, ExperimentParams &xParams) {

    int opt; // current option

    std::string paramsStr = std::string(xParams.isBitsPerPackedEnabled()?"b:":"") +
            std::string(xParams.isInBinaryMode()?"":"ciI") +
            std::string("r:a:pPngAvq?");

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
            case 'n':
                xParams.shortCircuitMode = false;
                break;
            case 'P':
                xParams.pivotsElectionMode = true;
            case 'p':
                xParams.pivotsFilterMode = true;
                break;
            case 'g':
                xParams.groupedBruteMode = true;
                break;
            case 'c':
                xParams.compactMode = true;
                break;
            case 'I':
                xParams.lazyInterleaveBitsMode = true;
            case 'i':
                xParams.interleaveBitsMode = true;
                break;
            case 'A':
                xParams.alignSequencesTo256bits = false;
                break;
            case 'q':
                bParams.verbose = false;
                xParams.verbose = false;
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
                fprintf(stderr, "Usage: %s [-a algorithmID] [-n] [-p] [-P] [-g] \n"
                                "\t\t[-r noOfRepeats] [-v] [-q] [-A] ", argv[0]);
                if (!xParams.isInBinaryMode())
                    fprintf(stderr, "[-c] [-i] [-I]");
                if (xParams.isBitsPerPackedEnabled())
                    fprintf(stderr, "[-b bitsPerPacked] ");
                if (xParams.isOnesInPromilesEnabled())
                    fprintf(stderr, "onesInPromiles ");
                if (xParams.isInDatasetMode())
                    fprintf(stderr, "datasetFileName ");
                else
                    fprintf(stderr, "m d ");
                fprintf(stderr, "k\n\n");
                fprintf(stderr, "algorithm ID : core algorithm name\n");
                for(pair<string, PwHammDistAlgorithmFactory*> p: pwHammDistAlgorithmTypesMap) {
                    fprintf(stderr, "%s : %s\n", p.first.c_str(), p.second->getAlgorithmName().c_str());
                }
                fprintf(stderr, "\n-p pivot filter mode (or -P with pivot election)"
                                "\n-n disable short-circuit"
                                "\n-g sequences matched in groups (for bf algorithm)");
                if (!xParams.isInBinaryMode()) {
                    fprintf(stderr, "\n-c compact mode processing");
                    fprintf(stderr, "\n-i interleaved bits mode (or -I with lazy evaluation)");
                }
                fprintf(stderr, "\n\n-A ignore aligning sequences to 256-bit");
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

    if (xParams.interleaveBitsMode && xParams.compactMode) {
        fprintf(stderr, "Error: modes compact and interleaved bits cannot be used together.\n");

        exit(EXIT_FAILURE);
    }

    if (xParams.lazyInterleaveBitsMode && xParams.pivotsFilterMode) {
        fprintf(stderr, "Error: mode pivots and lazy interleaved bits cannot be used together.\n");

        exit(EXIT_FAILURE);
    }

    if (xParams.lazyInterleaveBitsMode && !xParams.shortCircuitMode) {
        fprintf(stderr, "Error: lazy interleaved bits cannot be used without short-circuit mode.\n");

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