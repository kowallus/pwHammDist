#include "cli.h"

#include <unistd.h>
#include <iostream>

void parseArgs(int argc, char *argv[], BenchmarkParams &bParams, ExperimentParams &xParams) {

    int opt; // current option

    std::string paramsStr = std::string(xParams.isDisabledBitsPerPacked()?"":"b:") + std::string("r:a:vq?");

    while ((opt = getopt(argc, argv, paramsStr.c_str())) != -1) {
        switch (opt) {
            case 'b':
                xParams.bits_per_packed = atoi(optarg);
                break;
            case 'a':
                xParams.solverID = std::string(optarg);
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
                fprintf(stderr, "Usage: %s [-a algorithmID] [-r noOfRepeats] [-v] [-q] ", argv[0]);
                if (!xParams.isDisabledOnesInPromiles())
                    fprintf(stderr, "ones_in_promiles ");
                fprintf(stderr, "m d k\n\n");
                fprintf(stderr, "-a algorithmIDs: nbf, sbf\n");
                fprintf(stderr, "\n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
                exit(EXIT_FAILURE);
        }
    }

    int fixedArgsCount = 4 - xParams.isDisabledOnesInPromiles()?1:0;

    if (optind != (argc - 4)) {
        fprintf(stderr, "%s: Expected %d arguments after options (found %d)\n", argv[0], fixedArgsCount, argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    if (!xParams.isDisabledOnesInPromiles())
        xParams.ones_in_promiles = atoi(argv[optind++]);
    xParams.m = atoi(argv[optind++]);
    xParams.d = atoi(argv[optind++]);
    xParams.k = atoi(argv[optind++]);

}