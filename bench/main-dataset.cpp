#include "../utils/helper.h"
#include "../utils/testdata.h"
#include "cli.h"
#include "../PwHammDistAlgorithmFactory.h"
#include "bench.h"

using namespace std;

void scanDNAdataset(ifstream &src, ExperimentParams &xParams, bool mapSymbols2Values) {
    string seq;
    int length;
    bool alphaPresent[UINT8_MAX] = { false };
    while (getline( src, seq )) {
        xParams.d++;
        for (length = 0; length < seq.length(); length++) {
            if (!isalpha(seq[length]))
                break;
            if (!alphaPresent[seq[length]]) {
                alphaPresent[seq[length]] = true;
                xParams.alphabetSizeUpperBound++;
            }
        }
        if (xParams.m == 0)
            xParams.m = length;
        if (xParams.m != length) {
            fprintf(stderr, "ERROR: variable length of sequences.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (mapSymbols2Values) {
        uint8_t value = 0;
        for (int i = 0; i < UINT8_MAX; i++) {
            if (alphaPresent[i])
                xParams.symbol2value[i] = value++;
        }
    } else {
        for (int i = 0; i < UINT8_MAX; i++)
            xParams.symbol2value[i] = alphaPresent[i]?i:0;
    }
}

void scanIntegerDataset(ifstream &src, ExperimentParams &xParams) {
    string seq;
    int value;
    int length;
    while (getline( src, seq )) {
        xParams.d++;
        length = 0;
        istringstream seqStr( seq );
        while( seqStr >> value ) {
            if (value < 0) {
                fprintf(stderr, "ERROR: unsupported negative values in sequences: %d.\n", value);
                exit(EXIT_FAILURE);
            }
            length++;
            if (xParams.alphabetSizeUpperBound < value + 1)
                xParams.alphabetSizeUpperBound = value + 1;
        }
        if (xParams.m == 0)
            xParams.m = length;
        if (xParams.m != length) {
            fprintf(stderr, "ERROR: variable length of sequences.\n");
            exit(EXIT_FAILURE);
        }
    }
}

void scanSequences(ifstream &src, ExperimentParams &xParams) {
    src.clear();
    src.seekg(0, src.beg);
    xParams.dnaDataMode = isalpha((char) src.peek());
    if (xParams.dnaDataMode) {
        scanDNAdataset(src, xParams, false);
    } else {
        scanIntegerDataset(src, xParams);
    }
}

void loadDNAsequences(ifstream &src, uint8_t* sequences, ExperimentParams &xParams) {
    src.clear();
    src.seekg(0, src.beg);
    string seq;
    uint8_t* cur = sequences;
    while (getline( src, seq )) {
        for(int i = 0; i < xParams.m; i++)
            cur[i] = xParams.symbol2value[seq[i]];
        cur += xParams.bytesPerSequence;
    }
}

void loadIntegerSequences(ifstream &src, uint8_t* sequences, ExperimentParams &xParams) {
    src.clear();
    src.seekg(0, src.beg);
    string seq;
    int value;
    uint8_t* cur = sequences;
    while (getline( src, seq )) {
        istringstream seqStr( seq );
        for(int i = 0; i < xParams.m; i++) {
            seqStr >> value;
            void* dest = cur + (i * xParams.bytesPerElement);
            switch(xParams.bytesPerElement) {
                case 1: *((uint8_t*) dest) = value;
                    break;
                case 2: *((uint16_t*) dest) = value;
                    break;
                case 4: *((uint32_t*) dest) = value;
                    break;
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
        cur += xParams.bytesPerSequence;
    }
}

uint8_t* loadSequences(BenchmarkParams &bParams, ExperimentParams &xParams) {
    const char* srcFile = xParams.datasetFileName.c_str();
    ifstream src(srcFile, ios_base::in | ios_base::binary);
    if (src.fail()) {
        fprintf(stderr, "cannot open reads file %s\n", srcFile);
        exit(EXIT_FAILURE);
    }
    if (bParams.verbose) cout << "Scanning dataset " << xParams.datasetFileName << " ..." << endl;
    scanSequences(src, xParams);
    if (bParams.verbose) {
        cout << "m = " << xParams.m << "; ";
        cout << "d = " << xParams.d << "; ";
        cout << "sigma = " << xParams.alphabetSizeUpperBound << endl;
        cout << "Loading data..." << std::endl;
    }
    xParams.bytesPerElement =  xParams.alphabetSizeUpperBound - 1 <= UINT8_MAX?1:(xParams.alphabetSizeUpperBound - 1 <= UINT16_MAX?2:4);
    xParams.bytesPerSequence = (int) xParams.m * xParams.bytesPerElement;
    if (xParams.alignSequencesTo256bits)
        xParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(xParams.bytesPerSequence,
                ExperimentParams::ALINGMENT_IN_BYTES) * ExperimentParams::ALINGMENT_IN_BYTES;
    uint8_t* sequences = new uint8_t[(size_t) xParams.d * xParams.bytesPerSequence]();
    if (xParams.dnaDataMode) {
        loadDNAsequences(src, sequences, xParams);
    } else {
        loadIntegerSequences(src, sequences, xParams);
    }
    return sequences;
}

int main(int argc, char *argv[]) {

    BenchmarkParams bParams;
    ExperimentParams xParams;
    xParams.enableDatasetMode();

    parseArgs(argc, argv, bParams, xParams);

    uint8_t* sequences = loadSequences(bParams, xParams);
    PwHammDistAlgorithm* algorithm = getPwHammDistAlgorithmInstance(xParams);

    benchmark(sequences, algorithm, bParams, xParams);

    delete(algorithm);
    delete(sequences);

    return 0;
}

