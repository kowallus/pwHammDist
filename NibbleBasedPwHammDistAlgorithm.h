#ifndef PWHD_NIBBLEBASEDPWHAMMDISTALGORITHM_H
#define PWHD_NIBBLEBASEDPWHAMMDISTALGORITHM_H

#include "PwHammDistAlgorithmBase.h"

const string NIBBLE_NAIVE_BRUTE_FORCE_ID = "nnbf";
const string NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID = "nsbf";
const string GROUPED_NIBBLE_NAIVE_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + "nnbf";
const string GROUPED_NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + "nsbf";

template <bool naive, typename uint = uint8_t>
class NibbleBrutePwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    ExperimentParams nxParams;

    uint8_t* nibblify(const uint8_t* sequences) {
        nxParams.m = xParams.m;
        nxParams.d = xParams.d;
        nxParams.k = xParams.k;
        uint8_t *nibbles;
        if (sizeof(uint) == sizeof(uint8_t) && xParams.alphabetSize <= 8) {
            nxParams.bytesPerSequence = xParams.bytesPerSequence / 2;
            nibbles = new uint8_t[nxParams.d * nxParams.bytesPerSequence * 2]();
            uint8_t* x = (uint8_t*) sequences;
            uint8_t *y = nibbles;
            uint8_t *z = nibbles + nxParams.d * nxParams.bytesPerSequence;
            const uint32_t nibblePackedBytes = xParams.d * xParams.bytesPerSequence / 2;
            for(int i = 0; i < nibblePackedBytes; i++) {
                *y = *(x++);
                *y += *(x++)*16;
                *(z++) = *(y++) + 128 + 8;
            }
        } else if (sizeof(uint) == sizeof(uint16_t)) {
            nxParams.bytesPerSequence = xParams.bytesPerSequence;
            nibbles = new uint8_t[nxParams.d * nxParams.bytesPerSequence]();
            uint16_t* x = (uint16_t*) sequences;
            uint16_t* z = (uint16_t*) nibbles;
            const uint32_t length = xParams.d * xParams.bytesPerSequence / 2;
            for(int i = 0; i < length; i++) {
                *(z++) = *(x++) + 32768;
            }
        } else {
                fprintf(stderr, "ERROR: unsupported alphabet size for nibblification.\n");
                exit(EXIT_FAILURE);
        }
        if (xParams.verbose) cout << "nibbled... " << " (" << time_millis() << " msec)" << endl;
        return nibbles;
    }

    inline bool testNibblesSimilarity(const uint8_t* startX, const uint8_t* startY, uint16_t i, uint16_t j) {
        uint8_t* x = (uint8_t*) startX + (size_t) i * nxParams.bytesPerSequence;
        uint8_t* augY = (uint8_t*) startY + (size_t) j * nxParams.bytesPerSequence;
        return testNibblesSimilarity(x, augY);
    };

    inline bool testNibblesSimilarity(const void* nib1, const void* augNib2) {
        uint16_t dist;
        if (sizeof(uint) == sizeof(uint8_t)) {
            dist = naive ? hammingDistanceAugmentedNibble((uint64_t *) nib1, (uint64_t *) augNib2,
                                                          nxParams.bytesPerSequence * 2) :
                   hammingDistanceAugmentedNibble((uint64_t *) nib1, (uint64_t *) augNib2,
                                                  nxParams.bytesPerSequence * 2, xParams.k);
        } else {
            dist = naive ? hammingDistanceAugmented16bit((uint64_t *) nib1, (uint64_t *) augNib2,
                                                          nxParams.bytesPerSequence / 2) :
                   hammingDistanceAugmented16bit((uint64_t *) nib1, (uint64_t *) augNib2,
                                                  nxParams.bytesPerSequence / 2, xParams.k);
        }
        return dist <= xParams.k;
    }

public:
    NibbleBrutePwHammDistAlgorithm(ExperimentParams &xParams): PwHammDistAlgorithm(xParams) {
        if (xParams.bytesPerSequence % 16 != 0) {
            fprintf(stderr, "ERROR: brute algorithm does not support 128-bit unaligned data.\n");
            exit(EXIT_FAILURE);
        }
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        uint8_t* nibbles = nibblify(sequences);
        vector<pair<uint16_t, uint16_t>> res;
        uint8_t *x;
        uint8_t *startY;
        if (sizeof(uint) == sizeof(uint8_t) && xParams.alphabetSize <= 8) {
            x = nibbles;
            startY = nibbles + nxParams.d * nxParams.bytesPerSequence;
        } else {
            x = (uint8_t*) sequences;
            startY = nibbles;
        }
        for(int i = 0; i < nxParams.d - 1; i++) {
            uint8_t *augY = startY + (i + 1) * nxParams.bytesPerSequence;
            for (int j = i + 1; j < nxParams.d; j++) {
                if (testNibblesSimilarity(x, augY)) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
                augY += nxParams.bytesPerSequence;
            }
            x += nxParams.bytesPerSequence;
        }
        delete[] nibbles;
        return res;
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        uint8_t* nibbles = nibblify(sequences);
        uint8_t *x;
        uint8_t *startY;
        if (sizeof(uint) == sizeof(uint8_t) && xParams.alphabetSize <= 8) {
            x = nibbles;
            startY = nibbles + nxParams.d * nxParams.bytesPerSequence;
        } else {
            x = (uint8_t*) sequences;
            startY = nibbles;
        }
        vector<pair<uint16_t, uint16_t>> res;
        for (pair<uint16_t, uint16_t> pair: pairs) {
            if (testNibblesSimilarity(x, startY, pair.first, pair.second)) {
                res.push_back(pair);
            }
        }
        delete[] nibbles;
        return res;
    }

    inline bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) {
        uint8_t* x = (uint8_t*) sequences + (size_t) i * nxParams.bytesPerSequence;
        uint8_t* y = (uint8_t*) sequences + (size_t) j * nxParams.bytesPerSequence;
        return testSequencesSimilarity(x, y);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        uint16_t dist;
        dist = naive?hammingDistance((uint*) seq1, (uint*) seq2, xParams.m):
               hammingDistance((uint*) seq1, (uint*) seq2, xParams.m, xParams.k);
        return dist <= xParams.k;
    }

    string getName() { return (naive?NIBBLE_NAIVE_BRUTE_FORCE_ID:NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID); }
};

#endif //PWHD_NIBBLEBASEDPWHAMMDISTALGORITHM_H
