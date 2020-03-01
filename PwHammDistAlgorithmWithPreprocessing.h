#ifndef PWHD_DATAPREPROCESSORBASE_H
#define PWHD_DATAPREPROCESSORBASE_H

#include "xp-params.h"
#include "PwHammDistAlgorithmBase.h"

const string COLUMNS_SHUFFLER_ID = "shuf";

class PwHammDistAlgorithmWithPreprocessing: public PwHammDistAlgorithm {
private:
    PwHammDistAlgorithm* mainAlgorithm;

    static const uint8_t SHUFFLE_SOLID_BLOCK_IN_ULLS = 2;
    static const uint8_t SHUFFLE_SKIP_SOLID_BLOCKS = 16;

    uint8_t* defaultColumnsShuffle(uint8_t* sequences) {
        uint8_t* shuffle = new uint8_t[xParams.bytesPerSequence * xParams.d];
        uint64_t* shRowPtr = (uint64_t*) shuffle;
        uint64_t* shPtr = shRowPtr;
        uint16_t bytesToShuffle = xParams.binaryMode? (int) xParams.m / xParams.bitsPerPacked
                                  * ceilDivisionBySmallInteger(xParams.bitsPerPacked, 8):
                                  xParams.m * xParams.bytesPerElement;
        for(int i = 0; i < xParams.d; i++) {
            for (int j = 0; j < SHUFFLE_SKIP_SOLID_BLOCKS; j++) {
                uint16_t sePos = j * SHUFFLE_SOLID_BLOCK_IN_ULLS * 8;
                while(sePos + SHUFFLE_SOLID_BLOCK_IN_ULLS * 8 <= bytesToShuffle ) {
                    uint64_t* sePtr = (uint64_t*) (sequences + sePos);
                    for (int k = 0; k < SHUFFLE_SOLID_BLOCK_IN_ULLS; k++)
                        *shPtr++ = *sePtr++;
                    sePos += SHUFFLE_SKIP_SOLID_BLOCKS * SHUFFLE_SOLID_BLOCK_IN_ULLS * 8;
                }
            }
            uint64_t* sePtr = (uint64_t*) sequences + (shPtr - shRowPtr);
            shRowPtr += xParams.bytesPerSequence / 8;
            while(shRowPtr != shPtr)
                *shPtr++ = *sePtr++;
            sequences += xParams.bytesPerSequence;
        }
        return shuffle;
    };

    uint8_t *preprocessSequences(uint8_t *sequences) {
        if (xParams.shuffleColumnsMode) {
            sequences = defaultColumnsShuffle(sequences);
            if (xParams.verbose) cout << "shuffled columns in (" << time_millis() << " msec)" << endl;
        }
        return sequences;
    }

    void postprocess(uint8_t *sequences) {
        if (xParams.shuffleColumnsMode) {
            delete[] sequences;
        }
    }

public:
    PwHammDistAlgorithmWithPreprocessing(ExperimentParams xParams, PwHammDistAlgorithm* mainAlgorithm):
        PwHammDistAlgorithm(xParams), mainAlgorithm(mainAlgorithm) {};

    virtual ~PwHammDistAlgorithmWithPreprocessing() {
        delete(mainAlgorithm);
    }

    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t *sequences) {
        sequences = preprocessSequences(sequences);
        vector<pair<uint16_t, uint16_t>> res = mainAlgorithm->findSimilarSequences(sequences);
        postprocess(sequences);
        return res;
    }

    virtual vector<pair<uint16_t, uint16_t>>
    findSimilarSequences(uint8_t *sequences, const vector<pair<uint16_t, uint16_t>> pairs) {
        sequences = preprocessSequences(sequences);
        vector<pair<uint16_t, uint16_t>> res = mainAlgorithm->findSimilarSequences(sequences, pairs);
        postprocess(sequences);
        return res;
    }

    virtual string getName() {
        return (xParams.shuffleColumnsMode?(COLUMNS_SHUFFLER_ID + "+"):"") + mainAlgorithm->getName();
    }
};

#endif //PWHD_DATAPREPROCESSORBASE_H
