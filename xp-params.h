#ifndef PWHD_XP_PARAMS_H
#define PWHD_XP_PARAMS_H

#include <cstdint>
#include <string>

struct ExperimentParams {

    // basic params

    bool datasetMode = false;
    std::string datasetFileName;
    bool binaryMode = false;

    bool alignSequencesTo256bits = true;
    static const int ALINGMENT_IN_BYTES = 32;

    uint16_t m = 0; // sequence length
    uint16_t d = 0; // no of sequences
    uint16_t k = 0; // max hamming length

    std::string algorithmID = "bf";
    bool shortCircuitMode = true;
    bool compactMode = false;
    bool groupedBruteMode = false;
    bool pivotsFilterMode = false;
    bool pivotsElectionMode = false;

    // extra params

    bool verbose = true;

    static const uint8_t DISABLE_BITS_PER_PACKED = UINT8_MAX;
    static const uint16_t DISABLE_ONES_IN_PROMILES = UINT16_MAX;

    uint8_t bitsPerPacked = DISABLE_BITS_PER_PACKED; //bits in packed value
    uint16_t onesInPromiles = DISABLE_ONES_IN_PROMILES; // density of 1 bits in sequences (in promiles)

    // data properties

    bool dnaDataMode = false;
    uint8_t symbol2value[UINT8_MAX] = { 0 };
    uint32_t alphabetSize = 0;
    uint8_t bytesPerElement = 0;
    uint32_t bytesPerSequence = 0;

    void enableOnesInPromiles() { onesInPromiles = 500; }
    bool isOnesInPromilesEnabled() const { return onesInPromiles != DISABLE_ONES_IN_PROMILES; };

    void enableBitsPerPacked() { bitsPerPacked = 8; }
    bool isBitsPerPackedEnabled() const { return bitsPerPacked != DISABLE_BITS_PER_PACKED; };

    bool isInBinaryMode() const { return binaryMode; }
    bool isInDatasetMode() const { return datasetMode; }

    void enableBinaryMode() {
        binaryMode = true;
        compactMode = false;
        alphabetSize = 2;
        enableOnesInPromiles();
        enableBitsPerPacked();
    }

    void enableDatasetMode() {
        datasetMode = true;
    }
};

#endif //PWHD_XP_PARAMS_H
