#ifndef DBPIT_XP_PARAMS_H
#define DBPIT_XP_PARAMS_H

#include <cstdint>
#include <string>

struct ExperimentParams {

    // basic params

    bool datasetMode = false;
    std::string datasetFileName;
    bool binaryMode = false;

    bool alignSequences = true;

    uint16_t m = 512; // sequence length
    uint16_t d = 1024; // no of sequences
    uint16_t k = 32; // max hamming length
    std::string solverID = "nbf";

    // extra params

    const uint8_t DISABLE_BITS_PER_PACKED = UINT8_MAX;
    const uint16_t DISABLE_ONES_IN_PROMILES = UINT16_MAX;

    uint8_t bitsPerPacked = 8; //bits in packed value
    uint16_t onesInPromiles = 900; // density of 1 bits in sequences (in promiles)

    // data properties

    uint32_t bytesPerSequence = 0;


    void enableOnesInPromiles() { onesInPromiles = 900; }
    bool isOnesInPromilesEnabled() const { return onesInPromiles != DISABLE_ONES_IN_PROMILES; };

    void enableBitsPerPacked() { bitsPerPacked = 8; }
    bool isBitsPerPackedEnabled() const { return bitsPerPacked != DISABLE_BITS_PER_PACKED; };

    bool isInBinaryMode() const { return binaryMode; }
    bool isInDatasetMode() const { return datasetMode; }

    void enableBinaryMode() {
        binaryMode = true;
        enableOnesInPromiles();
        enableBitsPerPacked();
    }

    void enableDatasetMode() {
        datasetMode = true;
    }

};

#endif //DBPIT_XP_PARAMS_H
