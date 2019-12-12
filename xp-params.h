#ifndef DBPIT_XP_PARAMS_H
#define DBPIT_XP_PARAMS_H

#include <cstdint>
#include <string>

struct ExperimentParams {

    // basic params

    uint16_t m = 512; // sequence length
    uint16_t d = 1024; // no of sequences
    uint16_t k = 32; // max hamming length
    std::string solverID = "nbf";

    // extra params

    uint8_t bitsPerPacked = 8; //bits in packed value
    uint16_t onesInPromiles = 900; // density of 1 bits in sequences (in promiles)

    // data properties

    uint32_t bytesPerSequence = 0;

    const uint8_t DISABLE_BITS_PER_PACKED = UINT8_MAX;
    void disableOnesInPromiles() { onesInPromiles = DISABLE_ONES_IN_PROMILES; }
    bool isDisabledOnesInPromiles() { return onesInPromiles == DISABLE_ONES_IN_PROMILES; };

    const uint16_t DISABLE_ONES_IN_PROMILES = UINT16_MAX;
    void disableBitsPerPacked() { bitsPerPacked = DISABLE_BITS_PER_PACKED; }
    bool isDisabledBitsPerPacked() { return bitsPerPacked == DISABLE_BITS_PER_PACKED; }

};

#endif //DBPIT_XP_PARAMS_H
