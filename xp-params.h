#ifndef DBPIT_XP_PARAMS_H
#define DBPIT_XP_PARAMS_H

#include <cstdint>

struct ExperimentParams {

    // basic params

    uint16_t m = 512; // sequence length
    uint16_t d = 1024; // no of sequences
    uint16_t k = 32; // max hamming length

    // extra params

    uint8_t bits_per_packed = 12; //bits in packed value
    uint16_t ones_in_promiles = 900; // density of 1 bits in sequences (in promiles)

    const uint8_t DISABLE_BITS_PER_PACKED = UINT8_MAX;
    void disableOnesInPromiles() { ones_in_promiles = DISABLE_ONES_IN_PROMILES; }
    bool isDisabledOnesInPromiles() { return ones_in_promiles == DISABLE_ONES_IN_PROMILES; };

    const uint16_t DISABLE_ONES_IN_PROMILES = UINT16_MAX;
    void disableBitsPerPacked() { bits_per_packed = DISABLE_BITS_PER_PACKED; }
    bool isDisabledBitsPerPacked() { return bits_per_packed == DISABLE_BITS_PER_PACKED; }

};


#endif //DBPIT_XP_PARAMS_H
