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

    uint32_t totalPairsCount() { return (uint32_t) d * (d - 1) / 2; }

    std::string algorithmID = "bf";
    bool shortCircuitMode = true;
    bool compactMode = false;
    bool groupedBruteMode = false;
    bool pivotsFilterMode = false;
    bool pivotsRandomization = false;
    bool interleaveBitsMode = false;
    bool lazyInterleaveBitsMode = false;
    bool statsBasedQuantization = true;
    bool shuffleColumnsMode = false;
    bool perfectHashing = false; // INCOMPLETE

    void setModelAlgorithm() {
        algorithmID = "bf";
        shortCircuitMode = true;
        compactMode = false;
        groupedBruteMode = false;
        pivotsFilterMode = false;
        pivotsRandomization = false;
        interleaveBitsMode = false;
        lazyInterleaveBitsMode = false;
        statsBasedQuantization = true;
        shuffleColumnsMode = false;
        perfectHashing = false;
    }

    // experiment stats

    uint32_t pairsFound = 0; //no of similar pairs found
    uint32_t pairsVerified = 0; //no of pairs processed with sequences access
    uint32_t pairsFilterRejected = 0; //no of pairs individually rejected by filter (without sequences accessed)
    uint32_t pairsFilterAccepted = 0; //no of pairs accepted by filter (without sequences accessed)
    uint32_t pairsIgnored = 0; //no of pairs rejected without any access

    uint32_t preStageTimeInUsec = 0; //time of preprocessing (together with preFilter algorithm)

    inline void resetStats() {
#ifdef XP_STATS
        pairsFound = 0;
        pairsVerified = 0;
        pairsFilterRejected = 0;
        pairsFilterAccepted = 0;
        pairsIgnored = 0;
        preStageTimeInUsec = 0;
#endif
    }

    inline void statsIncPairsVerification() {
#ifdef XP_STATS
        pairsVerified++;
#endif
    }

    inline void statsIncPairsFilterAccepted(int acceptedCount = 1) {
#ifdef XP_STATS
        pairsFilterAccepted += acceptedCount;
#endif
    }

    inline void statsIncPairsFilterRejected(int rejectedCount = 1) {
#ifdef XP_STATS
        pairsFilterRejected += rejectedCount;
#endif
    }

    inline void statsIncPairsIgnored(int ignoredCount = 1) {
#ifdef XP_STATS
        pairsIgnored += ignoredCount;
#endif
    }

    // hashing-filter params

    uint16_t hbf_L = 10;

    // quantization params

    uint8_t quantizationBits = 3;

    // perfect-hashing

    uint16_t hashBlockLengthInULLs = 16;

    // extra params

    bool verbose = true;

    static const uint8_t DISABLE_BITS_PER_PACKED = UINT8_MAX;
    static const uint16_t DISABLE_ONES_IN_PROMILES = UINT16_MAX;

    uint8_t bitsPerPacked = DISABLE_BITS_PER_PACKED; //bits in packed value
    uint16_t onesInPromiles = DISABLE_ONES_IN_PROMILES; // density of 1 bits in sequences (in promiles)

    // data properties

    bool dnaDataMode = false;
    uint8_t symbol2value[UINT8_MAX] = { 0 };
    uint32_t alphabetSizeUpperBound = 0;
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
        interleaveBitsMode = false;
        alphabetSizeUpperBound = 2;
        enableOnesInPromiles();
        enableBitsPerPacked();
    }

    void enableDatasetMode() {
        datasetMode = true;
    }

    void disableFiltration() {
        pivotsFilterMode = false;
    }


};

#endif //PWHD_XP_PARAMS_H
