#ifndef PWHD_HASHINGPWHAMMDISTALGORITHM_H
#define PWHD_HASHINGPWHAMMDISTALGORITHM_H

#include "PwHammDistAlgorithmBase.h"
#include "hashes/Hashes.h"

const string HASHING_BASED_FILTER_PWHD_ID = "hbf";

template<typename uint>
class HashingBasedPwHammDistAlgorithm : public PwHammDistAlgorithm {
private:

    const static uint8_t HASH_BITS = 11;
    const static uint32_t HASH_SLOTS = ((uint32_t) 1) << (HASH_BITS - 1);
    const static uint32_t HASH_MASK = HASH_SLOTS - 1;
    int htSize;

    uint16_t h;
    int hInBytes;

    void preprocessing(const uint8_t *sequences) {
        htSize = xParams.d * HASH_SLOTS;
        h = xParams.m / (xParams.k + xParams.hbf_L);
        hInBytes = h * sizeof(uint);
    };

    void postProcessing() {
    }

    uint32_t hashFunc(uint8_t *x) {
        return XXH32((const void*)x, hInBytes, (uint32_t) 9876543210UL) & HASH_MASK;
    }

public:
    HashingBasedPwHammDistAlgorithm(ExperimentParams &xParams) : PwHammDistAlgorithm(
            xParams) {}

    virtual ~HashingBasedPwHammDistAlgorithm() { }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        preprocessing(sequences);
        vector<uint16_t> hT(htSize);
        vector<uint16_t> hTCounts(HASH_SLOTS);
        vector<uint32_t> hTKeys;
        vector<uint16_t> pairsCounter(xParams.d * xParams.d, 0);
        vector<pair<uint16_t, uint16_t>> res;
        for(int i = 0; i < xParams.m; i += h) {
            memset(&hTCounts[0], 0, HASH_SLOTS * sizeof(uint16_t));
            hTKeys.clear();
            uint8_t* x = (uint8_t*) sequences + (size_t) i * sizeof(uint);
            for(int j = 0; j < xParams.d; j++) {
                uint32_t hash = hashFunc(x);
                hT[xParams.d * hash + hTCounts[hash]++] = j;
                if (hTCounts[hash] == 2)
                    hTKeys.push_back(hash);
                x += xParams.bytesPerSequence;
            }
            for(uint32_t k = 0; k < hTKeys.size(); k++) {
                uint32_t key = hTKeys[k];
                uint16_t* hTSlot = &hT[xParams.d * key];
                for (int u = 0; u < hTCounts[key] - 1; u++) {
                    for (int v = u + 1; v < hTCounts[key]; v++) {
                        uint16_t i1 = hTSlot[u];
                        uint16_t i2 = hTSlot[v];
                        if (++pairsCounter[i1 * xParams.d + i2] == xParams.hbf_L)
                            res.push_back(pair<uint16_t, uint16_t>(i1, i2));
                    }
                }
            }
        }
        postProcessing();
        return res;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        fprintf(stderr, "ERROR: hashing-based filter supported only for all sequences.\n");
        exit(EXIT_FAILURE);
    }

    string getName() {
        return HASHING_BASED_FILTER_PWHD_ID;
    }
};


#endif //PWHD_HASHINGPWHAMMDISTALGORITHM_H
