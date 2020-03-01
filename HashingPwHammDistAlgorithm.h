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
        if (!h) {
            fprintf(stderr, "ERROR: sequence to short for hashing (elements = %d; k = %d; L = %d).\n", xParams.m, xParams.k, xParams.hbf_L);
            exit(EXIT_FAILURE);
        }
        hInBytes = h * sizeof(uint);
        if (xParams.verbose) {
            cout << "Hash bits: " << (int) HASH_BITS << "; ";
            cout << "L: " << xParams.hbf_L << "; ";
            cout << "h: " << h << endl;
        }
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

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences) {
        preprocessing(sequences);
        vector<vector<uint16_t>> hT(HASH_SLOTS);
        vector<uint32_t> hTKeys;
        vector<uint16_t> pairsCounter(xParams.d * xParams.d, 0);
        vector<pair<uint16_t, uint16_t>> res;
        if (xParams.verbose)
            cout << "hT initialization... " << " (" << time_millis() << " msec)" << endl;
        for(int i = 0; i < xParams.m; i += h) {
            for(int s = 0; s < HASH_SLOTS; s++)
                hT[s].clear();
            hTKeys.clear();
            uint8_t* x = (uint8_t*) sequences + (size_t) i * sizeof(uint);
            for(int j = 0; j < xParams.d; j++) {
                uint32_t hash = hashFunc(x);
                hT[hash].push_back(j);
                if (hT[hash].size() == 2)
                    hTKeys.push_back(hash);
                x += xParams.bytesPerSequence;
            }
            int tmp = 0;
            for(uint32_t k = 0; k < hTKeys.size(); k++) {
                uint32_t key = hTKeys[k];
                vector<uint16_t>& hTSlot = hT[key];
                for (int u = 0; u < hTSlot.size() - 1; u++) {
                    for (int v = u + 1; v < hTSlot.size(); v++) {
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

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        fprintf(stderr, "ERROR: hashing-based filter supported only for all sequences.\n");
        exit(EXIT_FAILURE);
    }

    string getName() {
        return HASHING_BASED_FILTER_PWHD_ID;
    }
};


#endif //PWHD_HASHINGPWHAMMDISTALGORITHM_H
