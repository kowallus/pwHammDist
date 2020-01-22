#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

const string NAIVE_BRUTE_FORCE_ID = "bf";
const string SHORT_CIRCUIT_BRUTE_FORCE_ID = "sbf";
const string GROUPED_PREFIX_ID = "g";
const string GROUPED_NAIVE_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + NAIVE_BRUTE_FORCE_ID;
const string GROUPED_SHORT_CIRCUIT_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + SHORT_CIRCUIT_BRUTE_FORCE_ID;
const string BINARY_MODE_ID_SUFFIX = "_bin";
const string PIVOT_FILTER_PREFIX_ID = "p";
const string PIVOT_FILTER_ID = PIVOT_FILTER_PREFIX_ID + NAIVE_BRUTE_FORCE_ID;
const string SHORT_CIRCUIT_PIVOT_FILTER_ID = PIVOT_FILTER_PREFIX_ID + SHORT_CIRCUIT_BRUTE_FORCE_ID;
const string GROUPED_PIVOT_FILTER_ID = PIVOT_FILTER_PREFIX_ID + GROUPED_NAIVE_BRUTE_FORCE_ID;
const string GROUPED_SHORT_CIRCUIT_PIVOT_FILTER_ID = PIVOT_FILTER_PREFIX_ID + GROUPED_SHORT_CIRCUIT_BRUTE_FORCE_ID;


class PwHammDistAlgorithm {
protected:
    ExperimentParams &xParams;
    PwHammDistAlgorithm(ExperimentParams &xParams): xParams(xParams) {};

public:
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t *sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
            const vector<pair<uint16_t, uint16_t>> pairs) = 0;
    virtual bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) = 0;
    virtual string getName() = 0;
};

class PwHammDistAlgorithmFactory {
private:
    string algorithmName;

public:
    PwHammDistAlgorithmFactory(string algorithmName): algorithmName(algorithmName) { }
    string getAlgorithmName() { return algorithmName; }
    virtual PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) = 0;
};

template <bool shortcircuit, bool binaryAlphabet, typename uint = uint8_t, bool grouped = false, bool pivot = false>
class BrutePwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    const uint16_t seqInULLs;
    static const uint8_t PIVOTS_COUNT = 8;
    uint16_t pivots[PIVOTS_COUNT];
    uint16_t* pivotDist = 0;

    // for debugging
    uint8_t* sequencesPtr;

    template<bool allowShortCircuit>
    uint16_t findSequencesDistance(const void* seq1, const void* seq2) {
        uint16_t dist;
        if (binaryAlphabet)
            dist = (allowShortCircuit && shortcircuit)?hammingDistanceBinary((uint64_t*) seq1, (uint64_t*) seq2, seqInULLs, xParams.k):
                   hammingDistanceBinary((uint64_t*) seq1, (uint64_t*)seq2, seqInULLs);
        else
            dist = (allowShortCircuit && shortcircuit)?hammingDistance((uint*) seq1, (uint*) seq2, xParams.m, xParams.k):
                   hammingDistance((uint*) seq1, (uint*) seq2, xParams.m);
        return dist;
    }

    void calculateDistancesToPivots(const uint8_t *sequences) {
        for(int p = 0; p < PIVOTS_COUNT; p++) {
            //pivots[p] = ((int) xParams.d) * p / PIVOTS_COUNT;
            pivots[p] = p;
            cout << pivots[p] << "\t";
        }
        cout << endl;
        if (pivotDist)
            delete[] pivotDist;
        pivotDist = new uint16_t[xParams.d * PIVOTS_COUNT];
        uint16_t* ptr = pivotDist;
        for(int p = 0; p < PIVOTS_COUNT; p++) {
            if (pivots[p] >= xParams.d) {
                fprintf(stderr, "ERROR: invalid sequence %u pivot (max: %u).\n", pivots[p], xParams.d - 1);
                exit(EXIT_FAILURE);
            }
            uint8_t* x = (uint8_t*) sequences + (size_t) pivots[p] * xParams.bytesPerSequence;
            uint8_t* y = (uint8_t*) sequences;
            for(int j = 0; j < xParams.d; j++) {
                ptr[j] = findSequencesDistance<false>(x, y);
                y += xParams.bytesPerSequence;
            }
            ptr += xParams.d;
        }
    }

    bool pivotFilter(uint16_t i, uint16_t j) {
        for(int p = 0; p < PIVOTS_COUNT; p++) {
            if (abs((int) pivotDist[p * xParams.d + i] - (int) pivotDist[p * xParams.d + j]) > xParams.k)
                return true;
        }
        return false;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingStandardBrute(const uint8_t *sequences) {
        sequencesPtr = (uint8_t*) sequences;
        vector<pair<uint16_t, uint16_t>> res;
        uint64_t *x = (uint64_t*) sequences;
        for(int i = 0; i < xParams.d - 1; i++) {
            uint64_t *y = x + seqInULLs;
            for (int j = i + 1; j < xParams.d; j++) {
                if (!(pivot && pivotFilter(i, j)) && testSequencesSimilarity(x, y)) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
                y += seqInULLs;
            }
            x += seqInULLs;
        }
        return res;
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingGroupedApproach(const uint8_t *sequences) {
        vector<pair<uint16_t, uint16_t>> res;
        const int GROUP_SIZE = 32;
        int gEnd = (xParams.d / GROUP_SIZE) * GROUP_SIZE;
        {
            uint8_t *g = (uint8_t *) sequences;
            for (int gStart = 0; gStart < gEnd; gStart += GROUP_SIZE) {
                uint8_t *x = g;
                for (int i = 0; i < GROUP_SIZE; i++) {
                    uint8_t *y = x + xParams.bytesPerSequence;
                    for (int j = i + 1; j < GROUP_SIZE; j++) {
                        if (!(pivot && pivotFilter(gStart + i, gStart + j)) && testSequencesSimilarity(x, y)) {
                            res.push_back(pair<uint16_t, uint16_t>(gStart + i, gStart + j));
                        }
                        y += xParams.bytesPerSequence;
                    }
                    x += xParams.bytesPerSequence;
                }
                g += xParams.bytesPerSequence * GROUP_SIZE;
            }
        }
        {
            uint8_t *g1 = (uint8_t*) sequences;
            for(int g1Start = 0; g1Start < gEnd; g1Start += GROUP_SIZE) {
                uint8_t *g2 = g1 + xParams.bytesPerSequence * GROUP_SIZE;
                for (int g2Start = g1Start + GROUP_SIZE; g2Start < gEnd; g2Start += GROUP_SIZE) {
                    uint8_t *x = g1;
                    for (int i = 0; i < GROUP_SIZE; i++) {
                        uint8_t *y = g2;
                        for (int j = 0; j < GROUP_SIZE; j++) {
                            if (!(pivot && pivotFilter(g1Start + i, g2Start + j)) && testSequencesSimilarity(x, y)) {
                                res.push_back(pair<uint16_t, uint16_t>(g1Start + i, g2Start + j));
                            }
                            y += xParams.bytesPerSequence;
                        }
                        x += xParams.bytesPerSequence;
                    }
                    g2 += xParams.bytesPerSequence * GROUP_SIZE;
                }
                g1 += xParams.bytesPerSequence * GROUP_SIZE;
            }
        }
        {
            uint64_t *x = (uint64_t*) sequences;
            for(int i = 0; i < xParams.d - 1; i++) {
                int j = i + 1;
                uint64_t *y = x + seqInULLs;
                if ((uint8_t*) y < sequences + gEnd * xParams.bytesPerSequence) {
                    j = gEnd;
                    y = (uint64_t*) (sequences + gEnd * xParams.bytesPerSequence);
                }
                for (; j < xParams.d; j++) {
                    if (!(pivot && pivotFilter(i, j)) && testSequencesSimilarity(x, y)) {
                        res.push_back(pair<uint16_t, uint16_t>(i, j));
                    }
                    y += seqInULLs;
                }
                x += seqInULLs;
            }
        }
        return res;
    };

public:
    BrutePwHammDistAlgorithm(ExperimentParams &xParams): PwHammDistAlgorithm(xParams), seqInULLs(xParams.bytesPerSequence / 8) {
        if (xParams.bytesPerSequence != seqInULLs * 8) {
            fprintf(stderr, "ERROR: brute algorithm does not support unaligned data.\n");
            exit(EXIT_FAILURE);
        }
    };

    virtual ~BrutePwHammDistAlgorithm() {
        if (pivotDist)
            delete[] pivotDist;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        if (pivot)
            calculateDistancesToPivots(sequences);
        if (grouped)
            return findSimilarSequencesUsingGroupedApproach(sequences);
        else
            return findSimilarSequencesUsingStandardBrute(sequences);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        if (pivot)
            calculateDistancesToPivots(sequences);
        vector<pair<uint16_t, uint16_t>> res;
        for (pair<uint16_t, uint16_t> pair: pairs) {
            if (testSequencesSimilarity(sequences, pair.first, pair.second)) {
                res.push_back(pair);
            }
        }
        return res;
    }

    inline bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) {
        if(pivot && pivotFilter(i, j))
            return false;
        uint8_t* x = (uint8_t*) sequences + (size_t) i * xParams.bytesPerSequence;
        uint8_t* y = (uint8_t*) sequences + (size_t) j * xParams.bytesPerSequence;
        return testSequencesSimilarity(x, y);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        return findSequencesDistance<true>(seq1, seq2) <= xParams.k;
    }

    string getName() { return (pivot?PIVOT_FILTER_PREFIX_ID:"") + (grouped?GROUPED_PREFIX_ID:"") + (shortcircuit?SHORT_CIRCUIT_BRUTE_FORCE_ID:NAIVE_BRUTE_FORCE_ID) + (binaryAlphabet?BINARY_MODE_ID_SUFFIX:""); }
};

#endif //ALGORITHM_BASE_H
