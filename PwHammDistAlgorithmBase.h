#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

const string BRUTE_FORCE_ID = "bf";
const string SHORT_CIRCUIT_PREFIX_ID = "s";
const string GROUPED_PREFIX_ID = "g";
const string PIVOT_FILTER_PREFIX_ID = "p";
const string ELECTION_PIVOT_FILTER_PREFIX_ID = "P";
const string BINARY_MODE_ID_SUFFIX = "_bin";


class PwHammDistAlgorithm {
protected:
    ExperimentParams &xParams;
    PwHammDistAlgorithm(ExperimentParams &xParams): xParams(xParams) {};

public:
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t *sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
            const vector<pair<uint16_t, uint16_t>> pairs) = 0;
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
    static const uint8_t PIVOTS_COUNT_MAX = 6;
    uint8_t pivotsCount = PIVOTS_COUNT_MAX;
    uint16_t pivots[PIVOTS_COUNT_MAX];
    uint16_t* pivotDist = 0;

    // for debugging
    uint8_t* sequencesPtr;

    void preprocessing(const uint8_t *sequences) {
        if (pivot) {
            if (xParams.verbose) cout << "pivots: ";
            if (xParams.pivotsElectionMode)
                electAndCalculateDistancesToPivots(sequences);
            else
                calculateDistancesToPivots(sequences);
            if (xParams.verbose) cout << endl;
        }
    }

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
        uint16_t candidate = 0;
        if (pivotDist)
            delete[] pivotDist;
        pivotDist = new uint16_t[xParams.d * PIVOTS_COUNT_MAX];
        uint16_t* ptr = pivotDist;
        for(int p = 0; p < pivotsCount; p++) {
            if (candidate >= xParams.d) {
                pivotsCount = p;
                break;
            }
            pivots[p] = candidate;
            if (xParams.verbose) cout << candidate << "\t";
            uint8_t* x = (uint8_t*) sequences + (size_t) pivots[p] * xParams.bytesPerSequence;
            uint8_t* y = (uint8_t*) sequences;
            for(int j = 0; j < xParams.d; j++) {
                ptr[j] = findSequencesDistance<false>(x, y);
                y += xParams.bytesPerSequence;
            }
            ptr += xParams.d;
            candidate++;
        }
    }

    void electAndCalculateDistancesToPivots(const uint8_t *sequences) {
        uint16_t candidate = 0;
        if (pivotDist)
            delete[] pivotDist;
        pivotDist = new uint16_t[xParams.d * pivotsCount];
        uint16_t minPivotDist[xParams.d];
        memset(&minPivotDist, UINT8_MAX, sizeof(minPivotDist));
        const uint16_t minThreshold = xParams.k;
        uint32_t sumPivotDist[xParams.d] = { 0 };
        uint16_t* ptr = pivotDist;
        for(int p = 0; p < pivotsCount; p++) {
            pivots[p] = candidate;
            if (xParams.verbose) cout << pivots[p] << "\t";
            uint32_t maxSumDist = 0;
            uint8_t* x = (uint8_t*) sequences + (size_t) pivots[p] * xParams.bytesPerSequence;
            uint8_t* y = (uint8_t*) sequences;
            for(int j = 0; j < xParams.d; j++) {
                ptr[j] = findSequencesDistance<false>(x, y);
                y += xParams.bytesPerSequence;
                if (minPivotDist[j] > ptr[j])
                    minPivotDist[j] = ptr[j];
                if (minPivotDist[j] > minThreshold) {
                    sumPivotDist[j] += ptr[j];
                    if (maxSumDist < sumPivotDist[j]) {
                        maxSumDist = sumPivotDist[j];
                        candidate = j;
                    }
                }
            }
            if (pivots[p] == candidate) {
                pivotsCount = p + 1;
                break;
            }
            ptr += xParams.d;
        }
    }

    enum filterResult { similar, inconclusive, different };

    filterResult pivotFilter(uint16_t i, uint16_t j) {
        for(int p = 0; p < pivotsCount; p++) {
            const int iDist = pivotDist[p * xParams.d + i];
            const int jDist = pivotDist[p * xParams.d + j];
            if (abs(iDist - jDist) > xParams.k)
                return different;
            if (iDist + jDist <= xParams.k)
                return similar;
            
        }
        return inconclusive;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingStandardBrute(const uint8_t *sequences) {
        sequencesPtr = (uint8_t*) sequences;
        vector<pair<uint16_t, uint16_t>> res;
        for(int i = 0; i < xParams.d - 1; i++) {
            for (int j = i + 1; j < xParams.d; j++) {
                if (testSequencesSimilarity(sequences, i, j)) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
            }
        }
        return res;
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingGroupedApproach(const uint8_t *sequences) {
        vector<pair<uint16_t, uint16_t>> res;
        const int GROUP_SIZE = 32;
        int gEnd = (xParams.d / GROUP_SIZE) * GROUP_SIZE;
        {
            for (int gStart = 0; gStart < gEnd; gStart += GROUP_SIZE) {
                for (int i = 0; i < GROUP_SIZE; i++) {
                    for (int j = i + 1; j < GROUP_SIZE; j++) {
                        if (testSequencesSimilarity(sequences, gStart + i, gStart + j)) {
                            res.push_back(pair<uint16_t, uint16_t>(gStart + i, gStart + j));
                        }
                    }
                }
            }
        }
        {
            for(int g1Start = 0; g1Start < gEnd; g1Start += GROUP_SIZE) {
                for (int g2Start = g1Start + GROUP_SIZE; g2Start < gEnd; g2Start += GROUP_SIZE) {
                    for (int i = 0; i < GROUP_SIZE; i++) {
                        for (int j = 0; j < GROUP_SIZE; j++) {
                            if (testSequencesSimilarity(sequences, g1Start + i, g2Start + j)) {
                                res.push_back(pair<uint16_t, uint16_t>(g1Start + i, g2Start + j));
                            }
                        }
                    }
                }
            }
        }
        {
            for(int i = 0; i < xParams.d - 1; i++) {
                int j = i < gEnd?gEnd:i + 1;
                for (; j < xParams.d; j++) {
                    if (testSequencesSimilarity(sequences, i, j)) {
                        res.push_back(pair<uint16_t, uint16_t>(i, j));
                    }
                }
            }
        }
        return res;
    };

public:
    BrutePwHammDistAlgorithm(ExperimentParams &xParams): PwHammDistAlgorithm(xParams), seqInULLs(xParams.bytesPerSequence / 8) {
        if (xParams.bytesPerSequence % xParams.ALINGMENT_IN_BYTES != 0) {
            fprintf(stderr, "ERROR: brute algorithm does not support unaligned data.\n");
            exit(EXIT_FAILURE);
        }
    };

    virtual ~BrutePwHammDistAlgorithm() {
        if (pivotDist)
            delete[] pivotDist;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        preprocessing(sequences);
        if (grouped)
            return findSimilarSequencesUsingGroupedApproach(sequences);
        else
            return findSimilarSequencesUsingStandardBrute(sequences);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        preprocessing(sequences);
        vector<pair<uint16_t, uint16_t>> res;
        for (pair<uint16_t, uint16_t> pair: pairs) {
            if (testSequencesSimilarity(sequences, pair.first, pair.second)) {
                res.push_back(pair);
            }
        }
        return res;
    }

    inline bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) {
        if(pivot) {
            filterResult res = pivotFilter(i, j);
            if (res == different)
                return false;
            else if (res == similar)
                return true;
        }
        uint8_t* x = (uint8_t*) sequences + (size_t) i * xParams.bytesPerSequence;
        uint8_t* y = (uint8_t*) sequences + (size_t) j * xParams.bytesPerSequence;
        return testSequencesSimilarity(x, y);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        return findSequencesDistance<true>(seq1, seq2) <= xParams.k;
    }

    string getName() {
        return (pivot?(xParams.pivotsElectionMode?ELECTION_PIVOT_FILTER_PREFIX_ID:PIVOT_FILTER_PREFIX_ID):"") + (grouped?GROUPED_PREFIX_ID:"") +
            (shortcircuit?SHORT_CIRCUIT_PREFIX_ID:"") + BRUTE_FORCE_ID +
            (binaryAlphabet?BINARY_MODE_ID_SUFFIX:"");
    }
};

#endif //ALGORITHM_BASE_H
