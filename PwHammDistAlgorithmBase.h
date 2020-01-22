#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

const string NAIVE_BRUTE_FORCE_ID = "nbf";
const string SHORT_CIRCUIT_BRUTE_FORCE_ID = "sbf";
const string GROUPED_PREFIX_ID = "g";
const string GROUPED_NAIVE_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + "nbf";
const string GROUPED_SHORT_CIRCUIT_BRUTE_FORCE_ID = GROUPED_PREFIX_ID + "sbf";
const string BINARY_MODE_ID_SUFFIX = "_bin";

class PwHammDistAlgorithm {
protected:
    ExperimentParams &xParams;
    PwHammDistAlgorithm(ExperimentParams &xParams): xParams(xParams) {};

public:
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t *sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
            const vector<pair<uint16_t, uint16_t>> pairs) = 0;
    virtual bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) = 0;
    virtual bool testSequencesSimilarity(const void* seq1, const void* seq2) = 0;
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

template <bool naive, bool binaryAlphabet, typename uint = uint8_t, bool grouped = false>
class BrutePwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    const uint16_t seqInULLs;

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingStandardBrute(const uint8_t *sequences) {
        vector<pair<uint16_t, uint16_t>> res;
        uint64_t *x = (uint64_t*) sequences;
        for(int i = 0; i < xParams.d - 1; i++) {
            uint64_t *y = x + seqInULLs;
            for (int j = i + 1; j < xParams.d; j++) {
                if (testSequencesSimilarity(x, y)) {
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
                        if (testSequencesSimilarity(x, y)) {
                            res.push_back(pair<uint16_t, uint16_t>(i, j));
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
                            if (testSequencesSimilarity(x, y)) {
                                res.push_back(pair<uint16_t, uint16_t>(i, j));
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
                    if (testSequencesSimilarity(x, y)) {
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

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        if (grouped)
            return findSimilarSequencesUsingGroupedApproach(sequences);
        else
            return findSimilarSequencesUsingStandardBrute(sequences);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                const vector<pair<uint16_t, uint16_t>> pairs) {
        vector<pair<uint16_t, uint16_t>> res;
        for (pair<uint16_t, uint16_t> pair: pairs) {
            if (testSequencesSimilarity(sequences, pair.first, pair.second)) {
                res.push_back(pair);
            }
        }
        return res;
    }

    inline bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) {
        uint8_t* x = (uint8_t*) sequences + (size_t) i * xParams.bytesPerSequence;
        uint8_t* y = (uint8_t*) sequences + (size_t) j * xParams.bytesPerSequence;
        return testSequencesSimilarity(x, y);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        uint16_t dist;
        if (binaryAlphabet)
            dist = naive? hammingDistanceBinary((uint64_t*) seq1, (uint64_t*)seq2, seqInULLs):
                    hammingDistanceBinary((uint64_t*) seq1, (uint64_t*) seq2, seqInULLs, xParams.k);
        else
            dist = naive?hammingDistance((uint*) seq1, (uint*) seq2, xParams.m):
                   hammingDistance((uint*) seq1, (uint*) seq2, xParams.m, xParams.k);
        return dist <= xParams.k;
    }

    string getName() { return (grouped?GROUPED_PREFIX_ID:"") + (naive?NAIVE_BRUTE_FORCE_ID:SHORT_CIRCUIT_BRUTE_FORCE_ID) + (binaryAlphabet?BINARY_MODE_ID_SUFFIX:""); }
};

#endif //ALGORITHM_BASE_H
