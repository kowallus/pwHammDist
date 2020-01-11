#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

const string NAIVE_BRUTE_FORCE_ID = "nbf";
const string SHORT_CIRCUIT_BRUTE_FORCE_ID = "sbf";
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

template <bool naive, bool binaryAlphabet, typename uint = uint8_t>
class BrutePwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    const uint16_t seqInULLs;

public:
    BrutePwHammDistAlgorithm(ExperimentParams &xParams): PwHammDistAlgorithm(xParams), seqInULLs(xParams.bytesPerSequence / 8) {
        if (xParams.bytesPerSequence != seqInULLs * 8) {
            fprintf(stderr, "ERROR: brute algorithm does not support unaligned data.\n");
            exit(EXIT_FAILURE);
        }
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
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

    string getName() { return (naive?NAIVE_BRUTE_FORCE_ID:SHORT_CIRCUIT_BRUTE_FORCE_ID) + (binaryAlphabet?BINARY_MODE_ID_SUFFIX:""); }
};

#endif //ALGORITHM_BASE_H
