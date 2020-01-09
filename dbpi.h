#ifndef DBPIT_DBPI_H
#define DBPIT_DBPI_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

const string NAIVE_BRUTE_FORCE_ID = "nbf";
const string SHORT_CIRCUIT_BRUTE_FORCE_ID = "sbf";
const string BINARY_MODE_ID_SUFFIX = "_bin";

class DBPISolver {
protected:
    ExperimentParams &xParams;
    DBPISolver(ExperimentParams &xParams): xParams(xParams) {};

public:
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t *sequences) = 0;
    virtual string getName() = 0;
};

template <bool naive, bool binaryAlphabet>
class Brute_DBPI_Solver : public DBPISolver {
public:
    Brute_DBPI_Solver(ExperimentParams &xParams): DBPISolver(xParams) {};

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        vector<pair<uint16_t, uint16_t>> res;
        uint16_t seqInULLs = xParams.bytesPerSequence / 8;
        assert(xParams.bytesPerSequence == seqInULLs * 8);
        uint64_t *x = (uint64_t*) sequences;
        for(int i = 0; i < xParams.d - 1; i++) {
            uint64_t *y = x + seqInULLs;
            for (int j = i + 1; j < xParams.d; j++) {
                uint64_t dist;
                if (binaryAlphabet)
                    dist = naive? hammingDistanceBinary(x, y, seqInULLs): hammingDistanceBinary(x, y, seqInULLs,
                                                                                                xParams.k);
                else
                    dist = naive?hammingDistance((uint8_t*) x, (uint8_t*) y, xParams.bytesPerSequence):
                           hammingDistance((uint8_t*) x, (uint8_t*) y, xParams.bytesPerSequence, xParams.k);
                if (dist <= xParams.k) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
                y += seqInULLs;
            }
            x += seqInULLs;
        }
        return res;
    };

    string getName() { return (naive?NAIVE_BRUTE_FORCE_ID:SHORT_CIRCUIT_BRUTE_FORCE_ID) + (binaryAlphabet?BINARY_MODE_ID_SUFFIX:""); };

};

#endif //DBPIT_DBPI_H
