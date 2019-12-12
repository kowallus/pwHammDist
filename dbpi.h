#ifndef DBPIT_DBPI_H
#define DBPIT_DBPI_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>

using namespace std;

class DBPISolver {
protected:
    ExperimentParams &xParams;
    DBPISolver(ExperimentParams &xParams): xParams(xParams) {};

public:

    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint16_t* sequences) = 0;

};

template <bool naive>
class Brute_DBPI_Solver : public DBPISolver {

private:

    template<typename t_packed>
    vector<pair<uint16_t, uint16_t>> solve(t_packed *sequences) {
        vector<pair<uint16_t, uint16_t>> res;
        uint16_t seqInBytes = ceilDivisionBySmallInteger(xParams.m, xParams.bits_per_packed) * sizeof(t_packed);
        uint16_t seqInULLs = seqInBytes / 8;
        assert(seqInBytes == seqInULLs * 8);
        uint64_t *x = (uint64_t*) sequences;
        for(int i = 0; i < xParams.d - 1; i++) {
            uint64_t *y = x + seqInULLs;
            for (int j = i + 1; j < xParams.d; j++) {
                uint64_t dist;
                if (naive)
                    dist = bit_cost(x, y, seqInULLs);
                else
                    dist = bit_cost(x, y, seqInULLs, xParams.k);
                if (dist <= xParams.k) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
                y += seqInULLs;
            }
            x += seqInULLs;
        }
        return res;
    };

public:
    Brute_DBPI_Solver(ExperimentParams &xParams): DBPISolver(xParams) {};

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences) { return solve(sequences); };
    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint16_t* sequences){ return solve(sequences); };

};

#endif //DBPIT_DBPI_H
