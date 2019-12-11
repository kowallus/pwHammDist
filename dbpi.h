#ifndef DBPIT_DBPI_H
#define DBPIT_DBPI_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>

using namespace std;

class DBPISolver {
protected:
    ExperimentParams &xParams;
    DBPISolver(ExperimentParams &xParams): xParams(xParams) {};

public:

    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint16_t* sequences) = 0;

};

class Brute_DBPI_Solver : public DBPISolver {

private:

    template<typename t_packed>
    vector<pair<uint16_t, uint16_t>> solve(t_packed *sequences) { return vector<pair<uint16_t, uint16_t>>(); };

public:
    Brute_DBPI_Solver(ExperimentParams &xParams): DBPISolver(xParams) {};

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint8_t* sequences) { return solve(sequences); };
    vector<pair<uint16_t, uint16_t>> findSimilarSequences(uint16_t* sequences){ return solve(sequences); };

};

#endif //DBPIT_DBPI_H
