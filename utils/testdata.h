#ifndef TESTDATA_H
#define TESTDATA_H

#include <vector>
#include <random>
#include "helper.h"
#include "../xp-params.h"

using namespace std;

std::mt19937 randgenerator;

void getRandomValues(uint8_t* data, ExperimentParams& xParams) {
    randgenerator.seed(randgenerator.default_seed);
    uint16_t m = xParams.m;
    uint16_t d = xParams.d;
    uint8_t bitsPerPacked = xParams.bitsPerPacked;
    uint16_t bytesPerSeq = xParams.bytesPerSequence;
    uint16_t onesInPromiles = xParams.onesInPromiles;
    uint8_t* cur;
    for(uint16_t i = 0; i < d; i++) {
        cur = (data + i * bytesPerSeq) - 1;
        uint8_t b = 0;
        for (uint16_t j = 0; j < m; j++) {
            if (b % 8 == 0)
                *(++cur) = 0;
            if (randgenerator() % 1000 < onesInPromiles)
                *cur += 1 << (b % 8);
            if (++b == bitsPerPacked)
                b = 0;
        }
    }
}

#endif /* TESTDATA_H */

