#ifndef TESTDATA_H
#define TESTDATA_H

#include <vector>
#include <random>

using namespace std;

std::mt19937 randgenerator;

template<typename t_packed>
void getRandomValues(t_packed* data, uint16_t m, uint16_t d, uint8_t bits_packed, uint16_t ones_in_promiles) {
    randgenerator.seed(randgenerator.default_seed);
    uint16_t noOfPacked = ((m - 1) / bits_packed) + 1;
    t_packed* cur = data;
    for(uint16_t i = 0; i < d; i++) {
        for (uint16_t j = 0; j < noOfPacked; j++) {
            *cur = 0;
            for (uint8_t b = 0; b < bits_packed; b++) {
                if (randgenerator() % 1000 < ones_in_promiles)
                    *cur += 1 << b;
            }
            cur++;
        }
    }
}

#endif /* TESTDATA_H */

