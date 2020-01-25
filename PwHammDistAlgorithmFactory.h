#ifndef PWHD_ALGORITHMFACTORY_H
#define PWHD_ALGORITHMFACTORY_H

#include "QuantizationPwHammDistAlgorithm.h"
#include "xp-params.h"
#include <map>

extern map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap;

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams);

class BrutePwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
private:
    static const int SHORTCIRCUIT_DISPATCH_BIT = 0;
    static const int BINARY_DISPATCH_BIT = 1;
    static const int GROUPED_DISPATCH_BIT = 2;
    static const int PIVOT_DISPATCH_BIT = 3;
    static const int NIBBLE_DISPATCH_BIT = 4;
    static const int DISPATCH_ENCODED_PARAMS_MAX = 32;

    template <bool shortcircuit, bool binaryAlphabet, typename uint, bool grouped, bool pivot, bool nibble>
    static PwHammDistAlgorithm* getBrutePwHammAlgorithmInstanceTemplate(ExperimentParams &xParams) {
        return new BrutePwHammDistAlgorithm<shortcircuit, binaryAlphabet, uint, grouped, pivot, nibble>(xParams);
    }

    template<typename uint, std::size_t...Is>
    PwHammDistAlgorithm* runtimeDispatchAlgorithmInstance(ExperimentParams &xParams, int encodedParams, std::index_sequence<Is...>) {
        using f_t = PwHammDistAlgorithm*(ExperimentParams &xParams);
        f_t* fs[] = {&getBrutePwHammAlgorithmInstanceTemplate<
                     (Is >> SHORTCIRCUIT_DISPATCH_BIT) & 1,
                     (Is >> BINARY_DISPATCH_BIT) & 1,
                     uint,
                     (Is >> GROUPED_DISPATCH_BIT) & 1,
                     (Is >> PIVOT_DISPATCH_BIT) & 1,
                     (Is >> NIBBLE_DISPATCH_BIT) & 1>...};
        return fs[encodedParams](xParams);
    }

    template<typename uint>
    PwHammDistAlgorithm* runtimeDispatchAlgorithmInstance(ExperimentParams &xParams) {
        int encodedParams =
            xParams.shortCircuitMode << SHORTCIRCUIT_DISPATCH_BIT |
            xParams.binaryMode << BINARY_DISPATCH_BIT |
            xParams.groupedBFMode << GROUPED_DISPATCH_BIT |
            xParams.pivotsFilterMode << PIVOT_DISPATCH_BIT |
            xParams.nibblesMode << NIBBLE_DISPATCH_BIT;
        return runtimeDispatchAlgorithmInstance<uint>(xParams, encodedParams, std::make_index_sequence<DISPATCH_ENCODED_PARAMS_MAX>());
    }


public:
    BrutePwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory("brute-force") {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        switch(xParams.bytesPerElement) {
            case 0:
            case 1: return runtimeDispatchAlgorithmInstance<uint8_t>(xParams);
            case 2: return runtimeDispatchAlgorithmInstance<uint16_t>(xParams);
            case 4: return runtimeDispatchAlgorithmInstance<uint32_t>(xParams);
            default:
                fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                exit(EXIT_FAILURE);
        }
    }
};


#endif //PWHD_ALGORITHMFACTORY_H
