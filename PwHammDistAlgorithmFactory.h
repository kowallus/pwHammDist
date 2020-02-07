#ifndef PWHD_ALGORITHMFACTORY_H
#define PWHD_ALGORITHMFACTORY_H

#include "QuantizationPwHammDistAlgorithm.h"
#include "xp-params.h"
#include <map>

extern map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap;

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams);

class ConfigurablePwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
private:
    static const int SHORTCIRCUIT_DISPATCH_BIT = 0;
    static const int COMPACT_DISPATCH_BIT = 1;
    static const int DISPATCH_ENCODED_PARAMS_MAX = 4;

    template <typename uint, bool shortcircuit, bool binaryAlphabet, bool compact>
    static PwHammDistAlgorithm* getBrutePwHammAlgorithmInstanceTemplate(ExperimentParams &xParams) {
        return new ConfigurablePwHammDistAlgorithm<uint, shortcircuit, binaryAlphabet, compact>(xParams);
    }

    template<typename uint, bool binaryMode, std::size_t...Is>
    PwHammDistAlgorithm* runtimeDispatchAlgorithmInstance(ExperimentParams &xParams, int encodedParams, std::index_sequence<Is...>) {
        using f_t = PwHammDistAlgorithm*(ExperimentParams &xParams);
        f_t* fs[] = {&getBrutePwHammAlgorithmInstanceTemplate<
                     uint, (Is >> SHORTCIRCUIT_DISPATCH_BIT) & 1,
                     binaryMode,
                     (Is >> COMPACT_DISPATCH_BIT) & 1>...};
        return fs[encodedParams](xParams);
    }

    template<bool binaryMode, typename uint>
    PwHammDistAlgorithm* runtimeDispatchAlgorithmInstance(ExperimentParams &xParams) {
        int encodedParams =
            xParams.shortCircuitMode << SHORTCIRCUIT_DISPATCH_BIT |
            xParams.compactMode << COMPACT_DISPATCH_BIT;
        return runtimeDispatchAlgorithmInstance<uint, binaryMode>(xParams, encodedParams, std::make_index_sequence<DISPATCH_ENCODED_PARAMS_MAX>());
    }


public:
    ConfigurablePwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory("brute-force") {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        if (xParams.binaryMode)
            return runtimeDispatchAlgorithmInstance<true, uint8_t>(xParams);
        switch(xParams.bytesPerElement) {
            case 1: return runtimeDispatchAlgorithmInstance<false, uint8_t>(xParams);
            case 2: return runtimeDispatchAlgorithmInstance<false, uint16_t>(xParams);
//            case 4: return runtimeDispatchAlgorithmInstance<uint32_t>(xParams);
            default:
                fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                exit(EXIT_FAILURE);
        }
    }
};


#endif //PWHD_ALGORITHMFACTORY_H
