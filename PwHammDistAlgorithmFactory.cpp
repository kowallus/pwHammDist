#include "PwHammDistAlgorithmFactory.h"

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams) {
    auto typeIt = pwHammDistAlgorithmTypesMap.find(xParams.algorithmID);
    if (typeIt == pwHammDistAlgorithmTypesMap.end()) {
        fprintf(stderr, "Invalid method type ID: %s\n", xParams.algorithmID.c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getAlgorithmInstance(xParams);
}


template<bool shortcircuit, bool nibble = false, bool grouped = false, bool pivot = false>
class BrutePwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    BrutePwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string(pivot?"pivot-based ":"") + string(grouped?"grouped ":"") + string(nibble?"nibble ":"") + string(shortcircuit?"short-circuit ":"naive ") + string("brute-force")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        if(xParams.isInBinaryMode())
            return new BrutePwHammDistAlgorithm<shortcircuit, true>(xParams);
        else if (nibble) {
            switch(xParams.bytesPerElement) {
                case 1: return new NibbleBrutePwHammDistAlgorithm<!shortcircuit, uint8_t>(xParams);
                case 2: return new NibbleBrutePwHammDistAlgorithm<!shortcircuit, uint16_t>(xParams);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        } else {
            switch(xParams.bytesPerElement) {
                case 1: return new BrutePwHammDistAlgorithm<shortcircuit, false, uint8_t, grouped, pivot>(xParams);
                case 2: return new BrutePwHammDistAlgorithm<shortcircuit, false, uint16_t, grouped, pivot>(xParams);
                case 4: return new BrutePwHammDistAlgorithm<shortcircuit, false, uint32_t, grouped, pivot>(xParams);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

class QuantizationBasedPwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    QuantizationBasedPwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string("quantization-based")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Quantization unsupported for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            PwHammDistAlgorithm* postAlgorithm = (xParams.alphabetSize < 8?
                    BrutePwHammDistAlgorithmFactory<true, true>().getAlgorithmInstance(xParams):
                    BrutePwHammDistAlgorithmFactory<true, false>().getAlgorithmInstance(xParams));
            PwHammDistAlgorithmFactory* binaryAlgorithmFactory = new BrutePwHammDistAlgorithmFactory<true>();
            switch(xParams.bytesPerElement) {
                case 1: return new QuantizationBasedPwHammDistAlgorithm<uint8_t>(xParams, new SimpleBinaryQuantizer<uint8_t>(), binaryAlgorithmFactory, postAlgorithm);
                case 2: return new QuantizationBasedPwHammDistAlgorithm<uint16_t>(xParams, new SimpleBinaryQuantizer<uint16_t>(), binaryAlgorithmFactory, postAlgorithm);
                case 4: return new QuantizationBasedPwHammDistAlgorithm<uint32_t>(xParams, new SimpleBinaryQuantizer<uint32_t>(), binaryAlgorithmFactory, postAlgorithm);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap =
        {{NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false>()},
         {SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true>()},
         {QUATIZATION_BASED_FILTER_PWHD_ID, new QuantizationBasedPwHammDistAlgorithmFactory()},
         {NIBBLE_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, true>()},
         {NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, true>()},
         {GROUPED_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, false, true>()},
         {GROUPED_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, false, true>()},
//         {GROUPED_NIBBLE_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, true, true>()},
//         {GROUPED_NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, true, true>()},
         {PIVOT_FILTER_ID, new BrutePwHammDistAlgorithmFactory<false, false, false, true>()},
         {SHORT_CIRCUIT_PIVOT_FILTER_ID, new BrutePwHammDistAlgorithmFactory<true, false, false, true>()},
         {GROUPED_PIVOT_FILTER_ID, new BrutePwHammDistAlgorithmFactory<false, false, true, true>()},
         {GROUPED_SHORT_CIRCUIT_PIVOT_FILTER_ID, new BrutePwHammDistAlgorithmFactory<true, false, true, true>()},
         };

