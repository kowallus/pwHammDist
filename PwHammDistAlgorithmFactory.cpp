#include "PwHammDistAlgorithmFactory.h"

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams) {
    auto typeIt = pwHammDistAlgorithmTypesMap.find(xParams.algorithmID);
    if (typeIt == pwHammDistAlgorithmTypesMap.end()) {
        fprintf(stderr, "Invalid method type: %s\n", typeIt->second->getAlgorithmName().c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getAlgorithmInstance(xParams);
}


template<bool naive, bool nibble = false, bool grouped = false>
class BrutePwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    BrutePwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string(nibble?"grouped ":"") + string(nibble?"nibble ":"") + string(naive?"naive ":"short-circuit ") + string("brute-force")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        if(xParams.isInBinaryMode())
            return new BrutePwHammDistAlgorithm<naive, true>(xParams);
        else if (nibble) {
            switch(xParams.bytesPerElement) {
                case 1: return new NibbleBrutePwHammDistAlgorithm<naive, uint8_t>(xParams);
                case 2: return new NibbleBrutePwHammDistAlgorithm<naive, uint16_t>(xParams);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        } else {
            switch(xParams.bytesPerElement) {
                case 1: return new BrutePwHammDistAlgorithm<naive, false, uint8_t, grouped>(xParams);
                case 2: return new BrutePwHammDistAlgorithm<naive, false, uint16_t, grouped>(xParams);
                case 4: return new BrutePwHammDistAlgorithm<naive, false, uint32_t, grouped>(xParams);
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
            PwHammDistAlgorithm* postAlgorithm = BrutePwHammDistAlgorithmFactory<true, true>().getAlgorithmInstance(xParams);
            PwHammDistAlgorithmFactory* binaryAlgorithmFactory = new BrutePwHammDistAlgorithmFactory<false>();
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
        {{NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true>()},
         {SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false>()},
         {QUATIZATION_BASED_FILTER_PWHD_ID, new QuantizationBasedPwHammDistAlgorithmFactory()},
         {NIBBLE_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, true>()},
         {NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, true>()},
         {GROUPED_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, false, true>()},
         {GROUPED_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, false, true>()},
         {GROUPED_NIBBLE_NAIVE_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<true, true, true>()},
         {GROUPED_NIBBLE_SHORT_CIRCUIT_BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory<false, true, true>()}
         };

