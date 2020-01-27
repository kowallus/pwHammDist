#include "PwHammDistAlgorithmFactory.h"

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams) {
    auto typeIt = pwHammDistAlgorithmTypesMap.find(xParams.algorithmID);
    if (typeIt == pwHammDistAlgorithmTypesMap.end()) {
        fprintf(stderr, "Invalid method type ID: %s\n", xParams.algorithmID.c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getAlgorithmInstance(xParams);
}

class QuantizationBasedPwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    QuantizationBasedPwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string("quantization-based")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Quantization unsupported for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            PwHammDistAlgorithm* postAlgorithm = BrutePwHammDistAlgorithmFactory().getAlgorithmInstance(xParams);
            PwHammDistAlgorithmFactory* binaryAlgorithmFactory = new BrutePwHammDistAlgorithmFactory();
            switch(xParams.bytesPerElement) {
                case 1: return new QuantizationBasedPwHammDistAlgorithm<uint8_t>(xParams, new SimpleBinaryQuantizer<uint8_t>(), binaryAlgorithmFactory, postAlgorithm);
                case 2: return new QuantizationBasedPwHammDistAlgorithm<uint16_t>(xParams, new BitShiftBinaryQuantizer<uint16_t>(), binaryAlgorithmFactory, postAlgorithm);
                case 4: return new QuantizationBasedPwHammDistAlgorithm<uint32_t>(xParams, new SimpleBinaryQuantizer<uint32_t>(), binaryAlgorithmFactory, postAlgorithm);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap =
        {{BRUTE_FORCE_ID, new BrutePwHammDistAlgorithmFactory()},
         {QUATIZATION_BASED_FILTER_PWHD_ID, new QuantizationBasedPwHammDistAlgorithmFactory()},
         };

