#include "PwHammDistAlgorithmFactory.h"

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams xParams) {
    auto typeIt = pwHammDistAlgorithmTypesMap.find(xParams.algorithmID);
    if (typeIt == pwHammDistAlgorithmTypesMap.end()) {
        fprintf(stderr, "Invalid method type ID: %s\n", xParams.algorithmID.c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getAlgorithmInstance(xParams);
}

class QuantizationBasedPwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    QuantizationBasedPwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string("quantization-based filter")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Quantization unsupported for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            PwHammDistAlgorithm* preFilterAlgorithm;
            PwHammDistAlgorithmFactory* binaryAlgorithmFactory = new ConfigurablePwHammDistAlgorithmFactory();
            switch(xParams.bytesPerElement) {
                case 1: preFilterAlgorithm = new QuantizationBasedPwHammDistAlgorithm<uint8_t>(xParams, new SimpleBinaryQuantizer<uint8_t>(), binaryAlgorithmFactory);
                    break;
                case 2: preFilterAlgorithm = new QuantizationBasedPwHammDistAlgorithm<uint16_t>(xParams, new BitShiftBinaryQuantizer<uint16_t>(), binaryAlgorithmFactory);
                    break;
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
            xParams.disableFiltration();
            PwHammDistAlgorithm* postVerificationAlgorithm = ConfigurablePwHammDistAlgorithmFactory().getAlgorithmInstance(xParams);
            return new TwoLevelFilterBasedPwHammDistAlgorithm(xParams, preFilterAlgorithm, postVerificationAlgorithm);
        }
    }
};

class HashingBasedPwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
public:
    HashingBasedPwHammDistAlgorithmFactory():PwHammDistAlgorithmFactory(string("hashing-based filter")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Hashing not implemented for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            PwHammDistAlgorithm* preFilterAlgorithm;
            switch(xParams.bytesPerElement) {
                case 1: preFilterAlgorithm = new HashingBasedPwHammDistAlgorithm<uint8_t>(xParams);
                    break;
                case 2: preFilterAlgorithm = new HashingBasedPwHammDistAlgorithm<uint16_t>(xParams);
                    break;
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
            xParams.disableFiltration();
            PwHammDistAlgorithm* postVerificationAlgorithm = ConfigurablePwHammDistAlgorithmFactory().getAlgorithmInstance(xParams);
            return new TwoLevelFilterBasedPwHammDistAlgorithm(xParams, preFilterAlgorithm, postVerificationAlgorithm);
        }
    }
};


map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap =
        {{BRUTE_FORCE_ID, new ConfigurablePwHammDistAlgorithmFactory()},
         {QUATIZATION_BASED_FILTER_PWHD_ID, new QuantizationBasedPwHammDistAlgorithmFactory()},
         {HASHING_BASED_FILTER_PWHD_ID, new HashingBasedPwHammDistAlgorithmFactory()}
         };

