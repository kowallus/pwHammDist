#include "PwHammDistAlgorithmFactory.h"

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams xParams) {
    auto typeIt = pwHammDistAlgorithmTypesMap.find(xParams.algorithmID);
    if (typeIt == pwHammDistAlgorithmTypesMap.end()) {
        fprintf(stderr, "Invalid method type ID: %s\n", xParams.algorithmID.c_str());
        exit(EXIT_FAILURE);
    }

    return typeIt->second->getAlgorithmInstance(xParams);
}

class QuantizationBasedPwHammDistPreFilterFactory: public PwHammDistAlgorithmFactory {
    PwHammDistAlgorithmFactory* binaryAlgorithmFactory;
    bool pivotsFilterMode;

public:
    QuantizationBasedPwHammDistPreFilterFactory(PwHammDistAlgorithmFactory* binaryAlgorithmFactory, bool pivotsFilterMode):
        PwHammDistAlgorithmFactory(string("quantization-based filter") + (pivotsFilterMode?" with pivots":"") +
            " and " + binaryAlgorithmFactory->getAlgorithmName()),
        binaryAlgorithmFactory(binaryAlgorithmFactory), pivotsFilterMode(pivotsFilterMode){};

    virtual ~QuantizationBasedPwHammDistPreFilterFactory() {
        delete(binaryAlgorithmFactory);
    }

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams xParams) {
        if(xParams.isInBinaryMode()) {
            fprintf(stderr, "Quantization unsupported for binary dataset.\n");
            exit(EXIT_FAILURE);
        } else {
            ExperimentParams qxParams = xParams;
            qxParams.pivotsFilterMode = pivotsFilterMode;
            switch(xParams.bytesPerElement) {
                case 1: return new QuantizationBasedPwHammDistAlgorithm<uint8_t>(qxParams, new SimpleBinaryQuantizer<uint8_t>(), binaryAlgorithmFactory);
                case 2: return new QuantizationBasedPwHammDistAlgorithm<uint16_t>(qxParams, new BitShiftBinaryQuantizer<uint16_t>(), binaryAlgorithmFactory);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

class HashingBasedPwHammDistPreFilterFactory: public PwHammDistAlgorithmFactory {
public:
    HashingBasedPwHammDistPreFilterFactory():PwHammDistAlgorithmFactory(string("hashing-based filter")) {};

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams xParams) {
        if(xParams.isInBinaryMode()) {
            ExperimentParams bxParams = xParams;
            bxParams.binaryMode = false;
            bxParams.alphabetSize = UINT8_MAX;
            bxParams.bytesPerElement = 1;
            bxParams.m = bxParams.bytesPerSequence;
            return new HashingBasedPwHammDistAlgorithm<uint8_t>(bxParams);
        } else {
            switch(xParams.bytesPerElement) {
                case 1: return new HashingBasedPwHammDistAlgorithm<uint8_t>(xParams);
                case 2: return new HashingBasedPwHammDistAlgorithm<uint16_t>(xParams);
                default:
                    fprintf(stderr, "ERROR: unsupported bytes per element: %d.\n", (int) xParams.bytesPerElement);
                    exit(EXIT_FAILURE);
            }
        }
    }
};

class FilterBasedPwHammDistAlgorithmFactory: public PwHammDistAlgorithmFactory {
private:
    PwHammDistAlgorithmFactory* preFilterAlgorithmFactory;
public:
    FilterBasedPwHammDistAlgorithmFactory(PwHammDistAlgorithmFactory* preFilterAlgorithmFactory):
        PwHammDistAlgorithmFactory(preFilterAlgorithmFactory->getAlgorithmName()),
        preFilterAlgorithmFactory(preFilterAlgorithmFactory){};

    virtual ~FilterBasedPwHammDistAlgorithmFactory() {
        delete(preFilterAlgorithmFactory);
    }

    PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams xParams) {
        PwHammDistAlgorithm* preFilterAlgorithm = preFilterAlgorithmFactory->getAlgorithmInstance(xParams);
        PwHammDistAlgorithm* postVerificationAlgorithm = ConfigurablePwHammDistAlgorithmFactory().getAlgorithmInstance(xParams);
        return new TwoLevelFilterBasedPwHammDistAlgorithm(xParams, preFilterAlgorithm, postVerificationAlgorithm);
    }
};

map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap =
        {{BRUTE_FORCE_ID, new ConfigurablePwHammDistAlgorithmFactory()}
         ,{QUATIZATION_BASED_FILTER_PWHD_ID, new FilterBasedPwHammDistAlgorithmFactory(
                 new QuantizationBasedPwHammDistPreFilterFactory(new ConfigurablePwHammDistAlgorithmFactory(), false))}
         ,{QUATIZATION_AND_PIVOTS_BASED_FILTER_PWHD_ID, new FilterBasedPwHammDistAlgorithmFactory(
                new QuantizationBasedPwHammDistPreFilterFactory(new ConfigurablePwHammDistAlgorithmFactory(), true))}
         ,{HASHING_BASED_FILTER_PWHD_ID, new FilterBasedPwHammDistAlgorithmFactory(new HashingBasedPwHammDistPreFilterFactory())}
         ,{QUATIZATION_HASHING_BASED_FILTER_PWHD_ID, new FilterBasedPwHammDistAlgorithmFactory(
                new QuantizationBasedPwHammDistPreFilterFactory(new HashingBasedPwHammDistPreFilterFactory(), false))}
         };

