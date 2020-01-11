#ifndef PWHD_ALGORITHMFACTORY_H
#define PWHD_ALGORITHMFACTORY_H

#include "QuantizationPwHammDistAlgorithm.h"
#include "xp-params.h"
#include <map>

extern map<string, PwHammDistAlgorithmFactory*> pwHammDistAlgorithmTypesMap;

PwHammDistAlgorithm* getPwHammDistAlgorithmInstance(ExperimentParams &xParams);

#endif //PWHD_ALGORITHMFACTORY_H
