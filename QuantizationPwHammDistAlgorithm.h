#ifndef PWHD_QUANTIZATIONALGORITHM_H
#define PWHD_QUANTIZATIONALGORITHM_H

#include "PwHammDistAlgorithmBase.h"

const string SIMPLE_BINARY_QUANTIZER_ID = "sbq";

const string QUATIZATION_BASED_FILTER_PWHD_ID = "qbf";

template<typename uint>
class BinaryQuantizer {
protected:
    uint8_t* qSequences = 0;
    ExperimentParams qxParams;

    virtual void quantizeSequence(uint8_t* dest, uint* src) = 0;

public:
    BinaryQuantizer() {}

    void quantize(const uint8_t* sequences, ExperimentParams &xParams) {
        qxParams = xParams;
        qxParams.enableBinaryMode();
        qxParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(qxParams.m, qxParams.bitsPerPacked)
                                   * ceilDivisionBySmallInteger(qxParams.bitsPerPacked, 8);
        if (qxParams.alignSequencesTo256bits)
            qxParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(qxParams.bytesPerSequence,
                    ExperimentParams::ALINGMENT_IN_BYTES) * ExperimentParams::ALINGMENT_IN_BYTES;
        qSequences = new uint8_t[(size_t) qxParams.d * qxParams.bytesPerSequence]();
        uint8_t* qCur;
        for(uint16_t i = 0; i < qxParams.d; i++) {
            uint* seq = (uint*) ((uint8_t*) sequences + (size_t) i * xParams.bytesPerSequence);
            qCur = (qSequences + i * qxParams.bytesPerSequence);
            quantizeSequence(qCur, seq);
        }
    };

    virtual string getName() = 0;

    virtual ~BinaryQuantizer() {
        if (qSequences)
            delete[] qSequences;
    }

    uint8_t *getQSequences() const {
        return qSequences;
    }

    ExperimentParams &getQxParams() {
        return qxParams;
    }
};

template<typename uint>
class SimpleBinaryQuantizer: public BinaryQuantizer<uint> {
protected:
    inline void quantizeSequence(uint8_t* dest, uint* src) {
        uint16_t j = 0;
        while(j < this->qxParams.m) {
            for(uint8_t b = 0; b < 8; b++)
                *dest += (src[j++] & 1) << b;
            ++dest;
        }
    }

public:
    SimpleBinaryQuantizer<uint>(): BinaryQuantizer<uint>() {};

    string getName() { return SIMPLE_BINARY_QUANTIZER_ID; };

};

template<typename uint>
class QuantizationBasedPwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    BinaryQuantizer<uint>* quantizer;
    PwHammDistAlgorithmFactory* quantizedFilterAlgorithmFactory;
    PwHammDistAlgorithm* postAlgorithm;
    PwHammDistAlgorithm* qAlgorithm = 0;

public:
    QuantizationBasedPwHammDistAlgorithm(ExperimentParams &xParams, BinaryQuantizer<uint> *quantizer,
            PwHammDistAlgorithmFactory *quantizedFilterAlgorithmFactory, PwHammDistAlgorithm *postAlgorithm) : PwHammDistAlgorithm(
            xParams), quantizer(quantizer), quantizedFilterAlgorithmFactory(quantizedFilterAlgorithmFactory), postAlgorithm(
            postAlgorithm) {}

    virtual ~QuantizationBasedPwHammDistAlgorithm() {
        delete(quantizedFilterAlgorithmFactory);
        delete(postAlgorithm);
        delete(quantizer);
        if (qAlgorithm)
            delete(qAlgorithm);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        if (xParams.verbose) cout << "checkpoint... " << " (" << time_millis() << " msec)" << endl;
        quantizer->quantize(sequences, xParams);
        if (xParams.verbose) cout << "quantized... " << " (" << time_millis() << " msec)" << endl;
        qAlgorithm = quantizedFilterAlgorithmFactory->getAlgorithmInstance(quantizer->getQxParams());
        auto qRes = qAlgorithm->findSimilarSequences(quantizer->getQSequences());
        if (xParams.verbose) cout << "filterCheck: " << qRes.size() << " (" << time_millis() << " msec)" << endl;
        vector<pair<uint16_t, uint16_t>> res = postAlgorithm->findSimilarSequences(sequences, qRes);
        return res;
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        return postAlgorithm->findSimilarSequences(sequences, pairs);
    }

    string getName() {
        return QUATIZATION_BASED_FILTER_PWHD_ID + "+" +
            quantizer->getName() + "+" +
            (qAlgorithm?qAlgorithm->getName():"UNKNOWN") + "+" + postAlgorithm->getName();
    }
};

#endif //PWHD_QUANTIZATIONALGORITHM_H
