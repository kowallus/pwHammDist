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

    void cleanup() {
        if (qSequences) {
            delete[] qSequences;
            qSequences = 0;
        }
    }

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
            for(uint8_t b = 0; b < 8; b += 4, j += 4) {
                *dest += (src[j] & 1) << b;
                *dest += (src[j + 1] & 1) << (b + 1);
                *dest += (src[j + 2] & 1) << (b + 2);
                *dest += (src[j + 3] & 1) << (b + 3);
            }
            ++dest;
        }
    }

public:
    SimpleBinaryQuantizer<uint>(): BinaryQuantizer<uint>() {};

    string getName() { return SIMPLE_BINARY_QUANTIZER_ID; };
};

template<typename uint>
class BitShiftBinaryQuantizer: public BinaryQuantizer<uint> {
protected:
    inline void quantizeSequence(uint8_t* dest, uint* src) {
        uint64_t* srcPtr64 = (uint64_t*)src;
        uint64_t* tmp = (uint64_t*) dest;
        int end = this->qxParams.m / 64;
        if(sizeof(uint) == sizeof(uint16_t)) {
            for(int j = 0; j < end; j++, tmp++) {
                for(int shift = 0; shift < 16; ++shift)
                    *tmp = ((*tmp) << 1) + ((*srcPtr64++) & 0x0001000100010001);
            }
            end = (this->qxParams.m - end * 64) / 4;
            for(int shift = 0; shift < end; ++shift)
                *tmp = ((*tmp) << 1) + ((*srcPtr64++) & 0x0001000100010001);
        }
    }

public:
    BitShiftBinaryQuantizer<uint>(): BinaryQuantizer<uint>() {};

    string getName() { return SIMPLE_BINARY_QUANTIZER_ID; };
};

template<typename uint>
class QuantizationBasedPwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    BinaryQuantizer<uint>* quantizer;
    PwHammDistAlgorithmFactory* quantizedFilterAlgorithmFactory;
    PwHammDistAlgorithm* qAlgorithm = 0;

    void preprocessing(const uint8_t *sequences) {
        quantizer->quantize(sequences, xParams);
        if (xParams.verbose) cout << "quantized... " << " (" << time_millis() << " msec)" << endl;
        if (qAlgorithm) {
            delete (qAlgorithm);
            qAlgorithm = 0;
        }
        qAlgorithm = quantizedFilterAlgorithmFactory->getAlgorithmInstance(quantizer->getQxParams());
    };

    void postProcessing() {
        quantizer->cleanup();
    }

public:
    QuantizationBasedPwHammDistAlgorithm(ExperimentParams &xParams, BinaryQuantizer<uint> *quantizer,
            PwHammDistAlgorithmFactory *quantizedFilterAlgorithmFactory) : PwHammDistAlgorithm(
            xParams), quantizer(quantizer), quantizedFilterAlgorithmFactory(quantizedFilterAlgorithmFactory) {}

    virtual ~QuantizationBasedPwHammDistAlgorithm() {
        delete(quantizedFilterAlgorithmFactory);
        delete(quantizer);
        if (qAlgorithm)
            delete(qAlgorithm);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        preprocessing(sequences);
        auto qRes = qAlgorithm->findSimilarSequences(quantizer->getQSequences());
        postProcessing();
        return qRes;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        preprocessing(sequences);
        auto qRes = qAlgorithm->findSimilarSequences(quantizer->getQSequences(), pairs);
        postProcessing();
        return qRes;
    }

    string getName() {
        return QUATIZATION_BASED_FILTER_PWHD_ID + "+" +
            quantizer->getName() + "+" +
            (qAlgorithm?qAlgorithm->getName():"UNKNOWN");
    }
};

#endif //PWHD_QUANTIZATIONALGORITHM_H
