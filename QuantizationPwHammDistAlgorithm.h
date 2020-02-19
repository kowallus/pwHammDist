#ifndef PWHD_QUANTIZATIONALGORITHM_H
#define PWHD_QUANTIZATIONALGORITHM_H

#include "PwHammDistAlgorithmBase.h"

const string SIMPLE_BINARY_QUANTIZER_ID = "sbq";
const string THRESHOLD_STATS_BASED_QUANTIZER_ID = "tsbq";

const string QUATIZATION_BASED_FILTER_PWHD_ID = "qbf";

const string QUATIZATION_AND_PIVOTS_BASED_FILTER_PWHD_ID = "qpbf";

const string QUATIZATION_HASHING_BASED_FILTER_PWHD_ID = "qhbf";

template<typename uint>
class BinaryQuantizer {
protected:
    uint8_t* qSequences = 0;
    ExperimentParams qxParams;

    virtual void preprocess(uint8_t* sequences, ExperimentParams &xParams) {}
    virtual void quantizeSequence(uint8_t* dest, uint* src) = 0;

public:
    BinaryQuantizer() {}

    void quantize(uint8_t* sequences, ExperimentParams &xParams) {
        preprocess(sequences, xParams);
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
            uint* seq = (uint*) (sequences + (size_t) i * xParams.bytesPerSequence);
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
class StatsBasedBinaryQuantizer: public BinaryQuantizer<uint> {
protected:
    static const uint8_t COLUMNS_GROUP = 64;
    static const uint8_t MAX_TESTED_BITS = 8;
    static const uint16_t SKIP_INTERVAL = 16;
    uint16_t counter[COLUMNS_GROUP * (1 << MAX_TESTED_BITS)];

    const uint8_t testedBits;
    const uint16_t sigma;
    const uint bitMask;

    virtual void initQuantizeRules(ExperimentParams &xParams) = 0;
    virtual void buildQuantizeRules(uint16_t begCol, uint16_t endCol) = 0;

    void preprocess(uint8_t *sequences, ExperimentParams &xParams) {
        BinaryQuantizer<uint>::preprocess(sequences, xParams);
        this->initQuantizeRules(xParams);
        for(uint16_t i = 0; i < xParams.m; i += COLUMNS_GROUP) {
            memset(counter, 0, COLUMNS_GROUP * sigma * sizeof(uint));
            int endCol = i + COLUMNS_GROUP < xParams.m?i + COLUMNS_GROUP:xParams.m;
            uint8_t* x = sequences + i * sizeof(uint);
            for(uint16_t j = 0; j < xParams.d; j += SKIP_INTERVAL) {
                uint* y = (uint*) x;
                uint16_t* curColCounter = counter;
                for (uint16_t c = i; c < endCol; c++) {
                    curColCounter[*y++ & bitMask]++;
                    curColCounter += sigma;
                }
                x += xParams.bytesPerSequence * SKIP_INTERVAL;
            }
            this->buildQuantizeRules(i, endCol);
        }
        if (xParams.verbose)
            cout << "column stats... " << " (" << time_millis() << " msec)" << endl;
    }

public:
    StatsBasedBinaryQuantizer(uint8_t testedBits): testedBits(testedBits),
        sigma(1 << testedBits), bitMask(sigma - 1) {}

};

template<typename uint>
class ThresholdStatsBasedBinaryQuantizer: public StatsBasedBinaryQuantizer<uint> {
protected:

    uint8_t* thresholds = 0;
    uint16_t halfD;

    void initQuantizeRules(ExperimentParams &xParams) {
        thresholds = new uint8_t[xParams.bytesPerSequence * sizeof(uint)];
        halfD = xParams.d / 2 / this->SKIP_INTERVAL;
    }

    void buildQuantizeRules(uint16_t begCol, uint16_t endCol) {
        uint16_t* curColCounter = this->counter;
        for(uint16_t i = begCol; i < endCol; i++) {
            uint16_t cumm = 0;
            uint16_t j;
            for(j = 0; j < this->sigma; j++) {
                uint16_t newCumm = cumm + curColCounter[j];
                if (newCumm >= halfD) {
                    if (j && (newCumm - halfD > halfD - cumm))
                        j--;
                    break;
                }
                cumm = newCumm;
            }
            thresholds[i] = j;
            curColCounter += this->sigma;
        }
    }

    inline void quantizeSequence(uint8_t* dest, uint* src) {
        uint16_t j = 0;
        while(j < this->qxParams.m) {
            for(uint8_t b = 0; b < 8; b += 4, j += 4) {
                *dest += ((src[j] & this->bitMask) > thresholds[j]) << b;
                *dest += ((src[j + 1] & this->bitMask) > thresholds[j + 1]) << (b + 1);
                *dest += ((src[j + 2] & this->bitMask) > thresholds[j + 2]) << (b + 2);
                *dest += ((src[j + 3] & this->bitMask) > thresholds[j + 3]) << (b + 3);
            }
            ++dest;
        }
    }

public:
    ThresholdStatsBasedBinaryQuantizer(uint8_t testedBits) : StatsBasedBinaryQuantizer<uint>(testedBits) {}

    virtual ~ThresholdStatsBasedBinaryQuantizer() {
        if (thresholds) {
            delete[] thresholds;
            thresholds = 0;
        }
    }

    string getName() { return THRESHOLD_STATS_BASED_QUANTIZER_ID + toString(this->testedBits); };

};

template<typename uint>
class QuantizationBasedPwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    BinaryQuantizer<uint>* quantizer;
    PwHammDistAlgorithmFactory* quantizedFilterAlgorithmFactory;
    PwHammDistAlgorithm* qAlgorithm = 0;

    void preprocessing(uint8_t *sequences) {
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
        preprocessing((uint8_t*) sequences);
        auto qRes = qAlgorithm->findSimilarSequences(quantizer->getQSequences());
        postProcessing();
        return qRes;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        preprocessing((uint8_t*) sequences);
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
