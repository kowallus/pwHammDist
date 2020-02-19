#ifndef PWHD_QUANTIZATIONALGORITHM_H
#define PWHD_QUANTIZATIONALGORITHM_H

#include "PwHammDistAlgorithmBase.h"

const string SIMPLE_BINARY_QUANTIZER_ID = "sbq";
const string THRESHOLD_STATS_BASED_QUANTIZER_ID = "tbq";
const string INTERVAL_STATS_BASED_QUANTIZER_ID = "ibq";

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
        } else {
            fprintf(stderr, "ERROR: unsupported bytes per element (%d) by simple quantization.\n", sizeof(uint));
            exit(EXIT_FAILURE);
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
    uint16_t half;

    void initQuantizeRules(ExperimentParams &xParams) {
        thresholds = new uint8_t[xParams.bytesPerSequence];
        half = xParams.d / 2 / this->SKIP_INTERVAL;
    }

    void buildQuantizeRules(uint16_t begCol, uint16_t endCol) {
        uint16_t* curColCounter = this->counter;
        for(uint16_t i = begCol; i < endCol; i++) {
            uint16_t cumm = 0;
            uint16_t j;
            for(j = 0; j < this->sigma; j++) {
                uint16_t newCumm = cumm + curColCounter[j];
                if (newCumm >= half) {
                    if (j && (newCumm - half > half - cumm))
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
        uint16_t i = 0;
        while(i < this->qxParams.m) {
            for(uint8_t b = 0; b < 8; b += 4, i += 4) {
                *dest += ((src[i] & this->bitMask) > thresholds[i]) << b;
                *dest += ((src[i + 1] & this->bitMask) > thresholds[i + 1]) << (b + 1);
                *dest += ((src[i + 2] & this->bitMask) > thresholds[i + 2]) << (b + 2);
                *dest += ((src[i + 3] & this->bitMask) > thresholds[i + 3]) << (b + 3);
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
class IntervalStatsBasedBinaryQuantizer: public StatsBasedBinaryQuantizer<uint> {
protected:

    pair<uint8_t, uint8_t>* intervals = 0;
    uint16_t half;

    void initQuantizeRules(ExperimentParams &xParams) {
        intervals = new pair<uint8_t, uint8_t>[xParams.bytesPerSequence];
        half = xParams.d / 2 / this->SKIP_INTERVAL;
    }

    void buildQuantizeRules(uint16_t begCol, uint16_t endCol) {
        uint16_t* curColCounter = this->counter;

        for(uint16_t i = begCol; i < endCol; i++) {
            uint16_t cumm = 0;
            uint16_t j;
            uint16_t start = 0;
            pair<uint8_t, uint8_t> interval(0, 0);
            uint16_t min = half;
            for(j = 0; j < this->sigma; j++) {
                uint16_t newCumm = cumm + curColCounter[j];
                if (newCumm >= half) {
                    if (min > half - cumm) {
                        min = half - cumm;
                        interval.first = start;
                        interval.second = j;
                    }
                    do {
                        cumm = newCumm;
                        newCumm -= curColCounter[start++];
                    } while (start < j && newCumm > half);
                    if (start <= j && min > cumm - half) {
                        min = cumm - half;
                        interval.first = start - 1;
                        interval.second = j + 1;
                    }
                }
                cumm = newCumm;
            }
            if (min > half - cumm) {
                min = half - cumm;
                interval.first = start;
                interval.second = j;
            }
            intervals[i] = interval;
            curColCounter += this->sigma;
        }
    }

    inline void quantizeSequence(uint8_t* dest, uint* src) {
        uint16_t i = 0;
        while(i < this->qxParams.m) {
            for(uint8_t b = 0; b < 8; b += 4, i += 4) {
                const uint8_t tmp0 = src[i] & this->bitMask;
                *dest += ((tmp0 >= intervals[i].first) && (tmp0 < intervals[i].second)) << b;
                const uint8_t tmp1 = src[i + 1] & this->bitMask;
                *dest += ((tmp1 >= intervals[i + 1].first) && (tmp1 < intervals[i + 1].second)) << (b + 1);
                const uint8_t tmp2 = src[i + 2] & this->bitMask;
                *dest += ((tmp2 >= intervals[i + 2].first) && (tmp2 < intervals[i + 2].second)) << (b + 2);
                const uint8_t tmp3 = src[i + 3] & this->bitMask;
                *dest += ((tmp3 >= intervals[i + 3].first) && (tmp3 < intervals[i + 3].second)) << (b + 3);
            }
            ++dest;
        }
    }

public:
    IntervalStatsBasedBinaryQuantizer(uint8_t testedBits) : StatsBasedBinaryQuantizer<uint>(testedBits) {}

    virtual ~IntervalStatsBasedBinaryQuantizer() {
        if (intervals) {
            delete[] intervals;
            intervals = 0;
        }
    }

    string getName() { return INTERVAL_STATS_BASED_QUANTIZER_ID + toString(this->testedBits); };

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
