#ifndef DBPIT_QUANTIZATIONSOLVER_H
#define DBPIT_QUANTIZATIONSOLVER_H

#include "SolverBase.h"

const string SIMPLE_BINARY_QUANTIZER_ID = "sbq";

const string QUATIZATION_BASED_SOLVER_ID = "qb";

template<typename uint>
class BinaryQuantizer {
protected:
    uint8_t* qSequences = 0;
    ExperimentParams qxParams;

    virtual void quantizeSequence(uint8_t* dest, uint* src) = 0;

public:
    BinaryQuantizer() {
        qxParams.enableBinaryMode();
    }

    void quantize(const uint8_t* sequences, ExperimentParams &xParams) {
        qxParams.m = xParams.m;
        qxParams.d = xParams.d;
        qxParams.k = xParams.k;
        qxParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(qxParams.m, qxParams.bitsPerPacked)
                                   * ceilDivisionBySmallInteger(qxParams.bitsPerPacked, 8);
        if (qxParams.alignSequences)
            qxParams.bytesPerSequence = (int) ceilDivisionBySmallInteger(qxParams.bytesPerSequence, 8) * 8;
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
class QuantizationBased_DBPI_Solver : public DBPISolver {
private:
    BinaryQuantizer<uint>* quantizer;
    SolverFactory* quantizedFilterSolverFactory;
    DBPISolver* postSolver;
    DBPISolver* qSolver = 0;

public:
    QuantizationBased_DBPI_Solver(ExperimentParams &xParams, BinaryQuantizer<uint> *quantizer,
                                  SolverFactory *quantizedFilterSolverFactory, DBPISolver *postSolver) : DBPISolver(
            xParams), quantizer(quantizer), quantizedFilterSolverFactory(quantizedFilterSolverFactory), postSolver(
            postSolver) {}

    virtual ~QuantizationBased_DBPI_Solver() {
        delete(quantizedFilterSolverFactory);
        delete(postSolver);
        delete(quantizer);
        if (qSolver)
            delete(qSolver);
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        cout << "checkpoint... " << " (" << time_millis() << " msec)" << endl;
        quantizer->quantize(sequences, xParams);
        cout << "quantized... " << " (" << time_millis() << " msec)" << endl;
        qSolver = quantizedFilterSolverFactory->getSolverInstance(quantizer->getQxParams());
        auto qRes = qSolver->findSimilarSequences(quantizer->getQSequences());
        cout << "filterCheck: " << qRes.size() << " (" << time_millis() << " msec)" << endl;
        vector<pair<uint16_t, uint16_t>> res = postSolver->findSimilarSequences(sequences, qRes);
        return res;
    };

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        return postSolver->findSimilarSequences(sequences, pairs);
    }

    inline bool testSequencesSimilarity(const uint8_t* sequences, uint16_t i, uint16_t j) {
        return postSolver->testSequencesSimilarity(sequences, i, j);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        return postSolver->testSequencesSimilarity(seq1, seq2);
    }

    string getName() {
        return QUATIZATION_BASED_SOLVER_ID + "+" +
            quantizer->getName() + "+" +
            (qSolver?qSolver->getName():"UNKNOWN") + "+" + postSolver->getName();
    }
};

#endif //DBPIT_QUANTIZATIONSOLVER_H
