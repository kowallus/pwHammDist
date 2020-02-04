#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "utils/helper.h"
#include "xp-params.h"
#include <vector>
#include <cassert>
#include <random>

using namespace std;

const string BRUTE_FORCE_ID = "bf";
const string NO_SHORT_CIRCUIT_PREFIX_ID = "n";
const string GROUPED_PREFIX_ID = "g";
const string COMPACT_PREFIX_ID = "c";
const string INTERLEAVE_BITS_PREFIX_ID = "i";
const string LAZY_INTERLEAVE_BITS_PREFIX_ID = "I";
const string PIVOT_FILTER_PREFIX_ID = "p";
const string ELECTION_PIVOT_FILTER_PREFIX_ID = "P";
const string BINARY_MODE_ID_SUFFIX = "_bin";


class PwHammDistAlgorithm {
protected:
    ExperimentParams &xParams;
    PwHammDistAlgorithm(ExperimentParams &xParams): xParams(xParams) {};

public:
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t *sequences) = 0;
    virtual vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
            const vector<pair<uint16_t, uint16_t>> pairs) = 0;
    virtual string getName() = 0;
};

class PwHammDistAlgorithmFactory {
private:
    string algorithmName;

public:
    PwHammDistAlgorithmFactory(string algorithmName): algorithmName(algorithmName) { }
    string getAlgorithmName() { return algorithmName; }
    virtual PwHammDistAlgorithm* getAlgorithmInstance(ExperimentParams &xParams) = 0;
};

template <typename uint, bool shortcircuit, bool binaryAlphabet, bool compact = false>
class BrutePwHammDistAlgorithm : public PwHammDistAlgorithm {
private:
    const uint16_t seqInULLs;

    uint8_t* origSeq = 0;
    uint16_t origBytesPerSequence;
    uint8_t* seq1 = 0;
    uint8_t* seq2 = 0;

    // pivot mode params
    static const uint8_t PIVOTS_COUNT_MAX = 6;
    uint8_t pivotsCount = PIVOTS_COUNT_MAX;
    uint16_t pivots[PIVOTS_COUNT_MAX];
    uint16_t* pivotDist = 0;
    uint8_t ctrlPivot = 0;

    // compact (or nibble) and interleaved modes
    uint8_t* seqAugmention = 0;
    uint16_t augSeqPerBitInULLs;
    uint8_t bitsPerElement;

    void preprocessing(const uint8_t *sequences) {
        origSeq = (uint8_t*) sequences;
        origBytesPerSequence = xParams.bytesPerSequence;
        if (compact) {
            compactation(sequences);
        } else if (xParams.interleaveBitsMode) {
            interleaveBits(sequences);
        } else {
            seq1 = (uint8_t*) sequences;
            seq2 = (uint8_t*) sequences;
        }
        if (xParams.pivotsFilterMode) {
            if (xParams.verbose) cout << "pivots: ";
            if (xParams.pivotsElectionMode)
                calculateDistancesToPivotsWithStats();
            else
                calculateDistancesToRandomPivots();
            if (xParams.verbose) cout << "..." << " (" << time_millis() << " msec)" << endl;
        }
    }

    void postprocessing() {
        xParams.bytesPerSequence = origBytesPerSequence;
        if (seqAugmention) {
            delete[] seqAugmention;
            seqAugmention = 0;
        }
        if (pivotDist) {
            delete[] pivotDist;
            pivotDist = 0;
        }
    }


    template<bool allowShortCircuit>
    uint16_t findSequencesDistance(const void* seq1, const void* seq2) {
        uint16_t dist;
        if (binaryAlphabet)
            dist = (allowShortCircuit && shortcircuit)?hammingDistanceBinary((uint64_t*) seq1, (uint64_t*) seq2, seqInULLs, xParams.k):
                   hammingDistanceBinary((uint64_t*) seq1, (uint64_t*)seq2, seqInULLs);
        else {
            if (compact) {
                if (sizeof(uint) == sizeof(uint8_t)) {
                    assert(xParams.alphabetSize <= 8);
                    dist = (allowShortCircuit && shortcircuit)? hammingDistanceAugmentedNibble((uint64_t *) seq1, (uint64_t *) seq2,
                                                          xParams.bytesPerSequence * 2, xParams.k):
                                        hammingDistanceAugmentedNibble((uint64_t *) seq1, (uint64_t *) seq2,
                                                          xParams.bytesPerSequence * 2);
                } else if (sizeof(uint) == sizeof(uint16_t)) {
                    dist = (allowShortCircuit && shortcircuit)? hammingDistanceAugmented16bit((uint64_t *) seq1, (uint64_t *) seq2,
                                                                 xParams.bytesPerSequence / 2, xParams.k):
                                        hammingDistanceAugmented16bit((uint64_t *) seq1, (uint64_t *) seq2,
                                                                 xParams.bytesPerSequence / 2);
                }
            } else {
                if (xParams.interleaveBitsMode) {
                    if (xParams.lazyInterleaveBitsMode)
                        dist = (allowShortCircuit && shortcircuit) ? lazyHammingInterleavedBitsDistance((uint64_t *) seq1, (uint64_t *) seq2,
                            augSeqPerBitInULLs, bitsPerElement, xParams.k):
                                lazyHammingInterleavedBitsDistance((uint64_t *) seq1, (uint64_t *) seq2, augSeqPerBitInULLs , bitsPerElement);
                    else
                        dist = (allowShortCircuit && shortcircuit) ? hammingInterleavedBitsDistance((uint64_t *) seq1, (uint64_t *) seq2,
                                                                                                    augSeqPerBitInULLs, bitsPerElement, xParams.k):
                               hammingInterleavedBitsDistance((uint64_t *) seq1, (uint64_t *) seq2, augSeqPerBitInULLs, bitsPerElement);
                } else {
                    dist = (allowShortCircuit && shortcircuit) ? hammingDistance((uint *) seq1, (uint *) seq2, xParams.m,
                                                                                 xParams.k) :
                           hammingDistance((uint *) seq1, (uint *) seq2, xParams.m);
                }
            }
        }
        return dist;
    }

    void compactation(const uint8_t *sequences) {
        if (sizeof(uint) == sizeof(uint8_t) && xParams.alphabetSize <= 8) {
            xParams.bytesPerSequence /= 2;
            seqAugmention = new uint8_t[xParams.d * xParams.bytesPerSequence * 2]();
            uint8_t* x = (uint8_t*) sequences;
            uint8_t *y = seqAugmention;
            uint8_t *z = seqAugmention + xParams.d * xParams.bytesPerSequence;
            const uint32_t nibblePackedBytes = xParams.d * xParams.bytesPerSequence;
            for(int i = 0; i < nibblePackedBytes; i++) {
                *y = *(x++);
                *y += *(x++)*16;
                *(z++) = *(y++) + 128 + 8;
            }
            seq1 = seqAugmention;
            seq2 = seqAugmention + xParams.d * xParams.bytesPerSequence;
            if (xParams.verbose) cout << "nibbled... " << " (" << time_millis() << " msec)" << endl;
        } else if (sizeof(uint) == sizeof(uint16_t)) {
            seqAugmention = new uint8_t[xParams.d * xParams.bytesPerSequence]();
            uint16_t* x = (uint16_t*) sequences;
            uint16_t* z = (uint16_t*) seqAugmention;
            const uint32_t length = xParams.d * xParams.bytesPerSequence / 2;
            for(int i = 0; i < length; i++) {
                *(z++) = *(x++) + 32768;
            }
            seq1 = (uint8_t*) sequences;
            seq2 = seqAugmention;
            if (xParams.verbose) cout << "augmented... " << " (" << time_millis() << " msec)" << endl;
        } else {
            fprintf(stderr, "ERROR: unsupported alphabet size for nibblification.\n");
            exit(EXIT_FAILURE);
        }
    }

    inline uint64_t lazyHammingInterleavedBitsDistance(const uint64_t *x, const uint64_t *y, int lengthInULLs, int bitsPerElement)
    {
        lengthInULLs--;
        assert(lengthInULLs % 4 == 0);
        x++; y++;
        for (int i = 0; i < lengthInULLs; i += 4) {
            mismatchFlags[i] = x[i] ^ y[i];
            mismatchFlags[i + 1] = x[i + 1] ^ y[i + 1];
            mismatchFlags[i + 2] = x[i + 2] ^ y[i + 2];
            mismatchFlags[i + 3] = x[i + 3] ^ y[i + 3];
        }
        for(int b = 1; b < bitsPerElement; b++) {
            x += lengthInULLs;
            y += lengthInULLs;
            lazyInterleaveBitsInSequence((uint8_t*) (x++), b);
            lazyInterleaveBitsInSequence((uint8_t*) (y++), b);
            for (int i = 0; i < lengthInULLs; i += 4) {
                mismatchFlags[i] |= x[i] ^ y[i];
                mismatchFlags[i + 1] |= x[i + 1] ^ y[i + 1];
                mismatchFlags[i + 2] |= x[i + 2] ^ y[i + 2];
                mismatchFlags[i + 3] |= x[i + 3] ^ y[i + 3];
            }
        }
        uint64_t res = 0;
        for (int i = 0; i < lengthInULLs; i += 4) {
            res += __builtin_popcountll(mismatchFlags[i]) + __builtin_popcountll(mismatchFlags[i + 1]) +
                   __builtin_popcountll(mismatchFlags[i + 2]) + __builtin_popcountll(mismatchFlags[i + 3]);
        }
        return res;
    }

    inline uint64_t lazyHammingInterleavedBitsDistance(const uint64_t *x, const uint64_t *y, int lengthInULLs, int bitsPerElement, int limit)
    {
        lengthInULLs--;
        assert(lengthInULLs % 4 == 0);
        x++; y++;
        uint64_t res = 0;
        for (int i = 0; i < lengthInULLs; i += 4) {
            mismatchFlags[i] = x[i] ^ y[i];
            mismatchFlags[i + 1] = x[i + 1] ^ y[i + 1];
            mismatchFlags[i + 2] = x[i + 2] ^ y[i + 2];
            mismatchFlags[i + 3] = x[i + 3] ^ y[i + 3];
            res += __builtin_popcountll(mismatchFlags[i]) + __builtin_popcountll(mismatchFlags[i + 1]) +
                   __builtin_popcountll(mismatchFlags[i + 2]) + __builtin_popcountll(mismatchFlags[i + 3]);
            if (res > limit)
                return res;
        }
        for(int b = 1; b < bitsPerElement; b++) {
            x += lengthInULLs;
            y += lengthInULLs;
            lazyInterleaveBitsInSequence((uint8_t*) (x++), b);
            lazyInterleaveBitsInSequence((uint8_t*) (y++), b);
            res = 0;
            for (int i = 0; i < lengthInULLs; i += 4) {
                mismatchFlags[i] |= x[i] ^ y[i];
                mismatchFlags[i + 1] |= x[i + 1] ^ y[i + 1];
                mismatchFlags[i + 2] |= x[i + 2] ^ y[i + 2];
                mismatchFlags[i + 3] |= x[i + 3] ^ y[i + 3];
                res += __builtin_popcountll(mismatchFlags[i]) + __builtin_popcountll(mismatchFlags[i + 1]) +
                       __builtin_popcountll(mismatchFlags[i + 2]) + __builtin_popcountll(mismatchFlags[i + 3]);
            }
            if (res > limit)
                return res;
        }
        return res;
    }

    void lazyInterleaveBitsInSequence(uint8_t* dest, int b) {
        if (*(uint64_t*) dest == 0)
            return;
        uint* seq = *(uint**) dest;
        uint64_t* z = (uint64_t*) dest;
        *z++ = 0;
        uint64_t mask = sizeof(uint) == 2 ? 0x0001000100010001 : 0x0101010101010101;
        mask <<= b;
        int end = this->xParams.m / 64;
        uint64_t *srcPtr64 = (uint64_t *) seq;
        for(int j = 0; j < end; j++, z++)
            for (int shift = 0; shift < 8 * sizeof(uint); ++shift)
                *z = ((*z) << 1) + (((*srcPtr64++) & mask) >> b);
        end = ceilDivisionBySmallInteger(this->xParams.m - end * 64, 8 / sizeof(uint));
        for(int shift = 0; shift < end; ++shift)
            *z = ((*z) << 1) + (((*srcPtr64++) & mask) >> b);
    }

    void interleaveBitsInSequence(uint8_t* dest, uint* seq) {
        uint64_t* z[32];
        z[0] = (uint64_t*) dest;
        uint64_t mask[32];
        mask[0] = sizeof(uint) == 2 ? 0x0001000100010001 : 0x0101010101010101;
        int up2BitsPerElement = xParams.lazyInterleaveBitsMode?1:bitsPerElement;
        for(int b = 1; b < up2BitsPerElement; b++) {
            z[b] = z[b - 1] + augSeqPerBitInULLs;
            mask[b] = mask[b - 1] * 2;
        }

        int end = this->xParams.m / 64;
        uint64_t *srcPtr64 = (uint64_t *) seq;
        if (xParams.lazyInterleaveBitsMode) {
            for (int b = 1; b < bitsPerElement; b++) {
                *((uint **) (z[0] + b * augSeqPerBitInULLs)) = seq;
            }
            z[0]++;
        }
        for(int j = 0; j < end; j++) {
            for (int shift = 0; shift < 8 * sizeof(uint); ++shift) {
                for (int b = 0; b < up2BitsPerElement; b++)
                    *z[b] = ((*z[b]) << 1) + (((*srcPtr64) & mask[b]) >> b);
                srcPtr64++;
            }
            for(int b = 0; b < up2BitsPerElement; b++)
                z[b]++;
        }

        end = ceilDivisionBySmallInteger(this->xParams.m - end * 64, 8 / sizeof(uint));
        for(int shift = 0; shift < end; ++shift) {
            for (int b = 0; b < up2BitsPerElement; b++)
                *z[b] = ((*z[b]) << 1) + (((*srcPtr64) & mask[b]) >> b);
            srcPtr64++;
        }
    }

    void interleaveBits(const uint8_t *sequences) {
        this->bitsPerElement = 32 - __builtin_clz(xParams.alphabetSize - 1);
        augSeqPerBitInULLs = ceilDivisionBySmallInteger(xParams.m, 64);
        augSeqPerBitInULLs = ceilDivisionBySmallInteger(augSeqPerBitInULLs, 4) * 4;
        if (xParams.lazyInterleaveBitsMode)
            augSeqPerBitInULLs++;
        xParams.bytesPerSequence = bitsPerElement * augSeqPerBitInULLs * 8;
        seqAugmention = new uint8_t[xParams.d * xParams.bytesPerSequence]();

        uint8_t* dest = seqAugmention;
        uint8_t* seq = (uint8_t*) sequences;
        for(uint16_t i = 0; i < xParams.d; i++) {
            interleaveBitsInSequence(dest, (uint*) seq);
            seq += origBytesPerSequence;
            dest += xParams.bytesPerSequence;
        }
        seq1 = seqAugmention;
        seq2 = seqAugmention;
        if (xParams.verbose) cout << "interleaved bits... " << " (" << time_millis() << " msec)" << endl;
    }

    void calculateDistancesToRandomPivots() {
        randgenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
        uint16_t candidate = randgenerator() % xParams.d;
        uint16_t maxUp2HalfKCount = 0;
        pivotDist = new uint16_t[xParams.d * PIVOTS_COUNT_MAX];
        uint16_t* ptr = pivotDist;
        for(int p = 0; p < pivotsCount; p++) {
            pivots[p] = candidate;
            if (xParams.verbose) cout << candidate << "\t";
            uint8_t* x = seq1 + (size_t) pivots[p] * xParams.bytesPerSequence;
            uint8_t* y = seq2;
            int up2HalfKCount = 0;
            for(int j = 0; j < xParams.d; j++) {
                ptr[j] = findSequencesDistance<false>(x, y);
                y += xParams.bytesPerSequence;
                if (ptr[j] < xParams.k / 2)
                    up2HalfKCount++;
            }
            cout << "(" << up2HalfKCount << ") ";
            if (maxUp2HalfKCount < up2HalfKCount) {
                maxUp2HalfKCount = up2HalfKCount;
                ctrlPivot = p;
            }
            if (p > xParams.d) {
                pivotsCount = p;
                break;
            }
            bool lackOfCandidate = true;
            while (lackOfCandidate) {
                candidate = randgenerator() % xParams.d;
                lackOfCandidate = false;
                for (int p2 = 0; p2 <= p; p2++)
                    if (pivots[p2] == candidate) {
                        lackOfCandidate = true;
                        break;
                    }
            }
            ptr += xParams.d;
        }
    }

    template<bool maxSumPivots = false>
    void calculateDistancesToPivotsWithStats() {
        randgenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
        uint16_t candidate = randgenerator() % xParams.d;
        pivotDist = new uint16_t[xParams.d * pivotsCount];
        uint16_t minPivotDist[xParams.d];
        memset(&minPivotDist, UINT8_MAX, sizeof(minPivotDist));
        const uint16_t minThreshold = xParams.k;
        uint32_t sumPivotDist[xParams.d] = { 0 };
        uint16_t maxUp2HalfKCount = 0;
        uint16_t* ptr = pivotDist;
        for(int p = 0; p < pivotsCount; p++) {
            pivots[p] = candidate;
            if (xParams.verbose) cout << pivots[p] << "\t";
            uint32_t maxSumDist = 0;
            uint8_t* x = seq1 + (size_t) pivots[p] * xParams.bytesPerSequence;
            uint8_t* y = seq2;
            int up2HalfKCount = 0;
            for(int j = 0; j < xParams.d; j++) {
                ptr[j] = findSequencesDistance<false>(x, y);
                y += xParams.bytesPerSequence;
                if (ptr[j] < xParams.k / 2)
                    up2HalfKCount++;
                if (minPivotDist[j] > ptr[j])
                    minPivotDist[j] = ptr[j];
                if (minPivotDist[j] > minThreshold) {
                    sumPivotDist[j] += ptr[j];
                    if (maxSumDist < sumPivotDist[j]) {
                        maxSumDist = sumPivotDist[j];
                        candidate = j;
                    }
                }
            }
            cout << "(" << up2HalfKCount << ") ";
            if (maxUp2HalfKCount < up2HalfKCount) {
                maxUp2HalfKCount = up2HalfKCount;
                ctrlPivot = p;
            }
            if (pivots[p] == candidate) {
                pivotsCount = p + 1;
                break;
            }
            if (!maxSumPivots && p + 1 < pivotsCount) {
                bool lackOfCandidate = true;
                while (lackOfCandidate) {
                    bool badCandidate = true;
                    do {
                        candidate = randgenerator() % xParams.d;
                        badCandidate = minPivotDist[candidate] < minThreshold;
                        if (badCandidate) cout << ".";
                    } while (badCandidate);
                    lackOfCandidate = false;
                    for (int p2 = 0; p2 <= p; p2++)
                        if (pivots[p2] == candidate) {
                            lackOfCandidate = true;
                            break;
                        }
                }
            }
            ptr += xParams.d;
        }
    }

    enum filterResult { similar, inconclusive, different };

    filterResult pivotFilter(uint16_t i, uint16_t j) {
        for(int p = 0; p < pivotsCount; p++) {
            const int iDist = pivotDist[p * xParams.d + i];
            const int jDist = pivotDist[p * xParams.d + j];
            if (abs(iDist - jDist) > xParams.k)
                return different;
            if (iDist + jDist <= xParams.k)
                return similar;
            
        }
        return inconclusive;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingControlPivot() {
        vector<uint16_t> pivotRank(xParams.d);
        for(int i = 0; i < xParams.d; i++)
            pivotRank[i] = i;
        uint16_t* ctrlPivotDist = pivotDist + ctrlPivot * xParams.d;
        std::sort(pivotRank.begin(), pivotRank.end(), [ctrlPivotDist](const uint16_t& idx1, const uint16_t& idx2) -> bool
            { return ctrlPivotDist[idx1] < ctrlPivotDist[idx2]; });

        vector<pair<uint16_t, uint16_t>> res;
        for(int ri = 0; ri < xParams.d - 1; ri++) {
            int rj = ri + 1;
            const int i = pivotRank[ri];
            uint16_t iDist[PIVOTS_COUNT_MAX];
            for(int p = 0; p < pivotsCount; p++)
                iDist[p] = pivotDist[p * xParams.d + i];
            while (rj < xParams.d && iDist[ctrlPivot] + ctrlPivotDist[pivotRank[rj]] <= xParams.k)
                res.push_back(pair<uint16_t, uint16_t>(i, pivotRank[rj++]));
            while (rj < xParams.d && ctrlPivotDist[pivotRank[rj]] - iDist[ctrlPivot] <= xParams.k) {
                const int j = pivotRank[rj];
                filterResult filterRes = inconclusive;
                for(int p = 0; p < pivotsCount; p++) {
                    if (p == ctrlPivot)
                        continue;
                    const int jDist = pivotDist[p * xParams.d + j];
                    if (abs(iDist[p] - jDist) > xParams.k) {
                        filterRes = different;
                        break;
                    } if (iDist[p] + jDist <= xParams.k) {
                        filterRes = similar;
                        break;
                    }
                }
                if (filterRes == similar || (filterRes == inconclusive && testSequencesSimilarity<false>(i, j)))
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                rj++;
            }
        }
        return res;
    };


    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingStandardBrute() {
        vector<pair<uint16_t, uint16_t>> res;
        for(int i = 0; i < xParams.d - 1; i++) {
            for (int j = i + 1; j < xParams.d; j++) {
                if (testSequencesSimilarity<false>(i, j)) {
                    res.push_back(pair<uint16_t, uint16_t>(i, j));
                }
            }
        }
        return res;
    };

    template<bool usePivotFilter>
    vector<pair<uint16_t, uint16_t>> findSimilarSequencesUsingGroupedApproach() {
        vector<pair<uint16_t, uint16_t>> res;
        const int GROUP_SIZE = 32;
        int gEnd = (xParams.d / GROUP_SIZE) * GROUP_SIZE;
        {
            for (int gStart = 0; gStart < gEnd; gStart += GROUP_SIZE) {
                for (int i = 0; i < GROUP_SIZE; i++) {
                    for (int j = i + 1; j < GROUP_SIZE; j++) {
                        if (testSequencesSimilarity<usePivotFilter>(gStart + i, gStart + j)) {
                            res.push_back(pair<uint16_t, uint16_t>(gStart + i, gStart + j));
                        }
                    }
                }
            }
        }
        {
            for(int g1Start = 0; g1Start < gEnd; g1Start += GROUP_SIZE) {
                for (int g2Start = g1Start + GROUP_SIZE; g2Start < gEnd; g2Start += GROUP_SIZE) {
                    for (int i = 0; i < GROUP_SIZE; i++) {
                        for (int j = 0; j < GROUP_SIZE; j++) {
                            if (testSequencesSimilarity<usePivotFilter>(g1Start + i, g2Start + j)) {
                                res.push_back(pair<uint16_t, uint16_t>(g1Start + i, g2Start + j));
                            }
                        }
                    }
                }
            }
        }
        {
            for(int i = 0; i < xParams.d - 1; i++) {
                int j = i < gEnd?gEnd:i + 1;
                for (; j < xParams.d; j++) {
                    if (testSequencesSimilarity<usePivotFilter>(i, j)) {
                        res.push_back(pair<uint16_t, uint16_t>(i, j));
                    }
                }
            }
        }
        return res;
    };

public:
    BrutePwHammDistAlgorithm(ExperimentParams &xParams): PwHammDistAlgorithm(xParams), seqInULLs(xParams.bytesPerSequence / 8) {
        if (xParams.bytesPerSequence % xParams.ALINGMENT_IN_BYTES != 0) {
            fprintf(stderr, "ERROR: brute algorithm does not support unaligned data.\n");
            exit(EXIT_FAILURE);
        }
    };

    virtual ~BrutePwHammDistAlgorithm() {
        if (pivotDist)
            delete[] pivotDist;
        if (seqAugmention)
            delete[] seqAugmention;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences) {
        preprocessing(sequences);
        vector<pair<uint16_t, uint16_t>> res;
        if (xParams.groupedBruteMode)
            res = xParams.pivotsFilterMode?findSimilarSequencesUsingGroupedApproach<true>()
                    :findSimilarSequencesUsingGroupedApproach<false>();
        else if (xParams.pivotsFilterMode)
            res = findSimilarSequencesUsingControlPivot();
        else
            res = findSimilarSequencesUsingStandardBrute();
        postprocessing();
        return res;
    }

    vector<pair<uint16_t, uint16_t>> findSimilarSequences(const uint8_t* sequences,
                                                          const vector<pair<uint16_t, uint16_t>> pairs) {
        preprocessing(sequences);
        vector<pair<uint16_t, uint16_t>> res;
        if (xParams.pivotsFilterMode)
            for (pair<uint16_t, uint16_t> pair: pairs) {
                if (testSequencesSimilarity<true>(pair.first, pair.second))
                    res.push_back(pair);
            }
        else
            for (pair<uint16_t, uint16_t> pair: pairs) {
                if (testSequencesSimilarity<false>(pair.first, pair.second))
                    res.push_back(pair);
            }
        postprocessing();
        return res;
    }

    template<bool usePivotFilter>
    inline bool testSequencesSimilarity(uint16_t i, uint16_t j) {
        if(usePivotFilter) {
            filterResult res = pivotFilter(i, j);
            if (res == different)
                return false;
            else if (res == similar)
                return true;
        }
        uint8_t* x = seq1 + (size_t) i * xParams.bytesPerSequence;
        uint8_t* y = seq2 + (size_t) j * xParams.bytesPerSequence;
        return testSequencesSimilarity(x, y);
    };

    inline bool testSequencesSimilarity(const void* seq1, const void* seq2) {
        return findSequencesDistance<true>(seq1, seq2) <= xParams.k;
    }

    string getName() {
        return (compact?COMPACT_PREFIX_ID:"") + (xParams.interleaveBitsMode?(
                    xParams.lazyInterleaveBitsMode?LAZY_INTERLEAVE_BITS_PREFIX_ID:INTERLEAVE_BITS_PREFIX_ID):"") +
            (xParams.pivotsFilterMode?(
                    xParams.pivotsElectionMode?ELECTION_PIVOT_FILTER_PREFIX_ID:PIVOT_FILTER_PREFIX_ID):"") +
            (xParams.groupedBruteMode?GROUPED_PREFIX_ID:"") +
            (shortcircuit?"":NO_SHORT_CIRCUIT_PREFIX_ID) + BRUTE_FORCE_ID +
            (binaryAlphabet?BINARY_MODE_ID_SUFFIX:"");
    }
};

#endif //ALGORITHM_BASE_H
