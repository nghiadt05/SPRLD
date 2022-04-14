#ifndef _PER_H_
#define _PER_H_

#include "polar.h"
#include <bitset>

extern uint32_t factorial(uint32_t n);

extern vector<vector<uint32_t>> FormPerN(vector<uint32_t> inPer);

extern void BitIdxMapping(vector<uint32_t> &PerIdxVec,
                          vector<uint32_t> &InvPerIdxVec,
                          vector<uint32_t> PerStage,
                          uint32_t N, uint32_t log2N);

extern void BitIdxMapping(vector<uint32_t> &PerIdxVec,
                          vector<uint32_t> PerStage,
                          uint32_t N, uint32_t log2N);

extern void FormPermutations(vector<uint32_t> FixedStages,
                             vector<uint32_t> PerStages,
                             uint32_t log2N,
                             uint32_t S_max,
                             vector<vector<uint32_t>> &StagePerSet);

extern void PerVec(vector<uint32_t> perIdxVec, vector<double> inVec, vector<double> &outVec);
extern void PerVec(vector<uint32_t> perIdxVec, vector<uint8_t> inVec, vector<uint8_t> &outVec);

#endif