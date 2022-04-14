#include "per.h"

uint32_t factorial(uint32_t n)
{
    // single line to find factorial
    return (n == 1 || n == 0) ? 1 : n * factorial(n - 1);
}

vector<vector<uint32_t>> FormPerN(vector<uint32_t> inPer)
{
    uint32_t inPerSize = inPer.size();
    vector<vector<uint32_t>> outPer(factorial(inPerSize),
                                    vector<uint32_t>(inPerSize));
    if (inPerSize == 2)
    {
        outPer[0][0] = inPer[0];
        outPer[0][1] = inPer[1];
        outPer[1][0] = inPer[1];
        outPer[1][1] = inPer[0];
    }
    else
    {
        for (auto i = 0; i < inPer.size(); i++)
        {
            // Form the current inPer
            vector<uint32_t> tmp_inPer;
            for (auto j = 0; j < inPerSize; j++)
                if (j != i)
                    tmp_inPer.push_back(inPer[j]);

            // Get all the permutations of tmp_inPer
            uint32_t factorial_inperSize_1 = factorial(inPerSize - 1);
            vector<vector<uint32_t>> tmp_outPer(factorial_inperSize_1,
                                                vector<uint32_t>(inPerSize - 1));
            tmp_outPer = FormPerN(tmp_inPer);

            // Concatenate the output to inPer[i]
            uint32_t curIdx = i * factorial_inperSize_1;
            for (auto j = 0; j < factorial_inperSize_1; j++)
            {
                outPer[curIdx + j][0] = inPer[i];
                for (auto k = 0; k < inPerSize - 1; k++)
                    outPer[curIdx + j][k + 1] = tmp_outPer[j][k];
            }
        }
    }
    return outPer;
}

void FormPermutations(vector<uint32_t> FixedStages,
                      vector<uint32_t> PerStages,
                      uint32_t log2N,
                      uint32_t S_max,
                      vector<vector<uint32_t>> &StagePerSet)
{
    // Get all the permutations of the permutable stages
    StagePerSet = FormPerN(PerStages);

    // Concatenate the fixed stages to the permutable stages
    if(FixedStages.size())
    {
        for (auto perIdx = 0; perIdx < StagePerSet.size(); perIdx++)
            StagePerSet[perIdx].insert(StagePerSet[perIdx].end(),
                                       FixedStages.begin(), FixedStages.end());
    }    
}

void BitIdxMapping(vector<uint32_t> &PerIdxVec,
                   vector<uint32_t> &InvPerIdxVec,
                   vector<uint32_t> PerStage,
                   uint32_t N, uint32_t log2N)
{
    // Form the stage-index permutation vector
    // to bit-index permutation vector
    vector<uint32_t> tmp(N);
    for (uint32_t idx = 0; idx < N; idx++)
    {
        vector<uint32_t> bin(log2N);
        vector<uint32_t> per_bin(log2N);
        uint32_t cur_idx = idx;
        uint32_t per_idx = 0;

        // Binary conversion
        for (uint32_t s = 0; s < log2N; s++)
        {
            if (cur_idx % 2)
                bin[s] = 1;
            cur_idx >>= 1;
        }

        // Permuting binary expression
        for (uint32_t s = 0; s < log2N; s++)
            per_bin[PerStage[s]] = bin[s];

        // Get the permuted indx
        for (uint32_t s = 0; s < log2N; s++)
            per_idx += per_bin[s] << s;

        PerIdxVec[idx] = per_idx;
        InvPerIdxVec[per_idx] = idx;
    }
}


void BitIdxMapping(vector<uint32_t> &PerIdxVec,
                   vector<uint32_t> PerStage,
                   uint32_t N, uint32_t log2N)
{
    // Form the stage-index permutation vector
    // to bit-index permutation vector
    PerIdxVec[0] = 0;
    PerIdxVec[N-1] = N-1;
    for (uint32_t idx = 1; idx < N-1; idx++)
    {
        // Binary conversion
        string str = std::bitset<16>(idx).to_string(); // string conversion
        str = str.substr(16 - log2N, log2N);
        string per_str = str;

        // Permuting binary expression
        for (uint32_t s = 0; s < log2N; s++)
            per_str[PerStage[s]] = str[s];

        // Get the permuted indx
        PerIdxVec[idx] = 0;
        for (uint32_t s = 0; s < log2N; s++)
            if (per_str[log2N - 1 - s] == '1')
                PerIdxVec[idx] += 1 << s;
    }
}

void PerVec(vector<uint32_t> perIdxVec, vector<double> inVec, vector<double> &outVec)
{
    for (auto i = 0; i < perIdxVec.size(); i++)
        outVec[perIdxVec[i]] = inVec[i];
}

void PerVec(vector<uint32_t> perIdxVec, vector<uint8_t> inVec, vector<uint8_t> &outVec)
{
    for (auto i = 0; i < perIdxVec.size(); i++)
        outVec[perIdxVec[i]] = inVec[i];
}