#ifndef __POLAR_H_
#define __POLAR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <inttypes.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <bitset>
#include <assert.h>
#include <chrono>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

// Define
#define CLLR 999 // Load channel LLR
#define PF 0 // Permutation function
#define FF 1 // F-function
#define PS 2 // partial sum update
#define GF 3 // G-function
#define SF 4 // Sign function
#define LF SF
#define R0 6
#define R1 7
#define REP 8
#define SPC 9
#define RR 10
#define HDM 11

#define BASELINE 0
#define NEW_ALGO 1
#define ALGO BASELINE

#define ZERO_CODEWORD 0

#define SIGNED true
#define UNSIGNED false

#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

// Structs
struct SCschedule
{
    vector<unsigned int> fspg_location;
    vector<unsigned int> fspg_stage;
    vector<vector<unsigned int>> fspg_branch;
};

struct FunctionType
{
    uint32_t name;
    uint32_t stage;
    uint32_t upBranch;
    uint32_t lowBranch;
};

struct SSC_Intruction
{
    uint32_t FID; //Function ID
    bool isChannelLLR;
    
    // alpha-router related control signals
    uint32_t in_LLR_idx_s;  // input LLR starting index needed by the decoder
    uint32_t in_LLR_idx_e;  // input LLR ending index needed by the decoder
    uint32_t out_LLR_idx_s; // output LLR starting index needed by the decoder
    uint32_t out_LLR_idx_e; // output LLR ending index needed by the decoder

    uint32_t in_beta_idx_s; // input beta memory starting index needed by the decoder
    uint32_t in_beta_idx_e; // input beta memory ending index needed by the decoder
    uint32_t out_beta_idx_s; // output beta memory starting index needed by the decoder
    uint32_t out_beta_idx_e; // output beta memory ending index needed by the decoder
};

struct StackPath
{
    vector<vector<double>> LLR;
    vector<vector<uint8_t>> pSum;
    double PM;
    uint32_t extended_cnt;
    uint32_t last_step;
    uint32_t last_info_idx;
    bool isFrozen;
    bool isDone;
};

// Marcos
inline uint8_t DIFF_BIT(uint8_t a, uint8_t b)
{
    if ((a - b) >= (b - a))
        return (a - b);
    else
        return (b - a);
}

inline uint8_t HARD_DECISION(double a)
{
    return (a > 0) ? 0 : 1;
}

template <typename T>
inline T sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T>
inline T clip(const T &n, const T &lower, const T &upper)
{
    return std::max(lower, std::min(n, upper));
}

inline void PolarCodeFFunction(double a, double b, double &c)
{
    c = sgn(a) * sgn(b) * min(abs(a), abs(b));
}

inline void PolarCodeGFunction(double a, double b, uint8_t s, double &c)
{
    c = b + (1 - 2 * s) * a;
}

inline void Int8_FFunction(int8_t a, int8_t b, int8_t &c)
{
    int8_t sign_val = ((a ^ b) & 0x80) >> 7;
    int8_t abs_a = (0x01 - ((a & 0x80) >> 6)) * a;
    int8_t abs_b = (0x01 - ((b & 0x80) >> 6)) * b;
    c = (1 - 2 * sign_val) * min(abs_a, abs_b);
}

inline void Int8_GFunction(int8_t a, int8_t b, uint8_t s, int8_t &c)
{
    c = b + (1 - 2 * s) * a;
}

inline double SIGMOID(double a)
{
    return 1 / (1 + exp(-a));
}

inline uint32_t nChoosek(unsigned n, unsigned k)
{
    if (k > n)
        return 0;
    if (k * 2 > n)
        k = n - k;
    if (k == 0)
        return 1;

    int result = n;
    for (int i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

inline double exp_approx(double x)
{
    return 1.0 + x + x * x / 2.0 + pow(x, 3) / 6.0 + pow(x, 4) / 24.0;
}

template <typename T>
inline void sort_indexes(const vector<T> &v, vector<uint32_t> &idx)
{
    // initialize original index locations
    //vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    sort(idx.begin(), idx.end(),
         [&v](uint32_t i1, uint32_t i2) { return v[i1] < v[i2]; });
}

template <typename T>
inline void sort_indexes_dec(const vector<T> &v, vector<uint32_t> &idx)
{
    // initialize original index locations
    //vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    sort(idx.begin(), idx.end(),
         [&v](uint32_t i1, uint32_t i2) { return v[i1] > v[i2]; });
}

template <typename T, typename A>
uint32_t arg_max(std::vector<T, A> const &vec)
{
    return static_cast<uint32_t>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename A>
uint32_t arg_min(std::vector<T, A> const &vec)
{
    return static_cast<uint32_t>(std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

extern double SPA(double a, double b);
extern double MIN_SUM(double a, double b);
extern double RPA_PROJECT(double a, double b);

extern void GetNodeType(vector<uint8_t> FrzVec, uint32_t &NodeType);
extern int GetNodeType(uint32_t N, uint32_t log2N, vector<uint8_t> FrzVec,
                       uint32_t begin_idx, uint32_t end_idx);
extern uint32_t GetNodeType(uint32_t r, uint32_t m);
extern void QUAN(double &x, double _2_to_m_, double min_quan, double max_quan);
extern void QUAN(double &x, uint64_t n, uint64_t m, bool s);
extern void sort_indexes_dec(vector<double> &v, vector<uint32_t> &idx);
extern void sort_indexes_double(vector<double> &v, vector<uint32_t> &idx);
extern double BitErrProb(vector<double> &Alpha, vector<uint32_t> &NonFrzIdx, uint32_t InfoIdx);
extern void PolarCodeEncoding(vector<uint8_t> &u, vector<uint8_t> &x);
extern void PC_PSumUdate(vector<uint8_t> &u, vector<uint8_t> &x);

extern void dec2bin(uint32_t a, vector<uint8_t> &b);
extern void bin2dec(uint32_t &a, vector<uint8_t> b);
extern uint32_t n_choose_k(uint32_t n, uint32_t k);

extern uint32_t vecdiff(vector<uint8_t> &a , vector<uint8_t> &b);

/* The polar code structure holds constant parameters of the code 
and all decoding algorithm*/
class PC
{
private:
    /* Functions */
    void Get_RelVec(void);
    void Get_CRCPoly(void);

public:
    /* Data structure */
    string code;    // type of code
    uint32_t SimID; // Simulation ID
    uint32_t N;     // code length
    uint32_t log2N; // code length
    uint32_t K;     // message length
    uint32_t C;     // CRC length
    uint32_t CS;
    uint32_t CE;    // THe number of channel errors
    uint32_t K_C;
    uint32_t L; // list size
    uint32_t S; // stack size
    uint32_t log2L;
    uint32_t I_max;      // number of BP iterations
    uint32_t I_min;      // number of BP iterations before early-termination
    uint32_t T_max;      // number of flipping attemps for bit-fliping decoding
    uint32_t P_max;      // number of random permutations
    uint32_t S_max;      // the index of the maximum fixed stage
    uint32_t k;          // the number of bandit arms
    uint32_t curSRN;     // The current value of the SNR times 10
    uint64_t Attempt, Latency, Complexity;
    double Memory;
    double LLR_Max;      // Maximum absolute LLR value
    bool isZeroCW;       // Use all-zero
    bool isTraining;     // Training flag of the bandit decoders
    bool isTrainDone;    // Training flag of the bandit decoders
    bool isFirstFrame;   // First frame flag
    bool isCS;           // Use critical set or not
    bool isGenie;        // Use a genie to pickup the correct codeword
    bool isEvalX;        // Validation using the codeword
    string Decoder;      // decoder
    string policy;       // rl policy
    string PMod;         // Power constrain mode
    SCschedule schedule; // decoding scheduling based on SC
    vector<FunctionType> RLDA_FunctionVec;
    vector<FunctionType> FastSC_FHT_FunctionVec;
    vector<FunctionType> FastFunctionVec;
    vector<FunctionType> OriginFunctionVec;
    vector<uint32_t> Rel;             // reliability vector
    vector<vector<uint32_t>> RelList; // list of reliability vector
    vector<uint8_t> Frz;              // frozen-bit flag
    vector<uint32_t> NonFrzIdx;       // non-frozen bit indices
    vector<uint32_t> NtoK_LUT;        // non-frozen bit indices
    vector<uint32_t> FrzIdx;          // frozen bit indices
    vector<uint8_t> CRCPoly;          // CRC polynomial
    vector<vector<uint8_t>> CRC_H;    // CRC parity-check matrix
    vector<vector<uint32_t>> CV_MAP;  // Check-to-Variable index map
    vector<vector<uint32_t>> VC_MAP;  // Variable-to-Check index map
    vector<uint8_t> m;                // information bits
    vector<uint8_t> crc;              // CRC bits
    bool syndrome;                    // Syndrome
    vector<uint8_t> u;                // message word
    vector<uint8_t> x;                // code word
    vector<vector<uint8_t>> pSum;     // correct partial sum
    uint8_t crc_sum;                  // CRC sum
    uint32_t RM_r;
    uint32_t RM_m;
    vector<vector<vector<uint32_t>>> QSpace;


    double beta;
    double beta_init;
    double r_avg;
    double lr;
    double dQt;
    uint32_t B;
    uint64_t TrainCnt;

    //------- Decoding time steps-------------
    uint32_t T_SCL, T_SSCL, T_FSSCL = 0;

    //------- Parameterized decoders ---------
    bool isGenData = false;

    //------- Quantization parameters --------
    uint8_t quan = 0;
    uint8_t LLR_n = 3;
    uint8_t LLR_m = 2;

    /* Functions */
    PC(void);
    void PC_Init(void);
    void SCScheduling(void);
    void RM_Form_QSpace(void);
    //void RPA_Form_QSpace(vector<uint32_t> &y_idx, uint32_t m, uint32_t i, vector<uint32_t> &y_b_idx);

    vector<FunctionType> FastSC_FHT_Scheduling(uint32_t offsetBranch,
                                               vector<uint8_t> &CurFrz);

    vector<FunctionType> RLDA_Scheduling(uint32_t offsetBranch, vector<uint8_t> &CurFrz);

    vector<FunctionType> FastSCScheduling(uint32_t offsetBranch,
                                          vector<uint8_t> &CurFrz);

    vector<FunctionType> OriginSCScheduling(uint32_t offsetBranch,
                                            vector<uint8_t> &CurFrz);

    void CRC_Gen(void);
    void CRC_Check(vector<uint8_t> &u_hat);
    void Get_CRC_H_Mat(void);
    void PC_Encode(default_random_engine &engine,
                   uniform_real_distribution<double> &uni_dist);
    void PC_HD_Decode(vector<double> &y, vector<uint8_t> &u_hat);
    bool PC_SyndromeCheck(vector<uint8_t> x);
};

#endif