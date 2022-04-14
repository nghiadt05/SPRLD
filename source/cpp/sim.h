#ifndef __SIM_H_
#define __SIM_H_

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <iostream>

using namespace std;

class SIM
{
public:
    string Channel; // Channel modelling
    string Noise;   // Noise modelling
    string Mod;     // Modulation technique

    uint32_t MAX_FRAMES;    // Maximum number of simulation frames
    uint32_t MIN_FRAMES;    // Minimum number of simulation frames
    uint32_t MIN_ERRORS;    // Minimum number of frame errors
    uint32_t STABLE_ERRORS; // Stable errors if num. frames exceed MAX_FRAMES
    uint32_t FC;            // Frame counts
    uint32_t FEC;           // Frame error counts
    uint32_t BEC;           // Bit error counts
    uint32_t TC;            // Total counts of decoding attempts
    uint64_t LatencyC;      // Total counts of decoding attempts
    uint64_t ComplexityC;   // Total counts of decoding attempts
    uint64_t TimeSteps;     // Decoding time steps
    uint32_t RandSeed;      // Random seed for simulation
    uint32_t SimID;         // Simulation ID for parallel simulation
    uint32_t UFrame;        // Frames to update

    double SNR_START;  // Start SNR for simulation
    double SNR_END;    // Stop SNR for simulation
    double SNR_STEP;   // SNR step for simulation
    double TARGET_FER; // Target FER
    double Millis;     // Total counts of time in milliseconds
    double FER;        // Frame error rate
    double BER;        // Bit error rate
    double TP;         // Decoding through puts
    double AvgT;       // Average number of decoding attempts
    double AvgLatency;
    double AvgComplexity;
    double AvgTimeSteps; // Average number of decoding time steps
    double AvgTransCnt;  // Average number of transcendental operations performed
    double AvgAddCnt;    // Average number of additions performed
    double TransCnt;     // Total number of transcendental operations performed
    double AddCnt;       // Total number of additions performed
    bool isChannelLLR;   // Use channel output or channel LLR for decoding

    SIM(void);
};

#endif