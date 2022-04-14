#ifndef _SCL_H_
#define _SCL_H_

#include "polar.h"

#define SPC_FLIP 1
#define GEN_SCLF_DATASET 0
#define HWF 1

struct ErrHist
{
	bool FoundFirstErr = false;
	uint32_t FirstErrStep = 0;
	uint32_t FirstErrInd = 0;
};

struct FirstErr
{
	uint32_t PathIdx = 0;
	uint32_t FirstErrStep = 0;
	uint32_t FirstErrInd = 0;
};

struct ErrPos
{
	double FlipMetric = (double)UINT32_MAX;
	uint32_t PathIdx = UINT32_MAX;
	uint32_t FirstErrStep = UINT32_MAX;
	uint32_t FirstErrInd = UINT32_MAX;
};

extern void FHT(vector<double> &y, uint32_t &stage, double &LLR_Max);
extern void FHT(vector<vector<uint32_t>> &Rel, vector<double> &y, uint32_t &stage,
				bool isReverseX, vector<uint8_t> &u_hat, vector<uint8_t> &x_hat);
extern void FHTL(vector<vector<uint32_t>> &Rel, vector<double> &y, uint32_t stage, uint32_t L_out,
				 bool isFHTDone, double abs_LLR, vector<vector<uint8_t>> &u_hat,
				 vector<vector<uint8_t>> &x_hat, vector<double> &dPM);
extern void RLDA(PC &myPC, vector<double> y, vector<uint8_t> &u_hat);
extern void SPRLD(PC &myPC, vector<double> y, vector<uint8_t> &u_hat, double &bestPM);
extern void Ens_SPRLD(PC &myPC, vector<double> y, vector<uint8_t> &u_hat);
extern void SSC_FHT(PC &myPC, vector<double> &y, vector<uint8_t> &u_hat);
extern void Aut_SSC_FHT(PC &myPC, vector<double> y, vector<uint8_t> &u_hat);

#endif