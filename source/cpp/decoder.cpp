#include "polar.h"
#include "decoder.h"

void FHT(vector<double> &y, uint32_t &stage, double &LLR_Max)
{
    // Perform the in-place FHT transform
    static double a, b;
    static uint32_t max_idx;
    static uint32_t s;
    max_idx = y.size() - 1;

    for (s = 0; s < stage; s++)
    {
        uint32_t seg_len = 1 << (stage - s);
        for (auto i = 0; i < (1 << s); i++)
        {
            uint32_t min_idx = i * seg_len;
            for (auto j = min_idx; j < min_idx + (seg_len >> 1); j++)
            {
                a = -y[j] + y[j + (seg_len >> 1)];
                b = +y[j] + y[j + (seg_len >> 1)];
                y[j] = a;
                y[j + (seg_len >> 1)] = b;
            }
        }
    }

    // Find the best [u0, u1,...,um]
    LLR_Max = -MAXFLOAT;
    uint32_t idx_max;
    for (auto i = 0; i < y.size(); i++)
    {
        if (LLR_Max < abs(y[i]))
        {
            LLR_Max = abs(y[i]);
            idx_max = i;
        }
    }
}

void FHT(vector<vector<uint32_t>> &Rel, vector<double> &y, uint32_t &stage,
         bool isReverseX, vector<uint8_t> &u_hat, vector<uint8_t> &x_hat)
{
    // Perform the in-place FHT transform
    static double a, b;
    static uint32_t max_idx;
    static uint32_t s;
    max_idx = y.size() - 1;

    for (s = 0; s < stage; s++)
    {
        uint32_t seg_len = 1 << (stage - s);
        for (auto i = 0; i < (1 << s); i++)
        {
            uint32_t min_idx = i * seg_len;
            for (auto j = min_idx; j < min_idx + (seg_len >> 1); j++)
            {
                a = -y[j] + y[j + (seg_len >> 1)];
                b = +y[j] + y[j + (seg_len >> 1)];
                y[j] = a;
                y[j + (seg_len >> 1)] = b;
            }
        }
    }

    // Find the best [u0, u1,...,um]
    static double LLR_max;
    static uint32_t idx_max;
    LLR_max = -MAXFLOAT;
    for (auto i = 0; i < y.size(); i++)
        if (LLR_max < abs(y[i]))
        {
            LLR_max = abs(y[i]);
            idx_max = i;
        }

    // Extract the reliability vector of the
    // sub RM(1,m) code from the mother RM code

    // Find u0 based on y[idx_max]
    if (y[idx_max] > 0)
        u_hat[Rel[stage][0]] = 0;
    else
        u_hat[Rel[stage][0]] = 1;

    // Find [u1,...,um]
    idx_max = max_idx - idx_max;
    for (auto s = 0; s < stage; s++)
    {
        u_hat[Rel[stage][s + 1]] = idx_max % 2;
        idx_max >>= 1;
    }

    // Get x_hat from u_hat
    x_hat = u_hat;
    for (uint32_t i = 0; i < stage; ++i)
        for (uint32_t j = 0; j < (y.size() >> (i + 1)); ++j)
            for (uint32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l)
                x_hat[l] ^= x_hat[l + (1 << i)];
}

void FHTL(vector<vector<uint32_t>> &Rel, vector<double> &y, uint32_t stage, uint32_t L_out,
          bool isFHTDone, double abs_LLR, vector<vector<uint8_t>> &u_hat,
          vector<vector<uint8_t>> &x_hat, vector<double> &dPM)
{
    // Perform the in-place FHT transform
    static uint32_t s;
    static double a, b;
    static uint32_t path_idx;
    vector<uint32_t> PLM_Idx(y.size());
    vector<double> org_y(y.begin(), y.end());
    vector<double> abs_y_transformed(y.begin(), y.end());

    if (!isFHTDone)
    {
        abs_LLR = 0.0;
        for (auto j = 0; j < y.size(); j++)
            abs_LLR += abs(org_y[j]);

        for (s = 0; s < stage; s++)
        {
            uint32_t seg_len = 1 << (stage - s);
            for (auto i = 0; i < (1 << s); i++)
            {
                uint32_t min_idx = i * seg_len;
                for (auto j = min_idx; j < min_idx + (seg_len >> 1); j++)
                {
                    a = y[j] + y[j + (seg_len >> 1)];
                    b = y[j] - y[j + (seg_len >> 1)];
                    y[j] = a;
                    y[j + (seg_len >> 1)] = b;
                }
            }
        }
    }

    for (auto j = 0; j < y.size(); j++)
        abs_y_transformed[j] = abs(y[j]);

    // Get the list of u_hat and x_hat corresponding to
    // the given channel LLR
    if (L_out < y.size())
        sort_indexes_dec(abs_y_transformed, PLM_Idx);

    for (auto l = 0; l < L_out; l++)
    {
        // Obtain the candidate path index
        if (L_out < y.size())
            path_idx = PLM_Idx[l];
        else
            path_idx = l;

        // Calculate the path metric of the current path by reusing the transformed y
		// Multiplication with 0.5 is right-shift by one bit, thus the complexity of this
		// computation is 1 unit.
        dPM[l] = 0.5 * (abs_LLR - abs_y_transformed[path_idx]);

        // Get the information bits
        u_hat[l][Rel[stage][0]] = y[path_idx] > 0 ? 0 : 1;
        for (auto s = 0; s < stage + 1; s++)
        {
            u_hat[l][Rel[stage][s + 1]] = path_idx % 2;
            path_idx >>= 1;
        }

        // Get x_hat from u_hat
        x_hat[l] = u_hat[l];
        for (uint32_t i = 0; i < stage; ++i)
            for (uint32_t j = 0; j < (y.size() >> (i + 1)); ++j)
                for (uint32_t k = (j << (i + 1)); k < (((j << 1) + 1) << i); ++k)
                    x_hat[l][k] ^= x_hat[l][k + (1 << i)];
    }
}

void RLDA(PC &myPC, vector<double> y, vector<uint8_t> &u_hat)
{
	static vector<vector<vector<uint32_t>>> GenPerIdxList;
	static vector<uint8_t> x_hat(myPC.N);

	static vector<uint32_t> orgPerIdx(myPC.L);
	static vector<vector<vector<double>>> LLR(myPC.L,
											  vector<vector<double>>(myPC.log2N + 1,
																	 vector<double>(myPC.N, 0)));
	static vector<vector<vector<uint8_t>>> Beta(myPC.L,
												vector<vector<uint8_t>>(myPC.log2N + 1,
																		vector<uint8_t>(myPC.N, 0)));
	static vector<double> PM(myPC.L << 1);
	static vector<vector<vector<uint32_t>>> perIdx(myPC.L, vector<vector<uint32_t>>(myPC.log2N + 1,
																					vector<uint32_t>(myPC.N)));

	static vector<uint32_t> tmp_orgPerIdx(myPC.L);
	static vector<vector<vector<double>>> tmp_LLR(myPC.L,
												  vector<vector<double>>(myPC.log2N + 1,
																		 vector<double>(myPC.N, 0)));
	static vector<vector<vector<uint8_t>>> tmp_Beta(myPC.L,
													vector<vector<uint8_t>>(myPC.log2N + 1,
																			vector<uint8_t>(myPC.N, 0)));
	static vector<double> tmp_PM(myPC.L);
	static vector<vector<vector<uint32_t>>> tmp_perIdx(myPC.L, vector<vector<uint32_t>>(myPC.log2N + 1,
																						vector<uint32_t>(myPC.N)));

	static uint32_t curL;
	static vector<uint32_t> PathIdx(myPC.L << 1);
	static vector<vector<uint32_t>> initial_perIndicesVec;

	static vector<uint32_t> left_child_vec;
	static vector<bool> isPerSelect;
	static bool init = false;
	bool isFHTDone = false;

	// Initialization
	if (!init)
	{
		init = true;

		// Read all the general permutations
		for (auto m = 0; m < myPC.log2N + 1; m++)
		{
			string line;
			string config_file = "../../config/rmcodes/m" + to_string(m) + ".per";
			ifstream file(config_file);
			vector<uint32_t> PerIdx(1 << m);
			vector<vector<uint32_t>> PerIdxList;

			if (file.is_open())
			{
				while (getline(file, line))
				{
					uint32_t char_cnt = 0;
					uint32_t rel = 0;
					char delimiter = ' ';
					for (uint32_t i = 0; i < (1 << m); i++)
					{
						sscanf(&line[char_cnt], "%d", &rel);
						PerIdx[i] = rel;
						// find the next reliability element position
						while (line[char_cnt] != delimiter)
							char_cnt++;
						char_cnt++;
					}
					PerIdxList.push_back(PerIdx);
				}
			}
			GenPerIdxList.push_back(PerIdxList);
		}

		// Memory in Kbytes
		if (myPC.L == 1)
			myPC.Memory = (2 * myPC.N) * 32.0 + myPC.N;
		else
			myPC.Memory = myPC.N * (myPC.L + 1) * 32.0 + 2 * myPC.N * myPC.L;
		myPC.Memory /= 8192.0;
	}

	// Step 0: start with L random permutations.

	// Note that this initialization allows for a fair comparison between the proposed SPRLD decoder
	// and that of RLDA in terms of memory consumption, where the same number of L initial permutations
	// is used. Furthermore, by using the general permutations, one in general does not need to prune duplicated permutaitons
	// as in the conventional RLDP algorithm, whose complexity is not accounted for in the paper.
	// Therefore, using L initial permutations for L decoding paths is a practical choice for RLDA.

	uint32_t FFidx = 0;
	std::fill(PM.begin(), PM.end(), 0);
	curL = myPC.L;
	while (initial_perIndicesVec.size() < myPC.L)
		initial_perIndicesVec.push_back(GenPerIdxList[myPC.log2N][rand() % GenPerIdxList[myPC.log2N].size()]);
	for (auto l = 0; l < myPC.L; l++)
	{
		orgPerIdx[l] = l;
		for (auto i = 0; i < myPC.N; i++)
			LLR[l][myPC.log2N][initial_perIndicesVec[l][i]] = y[i];
	}

	// Step 1+: decoding
	for (auto step = 0; step < myPC.RLDA_FunctionVec.size(); step++)
	{
		uint32_t Function = myPC.RLDA_FunctionVec[step].name;
		uint32_t StageIdx = myPC.RLDA_FunctionVec[step].stage;
		uint32_t UpBrach = myPC.RLDA_FunctionVec[step].upBranch;
		uint32_t DownBranch = myPC.RLDA_FunctionVec[step].lowBranch;
		uint32_t N_stage = DownBranch - UpBrach + 1;

		if (Function == FF)
		{
			// Normal F functions
			for (auto l = 0; l < curL; l++)
				for (auto i = UpBrach; i <= DownBranch; i++)
					PolarCodeFFunction(LLR[l][StageIdx + 1][i],
									   LLR[l][StageIdx + 1][i + N_stage],
									   LLR[l][StageIdx][i]);

			// F-functions costs
			myPC.Complexity += N_stage * curL;
			myPC.Latency += 1;

			// Increase the F function index
			FFidx++;
		}
		else if (Function == GF)
		{
			for (auto l = 0; l < curL; l++)
				for (auto i = UpBrach; i <= DownBranch; i++)
					PolarCodeGFunction(LLR[l][StageIdx + 1][i - N_stage],
									   LLR[l][StageIdx + 1][i],
									   Beta[l][StageIdx][i - N_stage],
									   LLR[l][StageIdx][i]);

			myPC.Complexity += N_stage * curL;
			myPC.Latency += 1;
		}
		else if (Function == REP)
		{
			for (auto l = 0; l < curL; l++)
			{
				// Calculate the LLR of the information bit
				LLR[l][0][DownBranch] = 0.0;
				for (auto i = UpBrach; i <= DownBranch; i++)
					LLR[l][0][DownBranch] += LLR[l][StageIdx][i];

				// Calculate the path metric of the current path
				// and the forked path
				PM[l + curL] = PM[l];
				for (auto i = UpBrach; i <= DownBranch; i++)
				{
					if (HARD_DECISION(LLR[l][0][DownBranch]) !=
						HARD_DECISION(LLR[l][StageIdx][i]))
						PM[l] += abs(LLR[l][StageIdx][i]);
					else
						PM[l + curL] += abs(LLR[l][StageIdx][i]);
				}
			}
			if (myPC.L > 1)
			{
				myPC.Complexity += (N_stage << 1) * curL;
				myPC.Latency += 1;
			}

			// Path killing
			vector<uint32_t> SurvPath;
			vector<uint32_t> DeadPath;
			sort_indexes(PM, PathIdx);
			if (myPC.L > 1)
			{
				myPC.Complexity += 2 * myPC.L * (1 + myPC.log2L);
				myPC.Latency += 1 + myPC.log2L;
			}

			for (auto l = 0; l < curL; l++)
			{
				if (PathIdx[l] >= curL)
					SurvPath.push_back(PathIdx[l]);
				if (PathIdx[l + curL] < curL)
					DeadPath.push_back(PathIdx[l + curL]);
			}

			for (auto l = 0; l < SurvPath.size(); l++)
			{
				uint32_t RootIdx = SurvPath[l] - curL;
				uint32_t ReplaceIdx = DeadPath[l];

				PM[ReplaceIdx] = PM[SurvPath[l]];
				LLR[ReplaceIdx] = LLR[RootIdx];
				Beta[ReplaceIdx] = Beta[RootIdx];
				orgPerIdx[ReplaceIdx] = orgPerIdx[RootIdx];

				LLR[ReplaceIdx][0][DownBranch] = -LLR[ReplaceIdx][0][DownBranch];
			}

			// Make hard decisions
			for (auto l = 0; l < curL; l++)
				for (auto i = UpBrach; i <= DownBranch; i++)
					Beta[l][StageIdx][i] = HARD_DECISION(LLR[l][0][DownBranch]);
		}
		else if (Function == SPC)
		{
			// Find the minimum abs-LLR bit index
			// and the parity check sum of all the
			// current paths
			vector<uint8_t> SPC_Sum(myPC.L, 0);
			vector<double> min_absLLR(myPC.L, 9999.0);
			vector<vector<double>> absLLR_list(myPC.L, vector<double>(N_stage));
			vector<vector<uint32_t>> index_list(myPC.L, vector<uint32_t>(N_stage));
			vector<vector<uint8_t>> cur_x(myPC.L, vector<uint8_t>(N_stage, 0));

			for (auto l = 0; l < curL; l++)
			{
				for (auto i = UpBrach; i <= DownBranch; i++)
				{
					absLLR_list[l][i - UpBrach] = abs(LLR[l][StageIdx][i]);
					cur_x[l][i - UpBrach] = HARD_DECISION(LLR[l][StageIdx][i]);
					SPC_Sum[l] ^= cur_x[l][i - UpBrach];
				}

				sort_indexes(absLLR_list[l], index_list[l]);
				min_absLLR[l] = absLLR_list[l][index_list[l][0]];

				// Penalize the frozen-bit estimation
				if (SPC_Sum[l])
					PM[l] += min_absLLR[l];
			}
			if (myPC.L > 1)
			{
				myPC.Complexity += N_stage * StageIdx * curL;
				myPC.Latency += StageIdx;
			}

			// Bit flipping decision for non-minimal abs-LLRs
			uint32_t min_path_splits = min(N_stage - 1, myPC.L - 1);
			for (auto i = 1; i <= min_path_splits; i++)
			{
				// Calculate the PM of the forked paths
				for (auto l = 0; l < curL; l++)
				{
					uint32_t BitIdx = index_list[l][i] + UpBrach;
					PM[l + curL] = PM[l] + abs(LLR[l][StageIdx][BitIdx]) +
								   (1 - 2 * SPC_Sum[l]) * min_absLLR[l];
				}
				if (myPC.L > 1)
				{
					myPC.Complexity += curL;
					myPC.Latency += 1;
				}

				// Path killing
				vector<uint32_t> SurvPath;
				vector<uint32_t> DeadPath;
				sort_indexes(PM, PathIdx);
				if (myPC.L > 1)
				{
					myPC.Complexity += 2 * myPC.L * (1 + myPC.log2L);
					myPC.Latency += 1 + myPC.log2L;
				}

				for (auto l = 0; l < curL; l++)
				{
					if (PathIdx[l] >= curL)
						SurvPath.push_back(PathIdx[l]);
					if (PathIdx[l + curL] < curL)
						DeadPath.push_back(PathIdx[l + curL]);
				}

				for (auto l = 0; l < SurvPath.size(); l++)
				{
					uint32_t RootIdx = SurvPath[l] - curL;
					uint32_t ReplaceIdx = DeadPath[l];

					PM[ReplaceIdx] = PM[SurvPath[l]];
					LLR[ReplaceIdx] = LLR[RootIdx];
					Beta[ReplaceIdx] = Beta[RootIdx];
					orgPerIdx[ReplaceIdx] = orgPerIdx[RootIdx];

					SPC_Sum[ReplaceIdx] = SPC_Sum[RootIdx];
					cur_x[ReplaceIdx] = cur_x[RootIdx];
					min_absLLR[ReplaceIdx] = min_absLLR[RootIdx];
					index_list[ReplaceIdx] = index_list[RootIdx];

					cur_x[ReplaceIdx][index_list[RootIdx][i]] ^= 1;
					SPC_Sum[ReplaceIdx] ^= SPC_FLIP;
				}
			}

			// Making hard decisions and make sure the even-parity
			// condision is satisfied
			for (auto l = 0; l < curL; l++)
			{
				cur_x[l][index_list[l][0]] = 0;
				for (auto i = 1; i < index_list[l].size(); i++)
					cur_x[l][index_list[l][0]] ^= cur_x[l][index_list[l][i]];

				// Assign the data
				for (auto i = UpBrach; i <= DownBranch; i++)
					Beta[l][StageIdx][i] = cur_x[l][i - UpBrach];
			}
		}
		else if (Function == PS)
		{
			vector<uint8_t> tmp_x(N_stage);
			for (auto l = 0; l < curL; l++)
			{
				for (auto i = UpBrach; i < UpBrach + (N_stage >> 1); i++)
					tmp_x[i - UpBrach] = Beta[l][StageIdx - 1][i] ^ Beta[l][StageIdx - 1][i + (N_stage >> 1)];
				for (auto i = UpBrach + (N_stage >> 1); i <= DownBranch; i++)
					tmp_x[i - UpBrach] = Beta[l][StageIdx - 1][i];
				for (auto i = UpBrach; i <= DownBranch; i++)
					Beta[l][StageIdx][i] = tmp_x[i - UpBrach];
			}
		}
	}

	// Selection of the best path, and re-permute to get the orignal codeword
	uint32_t best_metric_path = 0;
	for (auto l = 0; l < curL; l++)
		if (PM[l] < PM[best_metric_path])
			best_metric_path = l;
	myPC.Complexity += myPC.L;
	myPC.Latency += myPC.log2L;

	for (auto i = 0; i < myPC.N; i++)
		x_hat[i] = Beta[best_metric_path][myPC.log2N][initial_perIndicesVec[orgPerIdx[best_metric_path]][i]];
	PolarCodeEncoding(x_hat, u_hat);
}

void SPRLD(PC &myPC, vector<double> y, vector<uint8_t> &u_hat, double &bestPM)
{
	static vector<vector<vector<uint32_t>>> GenPerIdxList;
	static vector<uint8_t> x_hat(myPC.N);

	static vector<uint32_t> orgPerIdx(myPC.L);
	static vector<vector<vector<double>>> LLR(myPC.L,
											  vector<vector<double>>(myPC.log2N + 1,
																	 vector<double>(myPC.N, 0)));
	static vector<vector<vector<uint8_t>>> Beta(myPC.L,
												vector<vector<uint8_t>>(myPC.log2N + 1,
																		vector<uint8_t>(myPC.N, 0)));
	static vector<double> PM(myPC.L << 1);
	static vector<vector<vector<uint32_t>>> perIdx(myPC.L, vector<vector<uint32_t>>(myPC.log2N + 1,
																					vector<uint32_t>(myPC.N)));

	static vector<uint32_t> tmp_orgPerIdx(myPC.L);
	static vector<vector<vector<double>>> tmp_LLR(myPC.L,
												  vector<vector<double>>(myPC.log2N + 1,
																		 vector<double>(myPC.N, 0)));
	static vector<vector<vector<uint8_t>>> tmp_Beta(myPC.L,
													vector<vector<uint8_t>>(myPC.log2N + 1,
																			vector<uint8_t>(myPC.N, 0)));
	static vector<double> tmp_PM(myPC.L);
	static vector<double> abs_sum(myPC.L);
	static vector<vector<vector<uint32_t>>> tmp_perIdx(myPC.L, vector<vector<uint32_t>>(myPC.log2N + 1,
																						vector<uint32_t>(myPC.N)));

	static uint32_t curL;
	static vector<uint32_t> PathIdx(myPC.L << 1);
	static vector<vector<uint32_t>> initial_perIndicesVec;

	static vector<uint32_t> left_child_vec;
	static vector<bool> isPerSelect;
	static bool init = false;
	bool isFHTDone = false;

	// Initialization
	if (!init)
	{
		init = true;

		// Read all the general permutations
		for (auto m = 0; m < myPC.log2N + 1; m++)
		{
			string line;
			string config_file = "../../config/rmcodes/m" + to_string(m) + ".per";
			ifstream file(config_file);
			vector<uint32_t> PerIdx(1 << m);
			vector<vector<uint32_t>> PerIdxList;

			if (file.is_open())
			{
				while (getline(file, line))
				{
					uint32_t char_cnt = 0;
					uint32_t rel = 0;
					char delimiter = ' ';
					for (uint32_t i = 0; i < (1 << m); i++)
					{
						sscanf(&line[char_cnt], "%d", &rel);
						PerIdx[i] = rel;
						// find the next reliability element position
						while (line[char_cnt] != delimiter)
							char_cnt++;
						char_cnt++;
					}
					PerIdxList.push_back(PerIdx);
				}
			}
			GenPerIdxList.push_back(PerIdxList);
		}

		// Get the left child types
		for (auto step = 0; step < myPC.FastSC_FHT_FunctionVec.size(); step++)
		{
			uint32_t Function = myPC.FastSC_FHT_FunctionVec[step].name;
			uint32_t StageIdx = myPC.FastSC_FHT_FunctionVec[step].stage;
			uint32_t UpBrach = myPC.FastSC_FHT_FunctionVec[step].upBranch;
			uint32_t DownBranch = myPC.FastSC_FHT_FunctionVec[step].lowBranch;
			uint32_t N_stage = DownBranch - UpBrach + 1;

			if (Function == FF) // ok
			{
				static uint32_t log2_N_p_stage = 0;
				static uint32_t left_child;
				vector<uint8_t> LeftFrzVec;

				log2_N_p_stage = 0;
				while ((1 << log2_N_p_stage) < (N_stage << 1))
					log2_N_p_stage++;

				// Get the type of the left and right nodes
				for (auto i = 0; i < N_stage; i++)
					LeftFrzVec.push_back(myPC.Frz[UpBrach + i]);
				GetNodeType(LeftFrzVec, left_child);
				left_child_vec.push_back(left_child);
				isPerSelect.push_back(false);
			}
		}

		// Turn off factor graph selection for the first S(S>1) F-functions
		for (auto i = 0; i < min(myPC.S, (uint32_t)isPerSelect.size()); i++)
			isPerSelect[i] = true;

		// Memory consumption in kilobytes
		if (myPC.policy != "par" && myPC.policy != "seq")
			myPC.policy = "par";

		if (myPC.policy == "par")
			if (myPC.L > 1)
				myPC.Memory = myPC.N * (myPC.log2N * myPC.L + 1) * 32.0 + myPC.log2N * 32.0 + 2 * myPC.N * myPC.L;
			else
				myPC.Memory = (myPC.log2N + 1) * myPC.N * 32.0 + myPC.log2N * 32.0 + myPC.N;
		else if (myPC.policy == "seq")
			if (myPC.L > 1)
				myPC.Memory = myPC.N * (myPC.L + 1) * 32.0 + myPC.log2N * 32.0 + 2 * myPC.N * myPC.L;
			else
				myPC.Memory = 2 * myPC.N * 32.0 + myPC.log2N * 32.0 + myPC.N;
		myPC.Memory /= 8192.0;
	}

	// Step 0: start with L random permutations.
	uint32_t FFidx = 0;
	std::fill(PM.begin(), PM.end(), 0);
	curL = myPC.L;
	while (initial_perIndicesVec.size() < myPC.L)
		initial_perIndicesVec.push_back(GenPerIdxList[myPC.log2N][rand() % GenPerIdxList[myPC.log2N].size()]);
	for (auto l = 0; l < myPC.L; l++)
	{
		orgPerIdx[l] = l;
		for (auto i = 0; i < myPC.N; i++)
			LLR[l][myPC.log2N][initial_perIndicesVec[l][i]] = y[i];
	}

	// Step 1+: decoding
	for (auto step = 0; step < myPC.FastSC_FHT_FunctionVec.size(); step++)
	{
		uint32_t Function = myPC.FastSC_FHT_FunctionVec[step].name;
		uint32_t StageIdx = myPC.FastSC_FHT_FunctionVec[step].stage;
		uint32_t UpBrach = myPC.FastSC_FHT_FunctionVec[step].upBranch;
		uint32_t DownBranch = myPC.FastSC_FHT_FunctionVec[step].lowBranch;
		uint32_t N_stage = DownBranch - UpBrach + 1;

		if (Function == FF)
		{
			if (isPerSelect[FFidx])
			{
				// Get the type of the left child node
				static uint32_t log2_N_p_stage = 0;
				static uint32_t left_child;

				log2_N_p_stage = 0;
				while ((1 << log2_N_p_stage) < (N_stage << 1))
					log2_N_p_stage++;

				left_child = left_child_vec[FFidx];

				// Find the best permutations
				vector<vector<uint32_t>> BitIdxPer(log2_N_p_stage, vector<uint32_t>(N_stage << 1));

				// Select random general permutations
				for (auto i = 0; i < log2_N_p_stage; i++)
					BitIdxPer[i] = GenPerIdxList[log2_N_p_stage][rand() % GenPerIdxList[log2_N_p_stage].size()];

				for (auto l = 0; l < curL; l++)
				{
					double PerMetric = 0.0;
					double BestPerMetric = -MAXFLOAT;
					double tmp_abs_sum;
					uint32_t BestPerIdx;
					vector<double> tmpLLRvec(N_stage << 1);
					vector<double> alpha(N_stage), best_alpha(N_stage);

					for (auto i = 0; i < log2_N_p_stage; i++)
					{
						for (auto ii = UpBrach; ii < UpBrach + (N_stage << 1); ii++)
							tmpLLRvec[BitIdxPer[i][ii - UpBrach]] = LLR[l][StageIdx + 1][ii];

						tmp_abs_sum = 0.0;
						for (auto ii = 0; ii < N_stage; ii++)
						{
							PolarCodeFFunction(tmpLLRvec[ii], tmpLLRvec[ii + N_stage], alpha[ii]);
							tmp_abs_sum += abs(alpha[ii]);
						}

						// F-function costs
						myPC.Complexity += N_stage;
						if (i == 0 && l == 0 && myPC.policy == "par")
							myPC.Latency += 1;
						else if (l == 0 && myPC.policy == "seq")
							myPC.Latency += 1;

						// Calculate the reliability metric based on different cases
						if (left_child == HDM)
						{
							vector<uint8_t> tmp_u_hat(alpha.size());
							vector<uint8_t> tmp_x_hat(alpha.size());
							FHT(alpha, StageIdx, PerMetric);
							// FHT costs
							myPC.Complexity += N_stage * log2_N_p_stage;
							if (i == 0 && l == 0 && myPC.policy == "par")
								myPC.Latency += log2_N_p_stage;
							else if (l == 0 && myPC.policy == "seq")
								myPC.Latency += log2_N_p_stage;
						}
						else // Rate-r
						{
							PerMetric = tmp_abs_sum;
							// Rate-1 permutation metric's costs
							myPC.Complexity += N_stage;
							if (i == 0 && l == 0 && myPC.policy == "par")
								myPC.Latency += 1;
							else if (l == 0 && myPC.policy == "seq")
								myPC.Latency += 1;
						}

						// Select the best permutation
						if (BestPerMetric < PerMetric)
						{
							BestPerIdx = i;
							best_alpha = alpha;
							abs_sum[l] = tmp_abs_sum;
							BestPerMetric = PerMetric;
						}

						// Complexity and latency of the graph selection
						myPC.Complexity += log2_N_p_stage;
						if (i == 0 && l == 0)
							myPC.Latency += (uint32_t)log2(log2_N_p_stage);
					}

					for (auto j = 0; j < N_stage << 1; j++)
						perIdx[l][StageIdx + 1][UpBrach + j] = BitIdxPer[BestPerIdx][j];

					for (auto j = 0; j < N_stage << 1; j++)
						tmpLLRvec[BitIdxPer[BestPerIdx][j]] = LLR[l][StageIdx + 1][UpBrach + j];

					for (auto j = 0; j < N_stage << 1; j++)
						LLR[l][StageIdx + 1][UpBrach + j] = tmpLLRvec[j];

					for (auto i = UpBrach; i <= DownBranch; i++)
						LLR[l][StageIdx][i] = best_alpha[i - UpBrach];
				}

				if (left_child == HDM)
					isFHTDone = true;
			}
			else
			{
				// No SP, just by passing indices
				for (auto l = 0; l < curL; l++)
					for (auto j = 0; j < N_stage << 1; j++)
						perIdx[l][StageIdx + 1][UpBrach + j] = j;

				// Normal F functions
				for (auto l = 0; l < curL; l++)
					for (auto i = UpBrach; i <= DownBranch; i++)
						PolarCodeFFunction(LLR[l][StageIdx + 1][i],
										   LLR[l][StageIdx + 1][i + N_stage],
										   LLR[l][StageIdx][i]);

				// F-functions costs
				myPC.Complexity += N_stage * curL;
				myPC.Latency += 1;
			}

			// Increase the F function index
			FFidx++;
		}
		else if (Function == GF)
		{
			for (auto l = 0; l < curL; l++)
				for (auto i = UpBrach; i <= DownBranch; i++)
					PolarCodeGFunction(LLR[l][StageIdx + 1][i - N_stage],
									   LLR[l][StageIdx + 1][i],
									   Beta[l][StageIdx][i - N_stage],
									   LLR[l][StageIdx][i]);
			myPC.Complexity += N_stage * curL;
			myPC.Latency += 1;
		}
		else if (Function == SPC)
		{
			// Find the minimum abs-LLR bit index
			// and the parity check sum of all the
			// current paths
			vector<uint8_t> SPC_Sum(myPC.L, 0);
			vector<double> min_absLLR(myPC.L, 9999.0);
			vector<vector<double>> absLLR_list(myPC.L, vector<double>(N_stage));
			vector<vector<uint32_t>> index_list(myPC.L, vector<uint32_t>(N_stage));
			vector<vector<uint8_t>> cur_x(myPC.L, vector<uint8_t>(N_stage, 0));

			for (auto l = 0; l < curL; l++)
			{
				for (auto i = UpBrach; i <= DownBranch; i++)
				{
					absLLR_list[l][i - UpBrach] = abs(LLR[l][StageIdx][i]);
					cur_x[l][i - UpBrach] = HARD_DECISION(LLR[l][StageIdx][i]);
					SPC_Sum[l] ^= cur_x[l][i - UpBrach];
				}

				sort_indexes(absLLR_list[l], index_list[l]);
				min_absLLR[l] = absLLR_list[l][index_list[l][0]];

				// Penalize the frozen-bit estimation
				if (SPC_Sum[l])
					PM[l] += min_absLLR[l];
			}
			if (myPC.L > 1)
				myPC.Complexity += N_stage * StageIdx * curL;
			else
				myPC.Complexity += N_stage * curL;
			myPC.Latency += StageIdx;

			// Bit flipping decision for non-minimal abs-LLRs
			uint32_t min_path_splits = min(N_stage - 1, myPC.L - 1);
			for (auto i = 1; i <= min_path_splits; i++)
			{
				// Calculate the PM of the forked paths
				for (auto l = 0; l < curL; l++)
				{
					uint32_t BitIdx = index_list[l][i] + UpBrach;
					PM[l + curL] = PM[l] + abs(LLR[l][StageIdx][BitIdx]) +
								   (1 - 2 * SPC_Sum[l]) * min_absLLR[l];
				}				
				if(myPC.L>1)
				{
					myPC.Complexity += curL;
					myPC.Latency += 1;
				}

				// Path killing
				vector<uint32_t> SurvPath;
				vector<uint32_t> DeadPath;
				sort_indexes(PM, PathIdx);
				if (myPC.L > 1)
				{
					myPC.Complexity += 2 * myPC.L * (1 + myPC.log2L);
					myPC.Latency += 1 + myPC.log2L;
				}

				for (auto l = 0; l < curL; l++)
				{
					if (PathIdx[l] >= curL)
						SurvPath.push_back(PathIdx[l]);
					if (PathIdx[l + curL] < curL)
						DeadPath.push_back(PathIdx[l + curL]);
				}

				for (auto l = 0; l < SurvPath.size(); l++)
				{
					uint32_t RootIdx = SurvPath[l] - curL;
					uint32_t ReplaceIdx = DeadPath[l];

					PM[ReplaceIdx] = PM[SurvPath[l]];
					LLR[ReplaceIdx] = LLR[RootIdx];
					Beta[ReplaceIdx] = Beta[RootIdx];
					perIdx[ReplaceIdx] = perIdx[RootIdx];
					orgPerIdx[ReplaceIdx] = orgPerIdx[RootIdx];

					SPC_Sum[ReplaceIdx] = SPC_Sum[RootIdx];
					cur_x[ReplaceIdx] = cur_x[RootIdx];
					min_absLLR[ReplaceIdx] = min_absLLR[RootIdx];
					index_list[ReplaceIdx] = index_list[RootIdx];

					cur_x[ReplaceIdx][index_list[RootIdx][i]] ^= 1;
					SPC_Sum[ReplaceIdx] ^= SPC_FLIP;
				}
			}

			// Making hard decisions and make sure the even-parity
			// condision is satisfied
			for (auto l = 0; l < curL; l++)
			{
				cur_x[l][index_list[l][0]] = 0;
				for (auto i = 1; i < index_list[l].size(); i++)
					cur_x[l][index_list[l][0]] ^= cur_x[l][index_list[l][i]];

				// Assign the data
				for (auto i = UpBrach; i <= DownBranch; i++)
					if (myPC.K == myPC.N - 1)
					{
						uint32_t idx = perIdx[l][StageIdx][i - UpBrach];
						Beta[l][StageIdx][i] = cur_x[l][idx];
					}
					else
						Beta[l][StageIdx][i] = cur_x[l][i - UpBrach];
			}
		}
		else if (Function == HDM) // ok
		{
			// Perform FHT decoding on the current L active paths
			uint32_t L_out = min(myPC.L, N_stage);
			vector<vector<vector<uint8_t>>> tmp_u_hat(curL, vector<vector<uint8_t>>(L_out, vector<uint8_t>(N_stage, 0)));
			vector<vector<vector<uint8_t>>> tmp_x_hat(curL, vector<vector<uint8_t>>(L_out, vector<uint8_t>(N_stage, 0)));
			vector<double> tmp_y(N_stage);
			vector<double> dPM(L_out);
			vector<double> newPM(curL * L_out);
			vector<uint32_t> newPM_idx(curL * L_out);
			vector<vector<uint8_t>> new_u_hat_list(curL * L_out, vector<uint8_t>(myPC.N));
			uint32_t newPathCnt = min(myPC.L, curL * L_out);

			for (auto l = 0; l < curL; l++)
			{
				// Prepare the channel LLR for the sub punctured Hadamard (HDM) code
				for (auto i = UpBrach; i <= DownBranch; i++)
					tmp_y[DownBranch - i] = LLR[l][StageIdx][i];

				// Get the best L_out paths originated from the current path
				FHTL(myPC.RelList, tmp_y, StageIdx, L_out, isFHTDone, abs_sum[l], tmp_u_hat[l], tmp_x_hat[l], dPM);
				if (!isFHTDone) // FHT complexity
					myPC.Complexity += N_stage * StageIdx;
				if (L_out > 1) // PM complexity
					myPC.Complexity += L_out;
				if (L_out == 1) // Sorting complexity
					myPC.Complexity += N_stage;
				else if (L_out < y.size())
					myPC.Complexity += N_stage * StageIdx;

				// Calculate the corresponding PM of the L_out paths
				for (auto i = 0; i < L_out; i++)
					newPM[l * L_out + i] = PM[l] + dPM[i];
			}
			if (!isFHTDone) // FHT latency
				myPC.Latency += StageIdx;
			if (L_out < y.size()) // Comparison latency
				myPC.Latency += StageIdx;
			isFHTDone = false;

			// Sort all the path metrics of the FHTL decoders
			sort_indexes(newPM, newPM_idx);
			if (myPC.L > 1)
			{
				myPC.Complexity += newPM.size() * (uint32_t)log2(newPM.size());
				myPC.Latency += (uint32_t)log2(newPM.size());
			}

			for (auto l = 0; l < newPathCnt; l++)
			{
				uint32_t path_idx = newPM_idx[l];
				uint32_t origin_path_idx = (uint32_t)path_idx / L_out;
				uint32_t FHT_path_idx = (uint32_t)path_idx % L_out;

				tmp_PM[l] = newPM[path_idx];
				tmp_LLR[l] = LLR[origin_path_idx];
				tmp_Beta[l] = Beta[origin_path_idx];
				tmp_perIdx[l] = perIdx[origin_path_idx];
				tmp_orgPerIdx[l] = orgPerIdx[origin_path_idx];
				for (auto i = UpBrach; i <= DownBranch; i++)
					if (myPC.K == myPC.log2N + 1)
					{
						uint32_t idx = perIdx[origin_path_idx][StageIdx][i - UpBrach];
						tmp_Beta[l][StageIdx][i] = tmp_x_hat[origin_path_idx][FHT_path_idx][idx];
					}
					else
						tmp_Beta[l][StageIdx][i] = tmp_x_hat[origin_path_idx][FHT_path_idx][i - UpBrach];
			}

			for (auto l = 0; l < newPathCnt; l++)
			{
				PM[l] = tmp_PM[l];
				LLR[l] = tmp_LLR[l];
				Beta[l] = tmp_Beta[l];
				perIdx[l] = tmp_perIdx[l];
				orgPerIdx[l] = tmp_orgPerIdx[l];
			}
			curL = newPathCnt;
		}
		else if (Function == PS)
		{
			vector<uint8_t> tmp_x(N_stage);
			for (auto l = 0; l < curL; l++)
			{
				for (auto i = UpBrach; i < UpBrach + (N_stage >> 1); i++)
					tmp_x[i - UpBrach] = Beta[l][StageIdx - 1][i] ^ Beta[l][StageIdx - 1][i + (N_stage >> 1)];
				for (auto i = UpBrach + (N_stage >> 1); i <= DownBranch; i++)
					tmp_x[i - UpBrach] = Beta[l][StageIdx - 1][i];
				for (auto i = UpBrach; i <= DownBranch; i++)
					Beta[l][StageIdx][i] = tmp_x[perIdx[l][StageIdx][i]];
			}
		}
	}

	// Selection of the best path, and re-permute to get the orignal codeword
	uint32_t best_metric_path = 0;
	for (auto l = 0; l < curL; l++)
		if (PM[l] < PM[best_metric_path])
			best_metric_path = l;
	for (auto i = 0; i < myPC.N; i++)
		x_hat[i] = Beta[best_metric_path][myPC.log2N][initial_perIndicesVec[orgPerIdx[best_metric_path]][i]];
	bestPM = PM[best_metric_path];
	if (myPC.L > 1)
	{
		myPC.Complexity += myPC.L;
		myPC.Latency += myPC.log2L;
	}

	PolarCodeEncoding(x_hat, u_hat);
}

void Ens_SPRLD(PC &myPC, vector<double> y, vector<uint8_t> &u_hat)
{
	static double best_PM, tmp_PM;
	static vector<uint8_t> tmp_x_hat(myPC.N);
	static vector<uint8_t> tmp_u_hat(myPC.N);

	static bool isCompCal = false;
	static uint32_t complexity = 0;
	static uint32_t latency = 0;
	static double mem = 0.0;

	best_PM = MAXFLOAT;
	for (auto p = 0; p < myPC.P_max; p++)
	{
		myPC.Complexity = 0;
		myPC.Latency = 0;
		myPC.Memory = 0.0;
		SPRLD(myPC, y, tmp_u_hat, tmp_PM);
		PolarCodeEncoding(tmp_u_hat, tmp_x_hat);

		if (best_PM > tmp_PM)
		{
			best_PM = tmp_PM;
			u_hat = tmp_u_hat;
		}

		if (p==0)
		{
			latency = myPC.Latency + (uint32_t)ceil(log2(myPC.P_max));
			complexity = myPC.P_max * myPC.Complexity + (uint32_t)(myPC.P_max > 1 ? myPC.P_max : 0);
			if (myPC.L == 1 && myPC.P_max > 1)
				complexity += myPC.P_max * myPC.N;
			if (!isCompCal)
			{
				isCompCal = true;
				if (myPC.policy == "par")
					if (myPC.L > 1)
						mem = myPC.N * (myPC.log2N * myPC.L) * 32.0 + myPC.log2N * 32.0 + 2 * myPC.N * myPC.L;
					else
						mem = myPC.log2N * myPC.N * 32.0 + myPC.log2N * 32.0 + myPC.N;
				else if (myPC.policy == "seq")
					if (myPC.L > 1)
						mem = myPC.N * myPC.L * 32.0 + myPC.log2N * 32.0 + 2 * myPC.N * myPC.L;
					else
						mem = myPC.N * 32.0 + myPC.log2N * 32.0 + myPC.N;
				mem = (mem * myPC.P_max + myPC.N * 32.0) / 8192.0; // memory requirement in KBytes
			}
		}
	}
	myPC.Complexity = complexity;
	myPC.Latency = latency;
	myPC.Memory = mem;
}

void SSC_FHT(PC &myPC, vector<double> &y, vector<uint8_t> &u_hat)
{
	vector<vector<double>> LLR(myPC.log2N + 1, vector<double>(myPC.N));
	vector<uint8_t> tmp_x(myPC.N >> 1);
	uint32_t Function;
	uint32_t StageIdx;
	uint32_t UpBrach;
	uint32_t DownBranch;
	uint32_t N_stage;

	LLR[myPC.log2N] = y;
	for (uint32_t step = 0; step < myPC.FastSC_FHT_FunctionVec.size(); step++)
	{
		Function = myPC.FastSC_FHT_FunctionVec[step].name;
		StageIdx = myPC.FastSC_FHT_FunctionVec[step].stage;
		UpBrach = myPC.FastSC_FHT_FunctionVec[step].upBranch;
		DownBranch = myPC.FastSC_FHT_FunctionVec[step].lowBranch;
		N_stage = DownBranch - UpBrach + 1;

		if (Function == SPC)
		{
			// Decode the parent node
			double minABS_LLR = 9999.99;
			uint32_t minABS_LLR_idx = 0;
			bool SPC_sum = 0;

			for (uint32_t i = UpBrach; i <= DownBranch; i++)
			{
				tmp_x[i - UpBrach] = HARD_DECISION(LLR[StageIdx][i]);
				SPC_sum ^= tmp_x[i - UpBrach];
				if (abs(LLR[StageIdx][i]) < minABS_LLR)
				{
					minABS_LLR = abs(LLR[StageIdx][i]);
					minABS_LLR_idx = i;
				}
			}
			myPC.Complexity += N_stage;
			myPC.Latency += StageIdx;

			if (SPC_sum)
				tmp_x[minABS_LLR_idx - UpBrach] ^= 1;

			// Decode the information bits
			for (uint32_t i = 0; i < StageIdx; ++i)									  // stage index
				for (uint32_t j = 0; j < (N_stage >> (i + 1)); ++j)					  // XOR-group index
					for (uint32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l) // XOR index in each XOR Group
						tmp_x[l] ^= tmp_x[l + (1 << i)];

			for (uint32_t i = UpBrach; i <= DownBranch; i++)
				u_hat[i] = tmp_x[i - UpBrach];
		}
		else if (Function == HDM)
		{
			vector<double> tmp_y(N_stage);
			vector<uint8_t> tmp_u_hat(N_stage);
			vector<uint8_t> tmp_x_hat(N_stage);

			for (auto i = UpBrach; i <= DownBranch; i++)
				tmp_y[i - UpBrach] = LLR[StageIdx][i];

			// ML decoding for the constituent order-1 RM codes
			FHT(myPC.RelList, tmp_y, StageIdx, true, tmp_u_hat, tmp_x_hat);
			// fast Hadamard transform complexity and latency
			myPC.Complexity += N_stage * StageIdx;
			myPC.Latency += StageIdx;
			// sorting complexity and latency to get the best codeword
			myPC.Complexity += N_stage;
			myPC.Latency += StageIdx;

			for (auto i = UpBrach; i <= DownBranch; i++)
				u_hat[i] = tmp_u_hat[i - UpBrach];
		}
		else if (Function == FF)
		{
			for (uint32_t i = UpBrach; i <= DownBranch; i++)
				PolarCodeFFunction(LLR[StageIdx + 1][i],
								   LLR[StageIdx + 1][i + (1 << StageIdx)],
								   LLR[StageIdx][i]);
			myPC.Complexity += N_stage;
			myPC.Latency += 1;
		}
		else if (Function == GF)
		{
			for (uint32_t i = 0; i < N_stage; i++)
				tmp_x[i] = u_hat[UpBrach - N_stage + i];

			for (uint32_t i = 0; i < StageIdx; ++i)									  // stage index
				for (uint32_t j = 0; j < (N_stage >> (i + 1)); ++j)					  // XOR-group index
					for (uint32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l) // XOR index in each XOR Group
						tmp_x[l] ^= tmp_x[l + (1 << i)];

			for (uint32_t i = UpBrach; i <= DownBranch; i++)
				PolarCodeGFunction(LLR[StageIdx + 1][i - (1 << StageIdx)],
								   LLR[StageIdx + 1][i],
								   tmp_x[i - UpBrach],
								   LLR[StageIdx][i]);

			myPC.Complexity += N_stage;
			myPC.Latency += 1;
		}
		else if (Function == LF)
			u_hat[UpBrach] = HARD_DECISION(LLR[StageIdx][UpBrach]) * (1 - myPC.Frz[UpBrach]);
	}
}

void Aut_SSC_FHT(PC &myPC, vector<double> y, vector<uint8_t> &u_hat)
{
	static double best_PLM, tmp_PLM;
	static vector<uint32_t> PerStage(myPC.log2N);
	static vector<uint32_t> PerIdx(myPC.N);
	static vector<vector<uint32_t>> PerIdxList;
	static vector<double> y_per(myPC.N);
	static vector<uint8_t> u_hat_per(myPC.N, 0);
	static vector<uint8_t> u_hat_tmp(myPC.N, 0);
	static vector<uint8_t> x_hat_per(myPC.N);
	static vector<uint8_t> tmp_x_hat(myPC.N);
	static vector<uint8_t> tmp_u_hat(myPC.N);

	static bool init = false;

	static bool isCompCal = false;
	static uint32_t complexity = 0;
	static uint32_t latency = 0;
	static double mem = 0.0;

	if (!init)
	{
		init = true;

		string line;
		string config_file = "../../config/rmcodes/m" + to_string(myPC.log2N) + ".per";
		ifstream file(config_file);

		if (file.is_open())
			while (getline(file, line))
			{
				uint32_t char_cnt = 0;
				uint32_t rel = 0;
				char delimiter = ' ';
				for (uint32_t i = 0; i < myPC.N; i++)
				{
					sscanf(&line[char_cnt], "%d", &rel);
					PerIdx[i] = rel;
					// find the next reliability element position
					while (line[char_cnt] != delimiter)
						char_cnt++;
					char_cnt++;
				}
				PerIdxList.push_back(PerIdx);
			}

		// Here S is used as the number of parallel SSC-FHT decoders.
		myPC.L = clip(myPC.L, (uint32_t)1, myPC.P_max);
	}

	best_PLM = -MAXFLOAT;
	for (auto p = 0; p < myPC.P_max; p++)
	{
		PerIdx = PerIdxList[rand() % PerIdxList.size()];
		for (auto i = 0; i < myPC.N; i++)
			y_per[PerIdx[i]] = y[i];

		myPC.Latency = 0;
		myPC.Complexity = 0;
		SSC_FHT(myPC, y_per, u_hat_per);
		PolarCodeEncoding(u_hat_per, x_hat_per);

		tmp_PLM = 0;
		for (auto i = 0; i < myPC.N; i++)
			tmp_PLM += (1 - 2 * x_hat_per[i]) * y_per[i];
		myPC.Complexity += myPC.N;
		myPC.Latency += 1;

		if (tmp_PLM > best_PLM)
		{
			for (auto i = 0; i < myPC.N; i++)
				tmp_x_hat[i] = x_hat_per[PerIdx[i]];
			PolarCodeEncoding(tmp_x_hat, tmp_u_hat);
			best_PLM = tmp_PLM;
			u_hat = tmp_u_hat;
		}

		if (!isCompCal)
		{
			complexity = myPC.P_max * myPC.Complexity + myPC.P_max;
			latency = (uint32_t)ceil(log2(myPC.P_max)) + (uint32_t)(myPC.Latency * ceil(myPC.P_max / myPC.L));
			mem = ((myPC.N * 32.0 + myPC.N) * myPC.L + (myPC.P_max + myPC.N) * 32.0) / 8192.0; // memory requirement in KBytes
			isCompCal = true;
		}
	}
	myPC.Complexity = complexity;
	myPC.Latency = latency;
	myPC.Memory = mem;
}