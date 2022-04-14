#include "polar.h"
#include <assert.h>
using namespace std;

uint32_t n_choose_k(uint32_t n, uint32_t k)
{
    uint32_t tmp = 1;
    for (uint32_t i = 1; i <= k; i++)
        tmp = tmp * (n - k + i) / i;
    return tmp;
}

void dec2bin(uint32_t a, vector<uint8_t> &b)
{
    // Note that the binary expansion of $a$ is
    // stored in $b$ using little-edian convention.
    static uint32_t m, i;
    m = b.size();
    fill(b.begin(), b.end(), 0);
    i = 0;
    while (a > 0)
    {
        b[i] = a % 2;
        a >>= 1;
        i++;
        if (a > 0 && i == m)
        {
            cout << "\n b's length is too short for a" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void bin2dec(uint32_t &a, vector<uint8_t> b)
{
    a = 0;
    for (auto i = 0; i < b.size(); i++)
        a += b[i] << i;
}

double SPA(double a, double b)
{
    return 2 * atanh(tanh(a / 2) * tanh(b / 2));
}

double MIN_SUM(double a, double b)
{
    return 0.9375 * sgn(a) * sgn(b) * min(abs(a), abs(b));
}

double RPA_PROJECT(double a, double b)
{
    return log(exp(a + b) + 1) - log(exp(a) + exp(b));
}

void QUAN(double &x, double _2_to_m_, double min_quan, double max_quan)
{
    x = clip(round(x * _2_to_m_), min_quan, max_quan) / _2_to_m_;
}

void QUAN(double &x, uint64_t n, uint64_t m, bool s)
{
    uint32_t q = n + m;
    double _2_to_m_ = (double)(1 << m);
    double max, min = 0.0f;
    if (s)
    {
        double _2_to_q_1_ = (double)(1 << (q - 1));
        max = _2_to_q_1_ - 1;
        min = -max;
    }
    else
    {
        double _2_to_q_ = (double)(1 << q);
        max = _2_to_q_ - 1;
    }
    x = clip(round(x * _2_to_m_), min, max) / _2_to_m_;
}

void sort_indexes_dec(vector<double> &v, vector<uint32_t> &idx)
{
    // initialize original index locations
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
}

void sort_indexes_double(vector<double> &v, vector<uint32_t> &idx)
{
    // initialize original index locations
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
}

double BitErrProb(vector<double> &Alpha, vector<uint32_t> &NonFrzIdx, uint32_t InfoIdx)
{
    double error_prob = exp(-0.3 * abs(Alpha[NonFrzIdx[InfoIdx]]));
    for (uint32_t j = 0; j <= InfoIdx; j++)
        error_prob /= (1 + exp(-0.3 * abs(Alpha[NonFrzIdx[j]])));
    return error_prob;
}

void PolarCodeEncoding(vector<uint8_t> &u, vector<uint8_t> &x)
{
    uint32_t N = static_cast<uint32_t>(u.size());
    uint32_t log2N = 0;
    uint32_t N_tmp = N;
    while (N_tmp >>= 1)
        ++log2N;
    x = u;

    for (uint32_t i = 0; i < log2N; ++i)
        for (uint32_t j = 0; j < (N >> (i + 1)); ++j)
            for (uint32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l)
                x[l] ^= x[l + (1 << i)];
}

void PC_PSumUdate(vector<uint8_t> &u, vector<uint8_t> &x)
{
    unsigned int N = static_cast<unsigned int>(u.size());
    unsigned int log2N = 0;
    unsigned int N_tmp = N;
    while (N_tmp >>= 1)
        ++log2N;
    x = u;

    for (unsigned int i = 0; i < log2N; ++i)
        for (unsigned int j = 0; j < (N >> (i + 1)); ++j)
            for (unsigned int l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l)
                x[l] ^= x[l + (1 << i)];
}

void GetNodeType(vector<uint8_t> FrzVec, uint32_t &NodeType)
{
    bool isR0 = false, isR1 = false;
    bool isREP = false, isSPC = false;
    bool isHAM = false, isHMD = false;
    uint32_t K = 0;
    uint32_t N = FrzVec.size();
    uint32_t log2N = 0;

    while ((1 << log2N) < N)
        log2N++;

    for (uint32_t i = 0; i < N; i++)
        K += 1 - FrzVec[i];

    if (K == N)
        NodeType = R1;
    else if (K == 0)
        NodeType = R0;
    else if (K == N - 1)
        NodeType = SPC;
    else if (K == log2N+1)
        NodeType = HDM;
    else if (K == 1)
        NodeType = REP;
    else
        NodeType = RR;
}

int GetNodeType(uint32_t N, uint32_t log2N, vector<uint8_t> FrzVec,
                uint32_t begin_idx, uint32_t end_idx)
{
    bool isR0 = false, isR1 = false;
    bool isREP = false, isSPC = false;
    bool isHAM = false, isHMD = false;

    uint32_t K = 0;
    for (uint32_t i = begin_idx; i < end_idx; i++)
        K += 1 - FrzVec[i];

    if (K == N)
        return (int)R1;
    else if (K == 0)
        return (int)R0;
    else if (K == N - 1)
        return (int)SPC;
    else if (K == log2N+1)
        return (int)HDM;
    else if (K == 1)
        return (int)REP;
    else
        return (int)RR;
}

uint32_t vecdiff(vector<uint8_t> &a , vector<uint8_t> &b)
{
    static uint32_t diff = 0;
    diff = 0;
    for (auto i = 0; i < a.size(); i++)
        diff += (a[i] - b[i]) > 0 ? 1 : ((a[i] - b[i]) < 0 ? 1 : 0);
    return diff;
}

uint32_t GetNodeType(uint32_t r, uint32_t m)
{
    if (r == m)
        return R1;
    else if (m == -1)
        return R0;
    else if (r == m - 1)
        return SPC;
    // else if (r == 1)
    //     return HDM;
    else if (r == 0)
        return REP;
    else
        return FF;
}

PC::PC(void)
{
}

void PC::Get_RelVec(void)
{
    // Get the reliability vector
    string config_file;
    if (code == "pc")
        config_file = "../../config/pccodes/5G_N" + to_string(N) + ".pc";
    else if (code == "rm")
        config_file = "../../config/rmcodes/N" + to_string(N) + ".rm";
    string line;
    ifstream file(config_file);

    if (file.is_open())
    {
        getline(file, line);
        getline(file, line);
        getline(file, line);
    }
    else
    {
        cout << "Cannot open configuration file!" << endl;
        exit(EXIT_FAILURE);
    }

    uint32_t char_cnt = 0;
    uint32_t rel = 0;
    char delimiter = ' ';
    for (uint32_t i = 0; i < N; i++)
    {
        sscanf(&line[char_cnt], "%d", &rel);
        Rel.push_back(rel);
        // find the next reliability element position
        while (line[char_cnt] != delimiter)
            char_cnt++;
        char_cnt++;
    }

    file.close();
    file.clear();

    if (code == "rm")
    {
        vector<vector<uint32_t>> tmpRelList(log2N + 1);
        for (auto i = 1; i <= log2N; i++)
        {
            config_file = "../../config/rmcodes/N" + to_string(1 << i) + ".rm";
            file.open(config_file);
            if (file.is_open())
            {
                getline(file, line);
                getline(file, line);
                getline(file, line);
            }
            else
            {
                cout << "Cannot open configuration file!" << endl;
                exit(EXIT_FAILURE);
            }

            uint32_t char_cnt = 0;
            uint32_t rel = 0;
            char delimiter = ' ';
            for (uint32_t ii = 0; ii < (1 << i); ii++)
            {
                sscanf(&line[char_cnt], "%d", &rel);
                tmpRelList[i].push_back(rel);
                // find the next reliability element position
                while (line[char_cnt] != delimiter)
                    char_cnt++;
                char_cnt++;
            }

            file.close();
            file.clear();
        }
        RelList = tmpRelList;
    }

    // Get the frozen vector
    Frz.resize(N);
    for (uint32_t i = 0; i < N; i++)
    {
        if (i < (K + C))
            Frz[Rel[i]] = 0;
        else
            Frz[Rel[i]] = 1;
    }

    // Get the non-frozen bit-index vector
    for (uint32_t i = 0; i < N; i++)
    {
        if (Frz[i] == 0)
            NonFrzIdx.push_back(i);
        else
            FrzIdx.push_back(i);
    }
    std::sort(NonFrzIdx.begin(), NonFrzIdx.end());
    std::sort(FrzIdx.begin(), FrzIdx.end());

    // LUT from N to K domain
    std::fill(NtoK_LUT.begin(), NtoK_LUT.end(), 999);
    uint32_t k = 0;
    for (uint32_t n = 0; n < N; n++)
        if (!Frz[n])
        {
            NtoK_LUT[n] = k;
            k++;
        }
}

void PC::Get_CRCPoly(void)
{
    if (C)
    {
        CRCPoly.resize(C);
        switch (C)
        {
        case 1:
            CRCPoly[0] = 1;
            break;
        case 3:
            CRCPoly[0] = 1;
            CRCPoly[1] = 1;
            break;
        case 6:
            // 0b 0010 0001
            CRCPoly[5] = 1;
            CRCPoly[0] = 1;
            break;
        case 9:
            // 0b 01 0110 1111
            CRCPoly[8] = 1;
            CRCPoly[6] = 1;
            CRCPoly[5] = 1;
            CRCPoly[3] = 1;
            CRCPoly[2] = 1;
            CRCPoly[0] = 1;
            break;
        case 11:
            // 0b 0110 0010 0001
            CRCPoly[10] = 1;
            CRCPoly[9] = 1;
            CRCPoly[5] = 1;
            CRCPoly[0] = 1;
            break;
        case 16:
            // 0b 0001 0000 0010 0001
            CRCPoly[12] = 1;
            CRCPoly[9] = 1;
            CRCPoly[5] = 1;
            CRCPoly[0] = 1;
            break;
        case 24:
            // 0b 1011 0010 1011 0001 0001 0111;
            CRCPoly[23] = 1;
            CRCPoly[21] = 1;
            CRCPoly[20] = 1;
            CRCPoly[17] = 1;
            CRCPoly[15] = 1;
            CRCPoly[13] = 1;
            CRCPoly[12] = 1;
            CRCPoly[8] = 1;
            CRCPoly[4] = 1;
            CRCPoly[2] = 1;
            CRCPoly[1] = 1;
            CRCPoly[0] = 1;
            break;
        default:
            printf("Does not support %d-bit CRCs!\n", C);
            std::exit(EXIT_FAILURE);
            break;
        }
    }
}

void PC::SCScheduling(void)
{
    // fng_location: f = 0, nn = 1, ps = 2, g = 3, s = 4
    for (unsigned int i = 0; i <= log2N; i++)
    {
        // determine the function types and stages at each time location
        if (i == 0)
        {
            schedule.fspg_location.push_back(FF);
            schedule.fspg_location.push_back(SF);
            schedule.fspg_location.push_back(PS);
            schedule.fspg_stage.push_back(i);
            schedule.fspg_stage.push_back(i);
            schedule.fspg_stage.push_back(i);
        }
        else
        {
            // functions at time location
            schedule.fspg_location.insert(schedule.fspg_location.begin(), FF);
            schedule.fspg_location.push_back(GF);
            schedule.fspg_location.insert(schedule.fspg_location.end(), schedule.fspg_location.begin() + 2, schedule.fspg_location.end() - 1);
            // stage of the coressponding function at time location
            schedule.fspg_stage.insert(schedule.fspg_stage.begin(), i);
            schedule.fspg_stage.push_back(i - 1);
            schedule.fspg_stage.insert(schedule.fspg_stage.end(), schedule.fspg_stage.begin() + 2, schedule.fspg_stage.end() - 2);
            schedule.fspg_stage.push_back(i);
        }
    }

    // determine the branch positions of the function at each time location
    vector<vector<unsigned int>> tmpBranch(schedule.fspg_location.size());
    vector<unsigned int> f_cnt(log2N + 1, 0);
    vector<unsigned int> s_cnt(log2N + 1, 0);
    vector<unsigned int> p_cnt(log2N + 1, 0);
    vector<unsigned int> g_cnt(log2N + 1, 0);
    unsigned int upBranch = 0;
    unsigned int downBranch = 0;
    unsigned int curStage = 0;
    for (unsigned int i = 0; i < schedule.fspg_location.size(); i++)
    {
        // get the stage index
        curStage = schedule.fspg_stage[i];
        // check the function type
        if (schedule.fspg_location[i] == FF)
        {
            upBranch = (1 << (curStage + 1)) * f_cnt[curStage];
            downBranch = (1 << (curStage + 1)) * f_cnt[curStage] + (1 << curStage) - 1;
            tmpBranch[i].push_back(upBranch);
            tmpBranch[i].push_back(downBranch);
            f_cnt[curStage]++;
        }
        else if (schedule.fspg_location[i] == SF)
        {
            upBranch = s_cnt[curStage];
            downBranch = s_cnt[curStage];
            tmpBranch[i].push_back(upBranch);
            tmpBranch[i].push_back(downBranch);
            s_cnt[curStage]++;
        }
        else if (schedule.fspg_location[i] == PS)
        {
            upBranch = (1 << (curStage + 1)) * p_cnt[curStage];
            downBranch = (1 << (curStage + 1)) * p_cnt[curStage] + (1 << curStage) - 1;
            tmpBranch[i].push_back(upBranch);
            tmpBranch[i].push_back(downBranch);
            p_cnt[curStage]++;
        }
        else // GF
        {
            upBranch = (1 << (curStage + 1)) * g_cnt[curStage] + (1 << curStage);
            downBranch = (1 << (curStage + 1)) * (g_cnt[curStage] + 1) - 1;
            tmpBranch[i].push_back(upBranch);
            tmpBranch[i].push_back(downBranch);
            g_cnt[curStage]++;
        }
    }
    schedule.fspg_branch = tmpBranch;
}

void PC::RM_Form_QSpace(void)
{
    // Initialize the tmp SubSpace(space index, coset index, bit index)
    vector<vector<vector<uint32_t>>> tmpSubSpace(N - 1, vector<vector<uint32_t>>(N >> 1, vector<uint32_t>(2, 0)));
    vector<uint8_t> bin_diff(log2N, 0), a(log2N, 0);
    for (auto s = 0; s < N - 1; s++)
    {
        vector<bool> isInCoset(N, false);
        for (auto c = 0; c < N >> 1; c++)
        {
            if (c == 0)
            {
                tmpSubSpace[s][c][0] = 0;
                tmpSubSpace[s][c][1] = s + 1;
                isInCoset[0] = true;
                isInCoset[s + 1] = true;
                dec2bin(s + 1, bin_diff);
            }
            else
            {
                // Obtain the first element of the c-th coset
                for (auto i = 0; i < N; i++)
                    if (!isInCoset[i])
                    {
                        isInCoset[i] = true;
                        tmpSubSpace[s][c][0] = i;
                        break;
                    }
                // Obtain the second element of the c-th coset
                dec2bin(tmpSubSpace[s][c][0], a);
                for (auto i = 0; i < a.size(); i++)
                    a[i] ^= bin_diff[i];
                bin2dec(tmpSubSpace[s][c][1], a);
                isInCoset[tmpSubSpace[s][c][1]] = true;
            }
        }
    }
    QSpace = tmpSubSpace;

    RM_m = log2N;
    RM_r = 0;
    uint32_t tmp_K = 0;
    for (RM_r = 0; RM_r <= RM_m; RM_r++)
    {
        tmp_K += n_choose_k(RM_m, RM_r);
        if (tmp_K == K)
            break;
    }
}

vector<FunctionType> PC::FastSC_FHT_Scheduling(uint32_t offsetBranch,
                                              vector<uint8_t> &CurFrz)
{
    //
    vector<FunctionType> curFunctionVec;
    uint32_t curN = CurFrz.size();
    uint32_t log2_curN = (uint32_t)(log2(curN));

    // Check if the current node is a special node
    bool isR0 = true;
    bool isR1 = true;
    bool isREP = true;
    bool isSPC = true;
    bool isHDM = true;

    for (uint32_t i = 0; i < CurFrz.size(); i++)
    {
        // Detecting rate-0 node
        isR0 &= CurFrz[i];

        // Detecting rate-1 node
        isR1 &= !CurFrz[i];

        // Detecting repetition node
        if (i < CurFrz.size() - 1)
            isREP &= CurFrz[i];
        else
            isREP &= !CurFrz[i];

        // Detecting SPC node
        if (i == 0)
            isSPC &= CurFrz[i];
        else
            isSPC &= !CurFrz[i];

        // Detecting HDM code
        if (i < ((CurFrz.size() >> 1) - 1) || (i == (CurFrz.size() >> 1)))
            isHDM &= CurFrz[i];
    }

    // Uncomment these lines to turn off fast-node
    //isR0 = false;
    //isREP = false;
    //isR1 = false;
    //isSPC = false;
    //isHDM = false;

    // Return the node type if the current node is
    // a special node
    if (isR0 || isR1 || isSPC || isREP || isHDM)
    {
        FunctionType curFunction;
        if (isR0)
            curFunction.name = R0;
        else if (isREP)
            curFunction.name = REP;
        else if (isR1)
            curFunction.name = R1;
        else if (isSPC)
            curFunction.name = SPC;
        else if (isHDM)
            curFunction.name = HDM;

        curFunction.stage = log2_curN;
        curFunction.upBranch = offsetBranch;
        curFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(curFunction);
        return curFunctionVec;
    }
    else // Non-special node, need to go down the tree to see things
    {
        // Concat the F function to the list
        FunctionType Function;
        Function.name = FF;
        Function.stage = log2_curN - 1;
        Function.upBranch = offsetBranch;
        Function.lowBranch = offsetBranch + (curN >> 1) - 1;
        curFunctionVec.push_back(Function);

        // Leaf-node
        if (Function.stage == 0)
        {
            Function.name = LF;
            curFunctionVec.push_back(Function);
        }
        else
        {
            // Concat the functions needed for the left-child
            vector<FunctionType> curFVec;
            vector<uint8_t> FFrz = std::vector<uint8_t>(CurFrz.begin(),
                                                        CurFrz.begin() + (curN >> 1));
            curFVec = FastSC_FHT_Scheduling(offsetBranch, FFrz);
            for (uint32_t i = 0; i < curFVec.size(); i++)
                curFunctionVec.push_back(curFVec[i]);
        }

        // Concat the G function to the list
        FunctionType GFunction;
        GFunction.name = GF;
        GFunction.stage = log2_curN - 1;
        GFunction.upBranch = offsetBranch + (curN >> 1);
        GFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(GFunction);

        // Leaf-node
        if (GFunction.stage == 0)
        {
            GFunction.name = LF;
            curFunctionVec.push_back(GFunction);
        }
        else
        {
            // Concat the functions needed for the right-child
            vector<FunctionType> curGVec;
            vector<uint8_t> GFrz = vector<uint8_t>(CurFrz.begin() + (curN >> 1),
                                                   CurFrz.begin() + curN);
            curGVec = FastSC_FHT_Scheduling(offsetBranch + (curN >> 1), GFrz);
            for (uint32_t i = 0; i < curGVec.size(); i++)
                curFunctionVec.push_back(curGVec[i]);
        }

        // Concate the re-permute function to the list
        FunctionType RPFunction;
        RPFunction.name = PS;
        RPFunction.stage = log2_curN;
        RPFunction.upBranch = offsetBranch;
        RPFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(RPFunction);

        return curFunctionVec;
    }
}

vector<FunctionType> PC::RLDA_Scheduling(uint32_t offsetBranch, vector<uint8_t> &CurFrz)
{
    //
    vector<FunctionType> curFunctionVec;
    uint32_t curN = CurFrz.size();
    uint32_t log2_curN = (uint32_t)(log2(curN));

    // Check if the current node is a special node
    bool isREP = true;
    bool isSPC = true;

    for (uint32_t i = 0; i < CurFrz.size(); i++)
    {
        // Detecting repetition node
        if (i < CurFrz.size() - 1)
            isREP &= CurFrz[i];
        else
            isREP &= !CurFrz[i];

        // Detecting SPC node
        if (i == 0)
            isSPC &= CurFrz[i];
        else
            isSPC &= !CurFrz[i];
    }

    // Return the node type if the current node is
    // a special node
    if (isSPC || isREP)
    {
        FunctionType curFunction;
        if (isREP)
            curFunction.name = REP;
        else if (isSPC)
            curFunction.name = SPC;

        curFunction.stage = log2_curN;
        curFunction.upBranch = offsetBranch;
        curFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(curFunction);
        return curFunctionVec;
    }
    else // Non-special node, need to go down the tree to see things
    {
        // Concat the F function to the list
        FunctionType Function;
        Function.name = FF;
        Function.stage = log2_curN - 1;
        Function.upBranch = offsetBranch;
        Function.lowBranch = offsetBranch + (curN >> 1) - 1;
        curFunctionVec.push_back(Function);

        // Concat the functions needed for the left-child
        vector<FunctionType> curFVec;
        vector<uint8_t> FFrz = std::vector<uint8_t>(CurFrz.begin(),
                                                    CurFrz.begin() + (curN >> 1));
        curFVec = RLDA_Scheduling(offsetBranch, FFrz);
        for (uint32_t i = 0; i < curFVec.size(); i++)
            curFunctionVec.push_back(curFVec[i]);

        // Concat the G function to the list
        FunctionType GFunction;
        GFunction.name = GF;
        GFunction.stage = log2_curN - 1;
        GFunction.upBranch = offsetBranch + (curN >> 1);
        GFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(GFunction);

        // Leaf-node
        if (GFunction.stage == 0)
        {
            GFunction.name = LF;
            curFunctionVec.push_back(GFunction);
        }
        else
        {
            // Concat the functions needed for the right-child
            vector<FunctionType> curGVec;
            vector<uint8_t> GFrz = vector<uint8_t>(CurFrz.begin() + (curN >> 1),
                                                   CurFrz.begin() + curN);
            curGVec = RLDA_Scheduling(offsetBranch + (curN >> 1), GFrz);
            for (uint32_t i = 0; i < curGVec.size(); i++)
                curFunctionVec.push_back(curGVec[i]);
        }

        // Concate the re-permute function to the list
        FunctionType RPFunction;
        RPFunction.name = PS;
        RPFunction.stage = log2_curN;
        RPFunction.upBranch = offsetBranch;
        RPFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(RPFunction);

        return curFunctionVec;
    }
}

vector<FunctionType> PC::FastSCScheduling(uint32_t offsetBranch,
                                          vector<uint8_t> &CurFrz)
{
    //
    vector<FunctionType> curFunctionVec;
    uint32_t curN = CurFrz.size();
    uint32_t log2_curN = (uint32_t)(log2(curN));

    // Check if the current node is a special node
    bool isR0 = true;
    bool isR1 = true;
    bool isREP = true;
    bool isSPC = true;

    for (uint32_t i = 0; i < CurFrz.size(); i++)
    {
        // Detecting rate-0 node
        isR0 &= CurFrz[i];

        // Detecting rate-1 node
        isR1 &= !CurFrz[i];

        // Detecting repetition node
        if (i < CurFrz.size() - 1)
            isREP &= CurFrz[i];
        else
            isREP &= !CurFrz[i];

        // Detecting SPC node
        if (i == 0)
            isSPC &= CurFrz[i];
        else
            isSPC &= !CurFrz[i];
    }

    // Comment these lines to use fast-node
    // isR0 = false;
    // isREP = false;
    // isR1 = false;
    // isSPC = false;

    // Return the node type if the current node is
    // a special node
    if (isR0 || isR1 || isSPC || isREP)
    {
        FunctionType curFunction;
        if (isR0)
            curFunction.name = R0;
        else if (isREP)
            curFunction.name = REP;
        else if (isR1)
            curFunction.name = R1;
        else if (isSPC)
            curFunction.name = SPC;

        curFunction.stage = log2_curN;
        curFunction.upBranch = offsetBranch;
        curFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(curFunction);
        return curFunctionVec;
    }
    else // Non-special node, need to go down the tree to see things
    {
        // Concat the F function to the list
        FunctionType Function;
        Function.name = FF;
        Function.stage = log2_curN - 1;
        Function.upBranch = offsetBranch;
        Function.lowBranch = offsetBranch + (curN >> 1) - 1;
        curFunctionVec.push_back(Function);

        // Leaf-node
        if (Function.stage == 0)
        {
            Function.name = LF;
            curFunctionVec.push_back(Function);
            T_SCL++;
            T_SSCL++;
            T_FSSCL++;
        }
        else
        {
            // Concat the functions needed for the left-child
            vector<FunctionType> curFVec;
            vector<uint8_t> FFrz = std::vector<uint8_t>(CurFrz.begin(),
                                                        CurFrz.begin() + (curN >> 1));
            curFVec = FastSCScheduling(offsetBranch, FFrz);
            for (uint32_t i = 0; i < curFVec.size(); i++)
                curFunctionVec.push_back(curFVec[i]);
        }

        // Concat the G function to the list
        FunctionType GFunction;
        GFunction.name = GF;
        GFunction.stage = log2_curN - 1;
        GFunction.upBranch = offsetBranch + (curN >> 1);
        GFunction.lowBranch = offsetBranch + curN - 1;
        curFunctionVec.push_back(GFunction);

        // Leaf-node
        if (GFunction.stage == 0)
        {
            GFunction.name = LF;
            curFunctionVec.push_back(GFunction);
        }
        else
        {
            // Concat the functions needed for the right-child
            vector<FunctionType> curGVec;
            vector<uint8_t> GFrz = vector<uint8_t>(CurFrz.begin() + (curN >> 1),
                                                   CurFrz.begin() + curN);
            curGVec = FastSCScheduling(offsetBranch + (curN >> 1), GFrz);
            for (uint32_t i = 0; i < curGVec.size(); i++)
                curFunctionVec.push_back(curGVec[i]);
        }

        return curFunctionVec;
    }
}

vector<FunctionType> PC::OriginSCScheduling(uint32_t offsetBranch,
                                            vector<uint8_t> &CurFrz)
{
    vector<FunctionType> curFunctionVec;
    uint32_t curN = CurFrz.size();
    uint32_t log2_curN = (uint32_t)(log2(curN));

    // Concat the F function to the list
    FunctionType FFunction;
    FFunction.name = FF;
    FFunction.stage = log2_curN - 1;
    FFunction.upBranch = offsetBranch;
    FFunction.lowBranch = offsetBranch + (curN >> 1) - 1;
    curFunctionVec.push_back(FFunction);

    // Leaf-node
    if (FFunction.stage == 0)
    {
        FFunction.name = LF;
        curFunctionVec.push_back(FFunction);
    }
    else
    {
        // Concat the functions needed for the left-child
        vector<FunctionType> curFVec;
        vector<uint8_t> FFrz = std::vector<uint8_t>(CurFrz.begin(),
                                                    CurFrz.begin() + (curN >> 1));
        curFVec = OriginSCScheduling(offsetBranch, FFrz);
        for (uint32_t i = 0; i < curFVec.size(); i++)
            curFunctionVec.push_back(curFVec[i]);
    }

    // Concat the G function to the list
    FunctionType GFunction;
    GFunction.name = GF;
    GFunction.stage = log2_curN - 1;
    GFunction.upBranch = offsetBranch + (curN >> 1);
    GFunction.lowBranch = offsetBranch + curN - 1;
    curFunctionVec.push_back(GFunction);

    // Leaf-node
    if (GFunction.stage == 0)
    {
        GFunction.name = LF;
        curFunctionVec.push_back(GFunction);
    }
    else
    {
        // Concat the functions needed for the right-child
        vector<FunctionType> curGVec;
        vector<uint8_t> GFrz = vector<uint8_t>(CurFrz.begin() + (curN >> 1),
                                               CurFrz.begin() + curN);
        curGVec = OriginSCScheduling(offsetBranch + (curN >> 1), GFrz);
        for (uint32_t i = 0; i < curGVec.size(); i++)
            curFunctionVec.push_back(curGVec[i]);
    }

    return curFunctionVec;
}

void PC::Get_CRC_H_Mat(void)
{
    // CRC_H = [# parity checks, MSB->LSB]
    CRC_H.resize(C, (vector<uint8_t>)K_C);

    for (uint32_t i = 0; i < K; i++)
    {
        std::fill(m.begin(), m.end(), 0);
        m[i] = 1;
        CRC_Gen();
        for (uint32_t j = 0; j < C; j++)
            CRC_H[j][i] = crc[C - 1 - j];
    }

    for (uint32_t i = K; i < K + C; i++)
        CRC_H[K + C - 1 - i][i] = 1;

    // Get the V-C maps
    VC_MAP.resize(C);
    for (uint32_t j = 0; j < C; j++)
        for (uint32_t i = 0; i < K_C; i++)
            if (CRC_H[j][i])
                VC_MAP[j].push_back(i);

    // Get the C-V maps
    CV_MAP.resize(K_C);
    for (uint32_t i = 0; i < K_C; i++)
        for (uint32_t j = 0; j < C; j++)
            if (CRC_H[j][i])
                CV_MAP[i].push_back(j);
}

void PC::PC_Init(void)
{
    // Update encoder's and decoder's parameters
    // if (code == "rm")
    // {
    //     N = 1 << RM_m;
    //     K = 0;
    //     for(auto i=0; i<=RM_r; i++)
    //         K += nChoosek(RM_m, i);
    // }

    log2N = (uint32_t)log2((double)N);
    log2L = (uint32_t)log2((double)L);
    m.resize(K, 0);
    crc.resize(C);
    K_C = K + C;
    u.resize(N, 0);
    x.resize(N, 0);
    K_C = K + C;
    T_FSSCL = 0;
    T_SCL = 0;
    T_SSCL = 0;
    NtoK_LUT.resize(N, 0);
    pSum.resize(log2N + 1, vector<uint8_t>(N, 0));
    I_min = min(I_min, I_max);
    S_max = max(min(log2N, S_max), (uint32_t)0);

    if (code == "rm")
    {
        RM_Form_QSpace(); // Form the QSpace for an RM code
        // vector<uint32_t> y_idx(N), y_b_idx(N);
        // for (auto i = 0; i < N; i++)
        //     y_idx[i] = i;
        // for (auto i = 1; i <= N - 1; i++)
        //     RPA_Form_QSpace(y_idx, RM_m, i, y_b_idx);
    }

    Get_RelVec();  // Form the reliability vector
    Get_CRCPoly(); // Get the CRC polynomial vector
    if (C)
        Get_CRC_H_Mat(); // Get the CRC parity check matrix

    SCScheduling();
    RLDA_FunctionVec = RLDA_Scheduling(0, Frz);
    FastSC_FHT_FunctionVec = FastSC_FHT_Scheduling(0, Frz);
    FastFunctionVec = FastSCScheduling(0, Frz);
    OriginFunctionVec = OriginSCScheduling(0, Frz);
}

void PC::CRC_Gen(void)
{
    static vector<uint8_t> crc_dividend(K + C, 0);

    // Initialize the crc_dividend
    for (uint32_t i = 0; i < K; i++)
        crc_dividend[K + C - 1 - i] = m[i];

    for (uint32_t i = 0; i < C; i++)
        crc_dividend[i] = 0;

    // Generate the crc bits
    static uint32_t indx;
    indx = C + K - 1;
    while (indx >= C)
    {
        for (uint32_t i = indx; i >= C; i--)
        {
            if (crc_dividend[i] == 1)
            {
                crc_dividend[i] = 0;
                for (uint32_t j = 0; j < C; j++)
                    crc_dividend[i - 1 - j] ^= CRCPoly[C - 1 - j];
                break;
            }
        }
        indx--;
    }

    for (uint32_t i = 0; i < C; i++)
        crc[i] = crc_dividend[C - 1 - i];
}

void PC::CRC_Check(vector<uint8_t> &u_hat)
{
    vector<uint8_t> crc_dividend(K + C, 0);

    // Initialize the crc_dividend
    for (uint32_t i = 0; i < K_C; i++)
        crc_dividend[K_C - 1 - i] = u_hat[NonFrzIdx[i]];

    // Generate the crc bits
    static uint32_t indx;
    indx = C + K - 1;
    while (indx >= C)
    {
        for (uint32_t i = indx; i >= C; i--)
            if (crc_dividend[i] == 1)
            {
                crc_dividend[i] = 0;
                for (uint32_t j = 0; j < C; j++)
                    crc_dividend[i - 1 - j] ^= CRCPoly[C - 1 - j];
                break;
            }
        indx--;
    }

    crc_sum = 0;
    for (uint32_t i = 0; i < C; i++)
        crc_sum += crc_dividend[i];
    //printf("CRC sum: %d\n",crc_sum);
}

void PC::PC_Encode(default_random_engine &engine,
                   uniform_real_distribution<double> &uni_dist)
{
    // Generate message bits
    isZeroCW = false;
    for (uint32_t i = 0; i < K; i++)
        if (!isZeroCW)
            m[i] = (uint8_t)(uni_dist(engine) > 0.5 ? 1 : 0);
        else
            m[i] = 0;

    // Calculate the CRC bits
    if (C)
        CRC_Gen();

    // Construct the message word
    uint16_t m_cnt = 0, crc_cnt = 0;
    fill(u.begin(), u.end(), 0);
    for (uint32_t i = 0; i < NonFrzIdx.size(); i++)
    {
        if (code == "pc")
        {
            if (i < K)
            {
                u[NonFrzIdx[i]] = m[m_cnt];
                m_cnt++;
            }
            else
            {
                u[NonFrzIdx[i]] = crc[crc_cnt];
                crc_cnt++;
            }
        }
        else if (code == "rm")
            u[Rel[i]] = m[i];
    }

    // Construct the codeword
    x = u;
    pSum[0] = u;
    for (uint32_t i = 0; i < log2N; ++i)                                      // stage index
        for (uint32_t j = 0; j < (N >> (i + 1)); ++j)                         // XOR-group index
            for (uint32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l) // XOR index in each XOR Group
            {
                x[l] ^= x[l + (1 << i)];
                pSum[i + 1][l] = pSum[i][l] ^ pSum[i][l + (1 << i)];
                pSum[i + 1][l + (1 << i)] = pSum[i][l + (1 << i)];
            }
}

void PC::PC_HD_Decode(vector<double> &y, vector<uint8_t> &u_hat)
{
    // hard decision
    for (uint32_t i = 0; i < N; i++)
        u_hat[i] = HARD_DECISION(y[i]);

    // check for the syndrome
    syndrome = PC_SyndromeCheck(u_hat);

    // hard decoding
    for (int32_t i = log2N - 1; i >= 0; i--)
        for (int32_t j = 0; j < (N >> (i + 1)); ++j)
            for (int32_t l = (j << (i + 1)); l < (((j << 1) + 1) << i); ++l)
                u_hat[l] ^= u_hat[l + (1 << i)];

    for (uint32_t i = 0; i < FrzIdx.size(); i++)
        u_hat[FrzIdx[i]] = 0;

    // CRC verification
    if (C)
        CRC_Check(u_hat);
}

bool PC::PC_SyndromeCheck(vector<uint8_t> x_hat)
{
    vector<uint8_t> u_hat(N);
    PolarCodeEncoding(x_hat, u_hat);
    bool syndrome = false;
    for (auto i = 0; i < FrzIdx.size(); i++)
        if (u_hat[FrzIdx[i]])
        {
            syndrome = true;
            break;
        }
    return syndrome;
}