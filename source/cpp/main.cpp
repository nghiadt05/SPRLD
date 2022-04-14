#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
    /* Declare polar and simulation classes */
    PC myPC;
    SIM mySIM;
    myPC.code = "rm"; // polar codes: 'pc', reed-muller codes: 'rm'
    myPC.N = 256;     // Code length
    myPC.K = 37;      // Message length
    myPC.C = 0;       // CRC length
    myPC.L = 2;       // List size
    myPC.S = 200;
    myPC.CS = 0;
    myPC.CE = 0;
    myPC.I_max = 1;               // Number of maximum iterations for iterative decoders
    myPC.I_min = 1;               // Number of minimum iterations for iterative decoders
    myPC.T_max = 50;              // Number of maximum flipping attempts for bit-flipping decoders
    myPC.P_max = 10;              // Number of maximum permutations
    myPC.S_max = 0;               // Number of permutable last stages
    myPC.k = 100;                 // Number of bandit arms
    myPC.LLR_Max = 1e6;           // Maximum absolute LLR value
    myPC.Decoder = "Aut_SSC_FHT"; // Decoder
    myPC.policy = "per";          // RL policy for SSCF
    myPC.crc_sum = 0;             // CRC check sum
    myPC.isGenData = false;       // Gen dataset
    myPC.isZeroCW = false;        // Use all-zero codeword or not
    myPC.isCS = false;            // Critical set flag
    myPC.isEvalX = false;         // Validate x or u
    myPC.isGenie = false;
    myPC.PMod = "Eb"; // Power nomalization mode

    mySIM.RandSeed = 0;           // Random seed for noise generating
    mySIM.SNR_START = 0;          // Start SNR for simulation
    mySIM.SNR_END = 0;            // Stop SNR for simulation
    mySIM.SNR_STEP = 0.5;         // SNR step for simulation
    mySIM.Channel = "Linear";     // Channel modelling
    mySIM.Noise = "AWGN";         // Noise modelling
    mySIM.Mod = "BPSK";           // Modulation technique
    mySIM.MAX_FRAMES = 100000000; // Maximum simulated frames
    mySIM.MIN_FRAMES = 10000;     // Minimum simulated frames
    mySIM.MIN_ERRORS = 100;       // Minimum error frames
    mySIM.STABLE_ERRORS = 50;     // Stable errors if exceed MAX_FRAMES
    mySIM.TARGET_FER = 1e-4;      // Target FER
    mySIM.FC = 0;                 // Frame counts
    mySIM.FEC = 0;                // Frame error counts
    mySIM.BEC = 0;                // Bit error counts
    mySIM.FER = 0;                // Frame error rate
    mySIM.BER = 0;                // Bit error rate
    mySIM.TP = 0;                 // Decoding through puts
    mySIM.AvgT = 0;               // Average number of decoding attempts
    mySIM.isChannelLLR = 1;       // Convert the channel output to LLR or not
    mySIM.UFrame = 100;

    /* Parse new parameters */
    for (uint32_t i = 0; i < argc; i++)
    {
        // Print help messages
        if (!strcmp(argv[i], "--help"))
        {
            printf("--N: size of the code\n");
            printf("--K: number of information bits\n");
            printf("--L: list size for SCL decoding\n");
            printf("--SNR: start, stop, and step of the SNR samples\n");
            return EXIT_SUCCESS;
        }
        // Config polar code class (please ignore this part as many parameters are not used by SPRLD)
        else if (!strcmp(argv[i], "--code"))
            myPC.code = argv[i + 1];
        else if (!strcmp(argv[i], "--N"))
            myPC.N = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--K"))
            myPC.K = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--m"))
            myPC.RM_m = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--r"))
            myPC.RM_r = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--C"))
            myPC.C = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--L"))
            myPC.L = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--S"))
            myPC.S = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--CS"))
            myPC.CS = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--CE"))
            myPC.CE = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--I_max"))
            myPC.I_max = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--I_min"))
            myPC.I_min = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--T_max"))
            myPC.T_max = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--P_max"))
            myPC.P_max = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--S_max"))
            myPC.S_max = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--n_arm"))
            myPC.k = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--decoder"))
            myPC.Decoder = argv[i + 1];
        else if (!strcmp(argv[i], "--policy"))
            myPC.policy = argv[i + 1];
        else if (!strcmp(argv[i], "--Pmod"))
            myPC.PMod = argv[i + 1];
        else if (!strcmp(argv[i], "--quan"))
            myPC.quan = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--data"))
            myPC.isGenData = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--zero-cw"))
            myPC.isZeroCW = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--eval-x"))
            myPC.isEvalX = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--isCS"))
            myPC.isCS = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--isGenie"))
            myPC.isGenie = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--train"))
            myPC.isTraining = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--B"))
            myPC.B = (bool)atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--lr"))
            myPC.lr = (double)atof(argv[i + 1]);
        // Config simulation class
        else if (!strcmp(argv[i], "--Seed"))
            mySIM.RandSeed = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--SNR"))
            sscanf(argv[i + 1],
                   "%lf,%lf,%lf",
                   &mySIM.SNR_START,
                   &mySIM.SNR_END,
                   &mySIM.SNR_STEP);
        else if (!strcmp(argv[i], "--Channel"))
            mySIM.Channel = argv[i + 1];
        else if (!strcmp(argv[i], "--Noise"))
            mySIM.Noise = argv[i + 1];
        else if (!strcmp(argv[i], "--Mod"))
            mySIM.Mod = argv[i + 1];
        else if (!strcmp(argv[i], "--M_Frame"))
            mySIM.MAX_FRAMES = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--m_Frame"))
            mySIM.MIN_FRAMES = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--m_Error"))
            mySIM.MIN_ERRORS = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--s_Error"))
            mySIM.STABLE_ERRORS = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--UFrame"))
            mySIM.UFrame = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--T_FER"))
            mySIM.TARGET_FER = (double)atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--ChannelLLR"))
            mySIM.isChannelLLR = (bool)atoi(argv[i + 1]);
    }

    /* Config the polar class with new parameters */
    myPC.PC_Init();

    /* Configuration summary */
    printf("===============================================================================================================\n");
    printf("===========================================---Configuration summary---=========================================\n");
    printf("{code,N,K,C,L,S,P_max} = {%s,%d,%d,%d,%d,%d,%d}\n", myPC.code.c_str(), myPC.N, myPC.K, myPC.C, myPC.L, myPC.S, myPC.P_max);
    printf("Decoder: %s\n", myPC.Decoder.c_str());
    printf("Modulation: %s\n", mySIM.Mod.c_str());
    printf("Noise: %s\n", mySIM.Noise.c_str());
    printf("SNR = {%0.3f,%0.3f,%0.3f}\n", mySIM.SNR_START, mySIM.SNR_END, mySIM.SNR_STEP);
    printf("{Max #frames, Min #frames}={%d,%d}\n", mySIM.MAX_FRAMES, mySIM.MIN_FRAMES);
    printf("{Min #errors, Stable #errors}={%d,%d}\n", mySIM.MIN_ERRORS, mySIM.STABLE_ERRORS);
    printf("Target FER: %1.2E\n", mySIM.TARGET_FER);
    printf("==============================================================================================================\n");
    printf("==============================================================================================================\n");

    /* Channel modelling */
    vector<double> std_dev, SNR;
    double r = (double)myPC.K / (double)myPC.N;
    uint32_t Sample_Points = (uint32_t)((mySIM.SNR_END - mySIM.SNR_START) / mySIM.SNR_STEP + 1);
    if (mySIM.Noise == "AWGN")
    {
        for (uint32_t i = 0; i < Sample_Points; i++)
        {
            double tmpSNR = mySIM.SNR_START + mySIM.SNR_STEP * i;
            SNR.push_back(tmpSNR);
            if (myPC.PMod == "Eb")
                std_dev.push_back(1.0 / sqrt(pow(10, tmpSNR / 10.0) * 2 * r));
            else if (myPC.PMod == "Es")
                std_dev.push_back(1.0 / sqrt(pow(10, tmpSNR / 10.0)));
            else if (myPC.PMod == "Es_comp")
                std_dev.push_back(0.70710678118 / sqrt(pow(10, tmpSNR / 10.0)));
            else
            {
                printf("Power normalization mode not supported!");
                exit(EXIT_FAILURE);
            }
        }
    }

    /* Finding permutations for BP decoding */
    string tmpStr = myPC.PMod;
    tmpStr.append("/N0");
    printf("%5s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %11s \n",
           tmpStr.c_str(), "FC", "FEC", "BEC", "FER", "BER", "Complexity", "Time steps", "Mem[kBytes]");
    for (uint32_t i = 0; i < Sample_Points; i++)
    {
        // Reset error metrics
        mySIM.FC = 1;
        mySIM.FEC = 0;
        mySIM.BEC = 0;
        mySIM.TC = 0;
        mySIM.ComplexityC = 0;
        mySIM.LatencyC = 0;

        // Reset training flags if applicable
        myPC.curSRN = int(SNR[i] * 10);
        myPC.isFirstFrame = true;
        myPC.r_avg = 0.0;
        myPC.beta = myPC.beta_init;
        myPC.TrainCnt = 1;

        // Scalling factor for DSCF, SCLF decoders
        double var = pow(std_dev[i], 2);
        static double gamma = 0.4;
        if (!mySIM.isChannelLLR)
            gamma = 0.6 / var;

        // Random number generators initilization
        srandom(mySIM.RandSeed);
        srand(mySIM.RandSeed);
        default_random_engine m_gen(mySIM.RandSeed), noise_gen(mySIM.RandSeed);
        uniform_real_distribution<double> uniform(0.0, 1.0);
        normal_distribution<double> gaussian(0.0, 1.0);
        vector<double> y(myPC.N, 0.0);
        vector<uint8_t> u_hat(myPC.N, 0), x_hat(myPC.N, 0);

        while ((mySIM.FC < mySIM.MIN_FRAMES) || (mySIM.FEC < mySIM.MIN_ERRORS))
        {
            /* Performance variables */
            myPC.Attempt = 0;
            myPC.Latency = 0;
            myPC.Complexity = 0;

            /* Polar encoding */
            myPC.PC_Encode(m_gen, uniform);

            /* BPSK modulation, AWGN, and LLR calculation */
            for (uint32_t ii = 0; ii < myPC.N; ii++)
            {
                y[ii] = 1 - 2 * myPC.x[ii];                // BPSK
                y[ii] += std_dev[i] * gaussian(noise_gen); // AWGN
                if (mySIM.isChannelLLR)
                    y[ii] = (2 * y[ii]) / (pow(std_dev[i], 2)); // LLR
            }

            /* Explicitly run the decoder when only CE channel errors are presented */
            uint32_t error_cnt = 0;
            vector<uint8_t> tmp_x(myPC.N);
            for (auto i = 0; i < myPC.N; i++)
                if (HARD_DECISION(y[i]) != myPC.x[i])
                    error_cnt++;

            if (((error_cnt == myPC.CE) && (myPC.CE > 0)) || myPC.CE == 0)
            {
                /* Decoding algorithms */
                //----------------------------------
                //--------SCL-based decoders--------
                //----------------------------------
                if (myPC.Decoder == "Ens_SPRLD")
                    Ens_SPRLD(myPC, y, u_hat);
                else if (myPC.Decoder == "SPRLD")
                {
                    static double best_PM;
                    SPRLD(myPC, y, u_hat, best_PM);
                }
                else if (myPC.Decoder == "Aut_SSC_FHT")
                    Aut_SSC_FHT(myPC, y, u_hat);
                else if (myPC.Decoder == "RLDA")
                    RLDA(myPC, y, u_hat);
                else
                {
                    printf("Unsupported decoder!\n");
                    exit(EXIT_FAILURE);
                }

                /* Update throughput measurements */
                mySIM.TC += myPC.Attempt;
                mySIM.LatencyC += myPC.Latency;
                mySIM.ComplexityC += myPC.Complexity;

                /* Update BEC */
                uint32_t tmp_BEC = 0;
                if (myPC.isEvalX)
                    for (uint32_t ii = 0; ii < myPC.N; ii++)
                        tmp_BEC += DIFF_BIT(myPC.x[ii], x_hat[ii]);
                else
                    for (uint32_t ii = 0; ii < myPC.K; ii++)
                        tmp_BEC += DIFF_BIT(myPC.u[myPC.NonFrzIdx[ii]],
                                            u_hat[myPC.NonFrzIdx[ii]]);

                mySIM.BEC += tmp_BEC;

                /* Update FEC */
                if (tmp_BEC)
                    mySIM.FEC++;

                /* Update FER, BER and T_avg */
                mySIM.FER = (double)mySIM.FEC / (double)mySIM.FC;
                mySIM.BER = (double)mySIM.BEC / (long double)(mySIM.FC * myPC.K);
                mySIM.AvgComplexity = (long double)mySIM.ComplexityC / (double)mySIM.FC;
                mySIM.AvgLatency = (long double)mySIM.LatencyC / (double)mySIM.FC;

                /* Terminate simulation for the current SNR */
                if ((mySIM.FC >= mySIM.MAX_FRAMES) &&
                    (mySIM.FEC >= mySIM.STABLE_ERRORS))
                    break;

                /* Increase the frame count */
                mySIM.FC++;

                /* Print out the result */
                if ((mySIM.FC % mySIM.UFrame) == 0)
                {
                    printf("%5.3f | %10d | %10d | %10d | %10.4e | %10.4e | %10d | %10d | %11.2f \r",
                           SNR[i], mySIM.FC, mySIM.FEC, mySIM.BEC, mySIM.FER, mySIM.BER, (uint64_t)mySIM.AvgComplexity, (uint64_t)mySIM.AvgLatency, myPC.Memory);
                    fflush(stdout);
                }
            }
        }

        printf("%5.3f | %10d | %10d | %10d | %10.4e | %10.4e | %10d | %10d | %11.2f \n",
               SNR[i], mySIM.FC, mySIM.FEC, mySIM.BEC, mySIM.FER, mySIM.BER, (uint64_t)mySIM.AvgComplexity, (uint64_t)mySIM.AvgLatency, myPC.Memory);
        fflush(stdout);

        // Terminate simulation at high SNR
        if (mySIM.FER < mySIM.TARGET_FER)
        {
            printf("Simulation terminated!\n");
            exit(EXIT_SUCCESS);
        }
    }
    return EXIT_SUCCESS;
}
