# Successive-Permutation Recursive-List Decoding of Reed-Muller Codes

## Introduction
This package contains the C++ implementation of the successive-permutation (SP) recursive-list-decoding algorithm (SPRLD) of Reed-Muller (RM) codes, proposed in the ["Decoding Reed-Muller Codes with Successive Codeword Permutations"](https://arxiv.org/abs/2109.02122) paper.

Please consider citing this paper if the source code or the ideal of the paper is useful for your research. THANK YOU! :smile:

@ARTICLE{9906098,  author={Doan, Nghia and Hashemi, Seyyed Ali and Mondelli, Marco and Gross, Warren J.},  journal={IEEE Transactions on Communications},   title={Decoding Reed-Muller Codes with Successive Codeword Permutations},   year={2022},  volume={},  number={},  pages={1-1},  doi={10.1109/TCOMM.2022.3211101}}  

In addition, this package also provides the C++ implementations of state-of-the art RM decoders, namely:
- Automorphism simplified successive-cancellation decoding with fast Hadamard transform (Aut-SSC-FHT) [**Dumer'04**, **Geiselhart'21**] and,

- Recursive list decoding applied to the automorphism group of RM codes (RLDA) [**Dumer'06**, **Geiselhart'21**].

[**Dumer'04**] I. Dumer, "Recursive decoding and its performance for low-rate Reed-Muller codes," in IEEE Transactions on Information Theory, vol. 50, no. 5, pp. 811-823, May 2004, doi: [10.1109/TIT.2004.826632](https://ieeexplore.ieee.org/document/1291729).

[**Dumer'06**] I. Dumer and K. Shabunov, "Soft-decision decoding of Reed-Muller codes: Recursive lists," in IEEE Transactions on Information Theory, vol. 52, no. 3, pp. 1260-1266, March 2006, doi: [10.1109/TIT.2005.864443](https://ieeexplore.ieee.org/document/1603792).

[**Geiselhart'21**] M. Geiselhart, A. Elkelesh, M. Ebada, S. Cammerer and S. ten Brink, "Automorphism Ensemble Decoding of Reed–Muller Codes," in IEEE Transactions on Communications, vol. 69, no. 10, pp. 6424-6438, Oct. 2021, doi: [10.1109/TCOMM.2021.3098798](https://ieeexplore.ieee.org/document/9492151).

## Compilation
The C++ source codes are stored in ./source/cpp, while the python code used to generate the general codeword permutation is provided in ./source/python/GenPer.py. Here are the steps to compile the source codes for a Linux machine:
1. Navigate to ./source/cpp 
2. Compile and generate the executable binary file by running: 
"g++ -O3 -std=c++11 -o main *.cpp -lm" 
3. The codeword permutations are provided in the configuration folder ./config/rmcodes; however, you can generate them again by running the GenPer.py python script.

## Simulation of the proposed SPRLD decoder
Command to run simulation: "./main --N 512 --K 256 --SNR 1.5,3,0.5 --T_FER 1E-4 --decoder Ens_SPRLD --S 2 --L 4 --P_max 8 --policy seq".

Where  
--N N: code length,  
--K K: code dimension (number of information bits),  
--SNR start_SNR,end_SNR,step_SNR (by default SNR refers to Eb/N0, no space between the parameters),  
--T_FER T_FER: target FER to terminate the simulation,  
--S S: the number of the left-child RM codes that the SP technique is applied to (0&le;S&le;S_max), where S_max is the number of all the left-child nodes. If S>S_max the value of S will be automatically bounded by S_Max,  
--L L: list size,  
--P_max T (T&ge;1): the number of constituent SSP-RLD decoders that run in parallel,  
--policy policy: policy indicates either the sequential ("seq") or the parallel ("par") implementation of the proposed SP scheme, the specification of the policy will affect the memory and the decoding latency of the proposed decoder.  

## Simulation of the Aut-SSC-FHT decoder
Command to run simulation: "./main --N 512 --K 256 --SNR 1.5,3,0.5 --T_FER 1E-4 --decoder Aut_SSC_FHT --P_max 128 --L 32".  

Where  
--P_max P (P&ge;1): the total number of SSC-FHT decoding attempts.  
--L L (1&le;L&le;P): the number of constituent SSC-FHT decoders that can be run in parallel, the setting of L allows one to evaluate the latency and memory tradeoffs under serial (L=1), semi-parallel (1<L<P), and fully parallel (L=P) implementations of the Aut-SSC-FHT decoder. Note that the setting of L does not affect the FER performance of the Aut-SSC-FHT decoder.  

## Simulation of the RLDA decoder
Command to run simulation: "./main --N 512 --K 256 --SNR 1.5,3,0.5 --T_FER 1E-4 --decoder RLDA --L 64".

Where 
--L M (1&le;M): the list size used by the RLDA decoder.

