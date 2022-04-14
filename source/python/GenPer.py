"""Generate the general codeword permutations for RM codes"""
import numpy as np

path = '../../config/rmcodes/m'

m_list = [2, 3, 4, 5, 6, 7, 8, 9]

for m in m_list:
    file = open(path + str(m) +".per","w")
    P_max = np.minimum(1024, np.math.factorial(m))
    for p in range(0, P_max):
        isInvertible = False
        while(not isInvertible):
            A = np.identity(m)
            for i in range(100):
                elementary_mat = np.identity(m)
                elementary_mat[np.random.randint(0, m, 1), np.random.randint(0, m, 1)]=1
                A = np.matmul(A, elementary_mat)
            detA = np.linalg.det(A)
            isInvertible = (detA != 0)
        b = np.random.randint(0, 2, (m,1))
        N_per = []
        for i in range(0, 1 << m):
            bin_i = np.binary_repr(i, width=m)
            tmp = []
            for c in bin_i:
                if(c=='0'):
                    tmp.append(0)
                else:
                    tmp.append(1)
            bin_i = np.reshape(np.array(tmp), (m, 1))
            bin_i_pi = (np.matmul(A, bin_i)+b)%2
            dec_i_pi = 0
            for j in range(m):
                dec_i_pi += bin_i_pi[j] * (1<<j)
            N_per.append(dec_i_pi)
            file.write("%d " % dec_i_pi)
        file.write("\n")
    file.flush()
    file.close()


