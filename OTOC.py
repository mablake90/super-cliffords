import stim
import numpy as np
import gates
import entropy
import time
import matplotlib.pyplot as plt
import random


def REF_binary(A, signs, N):
    """Purpose: Take a binary matrix and reduce it to Row echelon form (REF) over gf2. While also updating the signs that go with it using the appropriate rowsum operation.
       Inputs: - A: binary matrix (size N x 2N).
               - signs: a binary vector (corresponding to the signs).
               - N: an integer, the size of the vector signs.
       Outputs: - A: binary matrix in REF.
                - signs: the correct signs to accompany this."""
    
    n_rows, n_cols = A.shape

    # Compute row echelon form (REF)
    current_row = 0
    for j in range(n_cols):  # For each column
        if current_row >= n_rows:
            break
            
        pivot_row = current_row
     
        # find the first row in this column with non-zero entry. 
        # this becomes the pivot row
        while pivot_row < n_rows and A[pivot_row, j] == 0:
            pivot_row += 1
        
        # if we reach the end, this column cannot be eliminated. 
        if pivot_row == n_rows:
            continue  
            
        # otherwise, swap current row with the pivot row
        A[[current_row, pivot_row]] = A[[pivot_row, current_row]]
        a = signs[current_row]
        signs[current_row] = signs[pivot_row]
        signs[pivot_row] = a
        

        pivot_row = current_row
        current_row += 1
        
        # Eliminate rows below
        for i in range(current_row, n_rows):
            # subtract current row from any other rows beneath with 
            # a non-zero entry in this current column 
            if A[i, j] == 1:
                A[i] = (A[i] +  A[pivot_row]) % 2 # subtracting is same as adding in GF(2)
                signs[i] = row_sum(A[i], A[pivot_row], signs[i], signs[pivot_row], N)
                
              
    return A, signs


def g(x1, z1, x2, z2):
    "Purpose: Computes the function g needed for the rowsum operation.
     Inputs: -x1, z1, x2, z2 all are a binary digit.
     Outputs: -g: also a binary digit."

    if (x1 ==0) and (z1 == 0):
        g = 0
    if (x1 == 0) and (z1 == 1):
        g = x2*(1 - 2*z2)
    if (x1 == 1) and (z1 == 0):
        g = z2*(2*x2 -1)
    if (x1 == 1) and (z1 == 1):
        g = z2 - x2
        
    return g    


def row_sum(h, i, rh, ri, N):
    """
    Purpose: Computes the full rowsum operation.
    Inputs: - h: a vector (the row that is being summed into).
            - i: a vector (the row that we are summing to h).
            - rh: binary digit (the sign of h).
            - ri: binary digit (the sign of i).
            - N: an integer (the size of the vector).
    Outputs: - rh: the updated sign of h, using the rowsum operation.
    """

    m = np.size(h)

    k = 0
    for j in range(N):
        k += g(i[j], i[j+N], h[j], h[j+N])

    f = 2*rh +2*ri + k
    
    if (f % 4) == 0:
        rh = 0
    else:
        rh = 1
        
    return rh    


    

def xs(binMat):
    """Purpose: Takes a stabilizer matrix and just returns the half corresponding to the xs.
       Inputs: -binMat: a binary N x 2N matrix.
       Outputs: -xs: a binary N x N matrix (just the xs).
    """

    N = len(binMat)
    xs = np.zeros((N, N))
    xs[:,:] = binMat[:, :N]

    return xs


def small_zs(binMat, size2, N):
    """
    Purpose: Takes a full binary matrix and returns only a small corner of the zs.
    Inputs: -binMat: binary matrix.
            -size2: a number (the size we want to cut it down to).
            -N: current size of matrix.
    Outputs: -small_zs: the small matrix we have cut out.


    """

    small_zs = np.zeros((size2, N))
    small_zs[:,:] = binMat[N-size2:, N:]
    return small_zs


def OTOC_FS3_Np(N, T, rep, res, slow, Op):
    """
        - Purpose: Run a fast scrambling circuit and use it to compute an OTOC. Each time step this circuit acts on N/slow qubits. On the qubits which it acts on, it applies Z.H on 1/4 of the qubits and C3 on the remaing ones. In order to compute the OTOC one needs to introduce another operator, this is given by Op.

        - Inputs: 
               - N: integer (number of qubits).
               - T: integer (number of timesteps).
               - rep: integer (number of repetitions).
               - res: integer (resolution - i.e. how often the OTOC gets computed).
               - slow: integer (determines how much we slow down the action of the circuit).
               - Op: stim.TableauSimulator (gives a new operator which should be a super-Clifford and which we will use for computing the OTOC
        - Outputs: - v: a vector (the result of the OTOC calculation).
                   - w: a vector (the result of the time evolution).

"""
    Size = int(T/res) #specify size of output vectors.
    v = np.zeros(Size)
    w = np.zeros(Size)
    for m in range(Size):
        w[m] = m*res

    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS3_NpStep(N, s, slow)
            if (i % res) == 0:
                tableau1: stim.Tableau = s.current_inverse_tableau()**-1
                n1 = len(tableau1)
                tableau2: stim.Tableau = Op.current_inverse_tableau()**-1
                n2 = len(tableau2)
                tableau3: stim.Tableau = s.current_inverse_tableau()
                tableau_tot: stim.Tableau = (tableau3*tableau2)*tableau1
                n = len(tableau_tot)
                zs = [tableau_tot.z_output(k) for k in range(n)]
                zs2 = np.array(zs)
                signs = [(zs[k].sign).real for k in range(n)]
                signs2 = np.array(signs)
                bMat = entropy.binaryMatrix(zs2)
            
                signs3 = entropy.convert_signs(signs2)
               
                REF, signs3 = REF_binary(bMat, signs3, N)
            

                
                xs1 = xs(REF)
                rows = entropy.rows(xs1)
                rank = entropy.gf2_rank(rows)
                size2 = N - rank
            

                small = small_zs(REF, size2, N)

                REF2 = small #RREF_binary(small)
                shape = np.shape(REF2)

                signs4 = signs3[rank:]
                Ans = 0

                for k in range(size2):
                    if (signs4[k] == 1):
                          Ans = 1

                        
                if (Ans == 1):
                    v[int(i/res)] += 0
                else:    
                    v[int(i/res)] +=(2**(-(rank)/2))/rep        
                    
    return v, w           

def Op(N):
    """
    Purpose: Generate an operator which will act as V(0) in the OTOC calculation.
    Inputs: - N: an integer, the size of the chain of spins.
    Outputs: - s: a stim circuit with V(0) on it.
    """
    s = stim.TableauSimulator()

    c1 = stim.Circuit()

    c = stim.Circuit()
    r = random.randint(0, N-3)
    c.append_operation("I", [N-1])
    c.append_operation("I", [0])
    s.do(c)
    s.do(gates.ZH(r))
    s.do(gates.C3(r, r+1, r+2))

    return s

    

def main():
    """
    Below to run the OTOC calculation, first set V(0) in the Op(N) function, above. Then fix the following parameter: N (size of chain), T (number of timesteps), rep (number of repititions), res (resolution), slow (how much to slow down the fast scrambling circuit).
    """


    
    startTime = time.time()

    for i in range(0, 10):
        for k in range(0, 4):
            
            N = 120 + k*120
            print(N)
            T = 80
            rep = 200
            res = 2
            slow = 10
            Op1 = Op(N)
            v, w = OTOC_FS3_Np(N, T, rep, res, slow, Op1)
            np.savez(f'random_C3TT_OTOC_FS3_N{N}_slow{slow}_rep{rep}_iteration{i}.npz', v)


    N = 120
    T = 80
    rep = 500
    res = 2
    slow = 2
    Op1 = Op(N)
    v, w = OTOC_FS3_Np(N, T, rep, res, slow, Op1)

    totTime = (time.time() - startTime)
    print('Execution time in seconds:' + str(totTime))


                 	      
    plt.plot(w, v)
    plt.xlabel('Time')
    plt.ylabel('OTOC')
    plt.title(f'OTOC for N = {N}')
    
    plt.show()  
    


    



if (__name__ == '__main__'):
    main()     
 
        
