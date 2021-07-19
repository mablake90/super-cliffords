import stim
import numpy as np
import gates
import entropy
import time
import matplotlib.pyplot as plt


def RREF_binary(A):
    """Purpose: Converts a matrix to reduced row echelon form (RREF) by Gaussian elimination over gf2.
    
    - Inputs: 
         - A (np.array): a binary matrix to be reduced to RREF
    -Outputs:
         - A (np.array): a binary matrix in RREF.    
    """
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

        pivot_row = current_row
        current_row += 1
        
        # Eliminate rows below
        for i in range(current_row, n_rows):
            # subtract current row from any other rows beneath with 
            # a non-zero entry in this current column 
            if A[i, j] == 1:
                A[i] = (A[i] +  A[pivot_row]) % 2 # subtracting is same as adding in GF(2)
                
              
    # Compute reduced row echelon form (RREF)
    # in the RREF form, there is only one non-zero entry in a column. 
    for i in reversed(range(current_row)):
        # Find pivot
        pivot_col = 0
        
        # find the column with the first non-zero entry. 
        while pivot_col < n_cols and A[i, pivot_col]==0:
            pivot_col += 1
        if pivot_col == n_cols:
            continue  # Skip this all-zero row
    
        # Eliminate this column in all the rows above
        for j in range(i):
            if A[j, pivot_col] == 1:
                A[j] = (A[j] +  A[i]) % 2

    return A

def xs(binMat):
    """Purpose: Takes a matrix of size N x 2N and returns a matrix of size N x N (the first half of the matrix corresponding to the xs.
    - Inputs: 
         - binMat (np.array): a N x 2N binary matrix to be reduced to RREF.
    -Outputs:
         - xs (np.array): a binary matrix in RREF.    
    """

    N = len(binMat)
    xs = np.zeros((N, N))
    xs[:,:] = binMat[:, :N]

    return xs




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
        - Outputs: ???

"""
    Size = int(T/res) #specify size of output vectors.
    v = np.zeros(Size)
    w = np.zeros(Size)
    for m in range(Size):
        w[m] = m


    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.Id_Step(N,s)
            if (i % res) == 0:
                tableau1: stim.Tableau = s.current_inverse_tableau()**-1
                n1 = len(tableau1)
                tableau2: stim.Tableau = Op.current_inverse_tableau()**-1
                n2 = len(tableau2)
                tableau3: stim.Tableau = s.current_inverse_tableau()
                tableau_tot: stim.Tableau = tableau3*(tableau2*tableau1)

                n = len(tableau_tot)
                zs = [tableau_tot.z_output(k) for k in range(n)]
                zs2 = np.array(zs)
                signs = [(zs[k].sign).real for k in range(n)]
                signs2 = np.array(signs)
                bMat = entropy.binaryMatrix(zs2)
                RREF = RREF_binary(bMat)
                xs1 = xs(RREF)
                rows = entropy.rows(xs1)
                rank = entropy.gf2_rank(rows)
                                
                v[int(i/res)] +=(2**(-(rank)/2))/rep        
                    
    return v, w           

def Op(N):
    """ 
    Purpose: Create an operator to act as V(0) in the OTOC calculation.
    - Inputs: 
            - N (integer): length of the chain, is needed to make sure that tableau multiplication works.
    - Outputs: 
            - s (stim.TableauSimulator): a stim Tableau representation of the desired V(0).
  
    """
    s = stim.TableauSimulator()
         

    c = stim.Circuit()
    c.append_operation("I", [N-1])
    c.append_operation("H", [1])
    s.do(c)

    return s

    


def main():
    startTime = time.time()
                   
    

    N = 120
    T = 300
    rep = 1
    res = 1
    slow = 30
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
 
        
