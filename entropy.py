from typing import List, Collection
import stim
import numpy as np
import gates
import matplotlib.pyplot as plt
import time
import math
import circuits

def sample_stabilizers(s):
    """
    - Purpose: Use stim to simulate a circuit.
    -Inputs:
         - s (stim.circuit): Any circuit you like which you wish to simulate.
    -Outputs:
         - zs2 (array): The result of conjugating the Z generators by the given stim circuit.     
    """

    tableau: stim.Tableau = s.current_inverse_tableau() ** -1
    n = len(tableau)
    zs = [tableau.z_output(k) for k in range(n)]
    zs2 = np.array(zs)
    signs = [(zs[k].sign).real for k in range(n)]
    signs2 = np.array(signs)
    return zs2, signs2
    
    
def binaryMatrix(zStabilizers):
    """
        - Purpose: Construct the binary matrix representing the stabilizer states.
        - Inputs:
            - zStabilizers (array): The result of conjugating the Z generators on the initial state.
        Outputs:
            - binaryMatrix (array of size (N, 2N)): An array that describes the location of the stabilizers in the tableau representation.
    """
    N = len(zStabilizers)
    binaryMatrix = np.zeros((N,2*N))
    r = 0 # Row number
    for row in zStabilizers:
        c = 0 # Column number
        for i in row:
            if i == 3: # Pauli Z
                binaryMatrix[r,N + c] = 1
            if i == 2: # Pauli Y
                binaryMatrix[r,N + c] = 1
                binaryMatrix[r,c] = 1
            if i == 1: # Pauli X
                binaryMatrix[r,c] = 1
            c += 1
        r += 1

    return binaryMatrix
    
def getCutStabilizers(binaryMatrix, cut):
    """
        - Purpose: Return only the part of the binary matrix that corresponds to the qubits we want to consider for a bipartition.
        - Inputs:
            - binaryMatrix (array of size (N, 2N)): The binary matrix for the stabilizer generators.
            - cut (integer): Location for the cut.
        - Outputs:
            - cutMatrix (array of size (N, 2cut)): The binary matrix for the cut on the left.
    """
    N = len(binaryMatrix)
    cutMatrix = np.zeros((N, 2*cut))

    cutMatrix[:,:cut] = binaryMatrix[:,:cut]
    cutMatrix[:,cut:] = binaryMatrix[:,N:N+cut]

    return cutMatrix 


def rows(binMatrix):
    """
    This probably does not work yet, but hopefully will give me a more direct way of computing entanglement entropy.

       
    """
    N = np.shape(binMatrix)[0]
    v = []
    for i in range(0, N):
         test_list = [int(a) for a in binMatrix[i] ]
         test_list.reverse()
         v.append(int("".join(str(x) for x in test_list), 2))
        
    return v

def rows2(zstab, N):

    v = []

    for rows in zstab:
        r = 0
        c = 0
        for i in rows:
            if i == 3:
                r+= 2**(2*N - c)
            if i == 2:
                r+= 2**(2*N - c) + 2**(N -c)
            if  i == 1:
                r+= 2**(N-c)
            c+= 1    
        v.append(r)

    return v    
        
    


def gf2_rank(rows):
    """
    - Purpose: Finds rank of a Binary matrix over F2.
     - Inputs: 
        - The rows of the matrix are given as nonnegative integers in an array, thought
    of as bit-strings.
     - Outputs:
        - an integer, the rank of the matrix.
        
    Note: This function modifies the input list. Use gf2_rank(rows.copy())
    instead of gf2_rank(rows) to avoid modifying rows.
    """
    rank = 0
    while rows:
        pivot_row = rows.pop()
        if pivot_row:
            rank += 1
            lsb = pivot_row & -pivot_row
            for index, row in enumerate(rows):
                if row & lsb:
                    rows[index] = row ^ pivot_row
    return rank





    
def main():
    """
    To compute operator entanglement of a circuit. Fix the following parameters:
    
    - N (integer): number of qubits.
    - T (integer): number of timesteps.
    - M (integer): number of times to repeat the simulation and average over.
    - l (integer): resolution (i.e. how often to compute the operator entanglement).
    
    Also make a choice of circuit from the circuits.py file.
    
    
    """


    startTime = time.time()

    N = 240 
    T = 200
    rep = 1
    res = 2
    cut = int(N/3)
    v, w = circuits.runFS3_Np(N, T, rep, res, 10)

       # np.savez(f'FS3_Entropy_M50_N{N}.npz', v)


    totTime = (time.time() - startTime)
    print('Execution time in seconds:' + str(totTime))
    
    
    A = np.min(np.argwhere(v > (cut - 2)))
    B = w[A]
    print(f'The Operator entanglement saturates after T={B} timesteps.')
                 	      
    plt.plot(w, v)
    plt.xlabel('Time')
    plt.ylabel('S')
    plt.title(f'Operator Entanglement for N = {N}')
    
    plt.show()  
    
if (__name__ == '__main__'):
    main()     
 
    

    
