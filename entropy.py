from typing import List, Collection
import stim
import numpy as np
import gates
import matplotlib.pyplot as plt
import time
import math
import circuits


"""
This file includes all functions used to compute the entropy. The command at the bottom of the file allows one to run the circuits and compute the entropy.

The function gf2_rank(rows) is taken from the page: https://stackoverflow.com/questions/56856378/fast-computation-of-matrix-rank-over-gf2
"""

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
    return zs2
    
    
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

def convert_signs(signs):
    
    
    n = np.size(signs)
    for r in range(n):    
        if (signs[r] == 1):
            signs[r] = 0
        else:
            signs[r] = 1
        
    return signs



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
    - Purpose: Take a binary matrix and convert it into a list of integers, corresponding to the rows of the binary matrix expressed as integers.
    - Inputs:
        - binMatrix: a binary matrix of any size.
    - Outputs:
        - v (array): length of array is number of rows of the binary matrix. All entries will be integers.  
    """
    
    N = np.shape(binMatrix)[0]
    v = []
    for i in range(0, N):
         test_list = [int(a) for a in binMatrix[i] ]
         test_list.reverse()
         v.append(int("".join(str(x) for x in test_list), 2))
        
    return v


        
    


def gf2_rank(rows):
    """
    - Purpose: Finds rank of a Binary matrix over F2.
     - Inputs: 
        - rows of Binary matrix are bit strings, expressed as integers.
     - Outputs:
        - an integer, the rank of the matrix.
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
    - rep (integer): number of times to repeat the simulation and average over.
    - res (integer): resolution (i.e. how often to compute the operator entanglement).
    - slow (integer): determines what proportion of total qubits to act on at each timestep.
    
    Also make a choice of circuit from the circuits.py file.
    
    """


    startTime = time.time()
    """
    To compute the entropy of a circuit, first make a choice of random circuit. The two choices that we have below are LocInt, a circuit which acts on O(N) qubits with nearest neighbour interactions. Or FS3_Np a circuit that acts on O(N) qubits with all-to-all interactions. In order to compute the entropy one needs to make the following choices:
    - N - integer, length of the chain of qubits.
    - T - integer, number of timesteps for the evolution of the circuit.
    - rep - integer, number of repititions to smooth over fluctutations.
    - res - integer, resolution that specifies how often to compute the entropy.
    - cut - integer (less than N), specifies the cut across which to compute the entanglement entropy.
    - slow - integer, fixes the proportion of the chain on which to act on at each timestep.

    """        

    N = 120
    T = 100
    rep = 10
    res = 2
    cut = int(N/3)
    slow = 1
    v, w = circuits.runLocInt(N, T, rep, res, slow, cut) #Circuit with nearest neighbour interactions.
#    v, w = circuits.runFS3_Np(N, T, rep, res, slow, cut) #Circuit with all to all interactions.
    
    totTime = (time.time() - startTime)
    print('Execution time in seconds:' + str(totTime))
    

                 	      
    plt.plot(w, v)
    plt.xlabel('Time')
    plt.ylabel('S')
    plt.title(f'Operator Entanglement for N = {N}')
    
    plt.show()  
    
if (__name__ == '__main__'):
    main()     
 
    

    
