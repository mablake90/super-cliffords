import stim
import numpy as np


"""
This file includes all functions used to compute the entropy. The command at the bottom of the file allows one to run the circuits and compute the entropy.

The function gf2_rank(rows) is taken from the page: https://stackoverflow.com/questions/56856378/fast-computation-of-matrix-rank-over-gf2
"""

def sample_stabilisers(s):
    """
    - Purpose: Use stim to simulate a circuit.
    -Inputs:
         - s (stim.circuit): Any circuit you like which you wish to simulate.
    -Outputs:
         - zs2 (np.ndarray): The result of conjugating the Z generators by the given stim circuit.     
    """
    tableau: stim.Tableau = s.current_inverse_tableau() ** -1
    n = len(tableau)
    zs = [tableau.z_output(k) for k in range(n)]
    zs2 = np.array(zs)
    return zs2
    
    
def binary_matrix(zStabilizers):
    """
        - Purpose: Construct the binary matrix representing the stabilizer states.
        - Inputs:
            - zStabilizers (np.ndarray): The result of conjugating the Z generators on the initial state.
        Outputs:
            - binaryMatrix (np.ndarray of size (N, 2N)): An array that describes the location of the stabilizers in the tableau representation.
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
    """
    Purpose: Convert signs into binary format.
    Inputs:
        - signs (np.ndarray): The signs of the stabilizer states.
    Outputs:
        - signs (np.ndarray): The signs of the stabilizer states in binary format.
    """
    n = np.size(signs)
    for r in range(n):    
        if (signs[r] == 1):
            signs[r] = 0
        else:
            signs[r] = 1
    return signs

def get_cut_stabilizers(binaryMatrix, cut):
    """
        - Purpose: Return only the part of the binary matrix that corresponds to the qubits we want to consider for a bipartition.
        - Inputs:
            - binaryMatrix (np.ndarray of size (N, 2N)): The binary matrix for the stabilizer generators.
            - cut (integer): Location for the cut.
        - Outputs:
            - cutMatrix (np.ndarray of size (N, 2cut)): The binary matrix for the cut on the left.
    """
    sh = binaryMatrix.shape
    assert sh[1] == 2*sh[0], "shape of array must be (N, 2N)"
    assert cut < sh[0], "cut must be less than N"
    N = len(binaryMatrix)
    cutMatrix = np.zeros((N, 2*cut))
    cutMatrix[:,:cut] = binaryMatrix[:,:cut]
    cutMatrix[:,cut:] = binaryMatrix[:,N:N+cut]
    return cutMatrix 


def rows(binMatrix):
    """
    - Purpose: Take a binary matrix and convert it into a list of integers, corresponding to the rows of the binary matrix expressed as integers
    conversion method: bigendian.
    - Inputs:
        - binMatrix (np.ndarray): a binary array of any size.
    - Outputs:
        - v (np.ndarray): length of array is number of rows of the binary matrix. All entries will be integers.  
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


def compute_entropy(s: stim.Circuit, cut: int):
    """
    - Purpose: Compute the entropy of a circuit.
    - Inputs:
        - s (stim.Circuit): The circuit you wish to compute the entropy of.
        - cut (integer): The cut across which to compute the entropy.
    - Outputs:
        - S (float): The entropy of the circuit.
    """
    zs2 = sample_stabilisers(s)
    mat = binary_matrix(zs2)
    b2 = get_cut_stabilizers(mat, cut)
    b3 = rows(b2)
    S = gf2_rank(b3.copy()) - cut
    return S