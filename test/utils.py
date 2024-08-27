import numpy as np
import math


"""
Purpose of code is to run a test of the OTOC code. We will test it by constructing complicated random circuits on three qubits, we can then compute the OTOC directly as well as using the stabilizer formalism (as we use in OTOC.py). The figures generated in these two ways should match exactly.
"""

"""
First, create the operators we need as matrices on the underlying Hilbert space.
"""
X = np.array([[0,1],[1,0]])

I = np.array([[1,0],[0,1]])

Y = np.array([[0, complex(0,-1)],[complex(0, 1), 0]])

Z = np.array([[1,0],[0,-1]])

XXX = np.kron(X, np.kron(X, X))

YXX = np.kron(Y, np.kron(X, X))

YYX = np.kron(Y, np.kron(Y, X))

YXY = np.kron(Y, np.kron(X, Y))

XYX = np.kron(X, np.kron(Y, X))

XXY = np.kron(X, np.kron(X, Y))

XYY = np.kron(X, np.kron(Y, Y))

YYY = np.kron(Y, np.kron(Y, Y))

CX12 = np.kron(np.array([[1, 0,0,0],[0, 1,0,0],[0,0,0,1],[0,0,1,0]]), I)

CZ12 =np.kron( np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]), I)

CX21 = np.kron(np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]), I)

CX31 = np.array([[1,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,1,0,0,0],[0,1,0,0,0,0,0,0],[0,0,0,0,0,0,1,0],[0,0,0,1,0,0,0,0]])


T = np.array([[1,0], [0, (complex(1,1))*(1/math.sqrt(2))]])

T_trip = np.kron(T, np.kron(T,T))

T1 = np.kron(T, np.kron(I,I))

T1_6 = np.linalg.matrix_power(T1, 6)

T2 = np.kron(I, np.kron(T,I))

T3 = np.kron(I, np.kron(I, T))

T2_6 = np.linalg.matrix_power(T2, 6)

SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0,1,0,0],[0,0,0,1]])

SWAP12 = np.kron(SWAP, I)

SWAP23 = np.kron(I, SWAP)

SWAP13 = np.matmul(SWAP12, np.matmul(SWAP23, SWAP12))

C3_123 = np.matmul(CX21, np.matmul(CX31, np.matmul(CZ12, np.matmul(T1_6, T2_6))))

C3_123_dag = np.conjugate(np.transpose(C3_123))

C3_213 = np.matmul(SWAP12, np.matmul(C3_123, SWAP12))

C3_132 = np.matmul(SWAP23, np.matmul(C3_123, SWAP23))

C3_312 = np.matmul(SWAP13, np.matmul(C3_132, SWAP13))

C3_321 = np.matmul(SWAP13, np.matmul(C3_123, SWAP13))

C3_231 = np.matmul(SWAP12, np.matmul(C3_132, SWAP12))

H = np.array([[1, 1], [1, -1]])/math.sqrt(2)
S = np.array([[1, 0], [0, 1j]])
HXY = S @ X