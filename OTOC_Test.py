import numpy as np
import math
import random
import OTOC
import matplotlib.pyplot as plt
import gates
import stim
import entropy


"""
Purpose of code is to run a test of the OTOC code. We will test it by constructing complicated random circuits on three qubits, we can then compute the OTOC directly as well as using the stabilizer formalism (as we use in OTOC.py). The figures generated in these two ways should match exactly.
"""

"""
First, create the operators we need as matrices on the underlying Hilbert space.
"""
X = np.array([[0,1],[1,0]])

I = np.array([[1,0],[0,1]])

Y = np.array([[0, complex(0,-1)],[complex(0, 1), 0]])

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


def dagger(U):
    """
    Takes the conjugate transpose of an operator.
    """
    return  np.conjugate(np.transpose(U))

def sqrt(n):
    return math.sqrt(n)



def Test_circuit_1Step(s, U):
    """
    Purpose: Move forward by one step of a random circuit.
    
    Inputs: - s: a stim circuit.
            - U: a matrix on the underlying Hilbert space, should be the same operator as that found in s.

    Outputs: -s: stim circuit, updated by one step of random circuit.
             - U: a matrix on underlying Hilbert space, updated by one step of random circuit.
    """


    r1 = random.randint(0, 2)
    s.do(gates.ZH(r1))


    if r1 == 0:
        U = np.matmul(U, T1)
    if r1 == 1:
        U = np.matmul(U, T2)
    if r1 == 2:
        U = np.matmul(U, T3)
        
   


    r3 = random.randint(0,5)
    
    if r3 == 0:
       s.do(gates.C3(0, 1, 2))
       U = np.matmul(U, C3_123)
    elif r3 ==1:
       s.do(gates.C3(1, 0, 2))
       U = np.matmul(U, C3_213)
    elif r3 == 2:
       s.do(gates.C3(1, 2, 0))
       U = np.matmul(U, C3_231)
    elif r3 == 3:
       s.do(gates.C3(0, 2, 1))
       U = np.matmul(U, C3_132)
    elif r3 ==4:
       s.do(gates.C3(2, 1, 0))
       U = np.matmul(U, C3_321)
    else:
       s.do(gates.C3(2, 0, 1))
       U = np.matmul(U, C3_312)

    return s, U     
    
def Op():
    """
    Purpose: Create an operator which will be the perturbation V(0) in the OTOC calculation.
    
    Inputs: - None.
    Outputs: - s: a stim circuit with operator V(0) encoded in it.
             - V0: a matrix on underlying Hilbert space representing the same operator.

    """
    
    s = stim.TableauSimulator()

    c1 = stim.Circuit()
 #   for i in range(int(N/3)):
 #       j = 3*i
 #       s.do(gates.C3(j, j+1, j+2))

    c = stim.Circuit()
    c.append_operation("I", [2])
    c.append_operation("I", [0])
    s.do(c)
    s.do(gates.C3(0,1,2))

    V0 = C3_123

    return s, V0

def M(U, W0, V0):
     return np.matmul(W0,np.matmul(U, np.matmul(dagger(V0), np.matmul(dagger(U), np.matmul(W0, np.matmul(U, np.matmul(V0, dagger(U))))))))

 

   

def OTOC_Test2(N, T, rep, res, slow, Op, V0):
    """
        - Purpose: Run a random circuit and use it to compute an OTOC using both the stab formalism and direct matrix multiplication.  In order to compute the OTOC one needs to introduce another operator, this is given by Op.

        - Inputs: 
               - N: integer (number of qubits).
               - T: integer (number of timesteps).
               - rep: integer (number of repetitions).
               - res: integer (resolution - i.e. how often the OTOC gets computed).
               - slow: integer (determines how much we slow down the action of the circuit).
               - Op: stim.TableauSimulator (gives a new operator which should be a super-Clifford and which we will use for computing the OTOC)
        - Outputs: - v: vector of length T/res (result of OTOC calculation using stab formalism).
                   - w: vector of length T/res (number of time steps).
                   - v2: vector of length T/res (result of OTOC calculation using matrix multiplication).

"""
    Size = int(T/res) #specify size of output vectors.
    v = np.zeros(Size)
    w = np.zeros(Size)
    v2 = np.zeros(Size)
    for m in range(Size):
        w[m] = m

    for k in range(0, rep):
        s = stim.TableauSimulator()
        U = np.kron(I, np.kron(I, I))
        for i in range(0, T):
            s, U = Test_circuit_1Step(s, U)
            if (i % res) == 0:
                tableau1: stim.Tableau = s.current_inverse_tableau()**-1
                n1 = len(tableau1)
                tableau2: stim.Tableau = Op.current_inverse_tableau()**-1
                n2 = len(tableau2)
                tableau3: stim.Tableau = s.current_inverse_tableau()
                tableau_tot: stim.Tableau =(tableau3*tableau2)*tableau1
                n = len(tableau_tot)
                zs = [tableau_tot.z_output(k) for k in range(n)]
                zs2 = np.array(zs)
                signs = [(zs[k].sign).real for k in range(n)]
                signs2 = np.array(signs)
                bMat = entropy.binaryMatrix(zs2)
                signs3 = entropy.convert_signs(signs2)
                REF, signs3 = OTOC.REF_binary(bMat, signs3, N)
                
                xs1 = OTOC.xs(REF)
                rows = entropy.rows(xs1)
                rank = entropy.gf2_rank(rows)
                size2 = N - rank
            

                small = OTOC.small_zs(REF, size2, N)

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


                W0 = XXX
                v2[i] = abs(np.trace(M(U, W0, V0)))/(2**3) 
                    
    return v, w, v2           

   

def main():

    """
    To run the below simply choose N: number of qubits, T: number of timesteps, rep: number of repititions, res: resolution. The figure should always give two lines that are precisely identical.
    """

    V0_s, V0 = Op()
    N = 3
    T = 50
    rep = 1
    res = 1
    slow = 1
    v, w, v2 = OTOC_Test2(N, T, rep, res, slow, V0_s, V0)




                     	      
    plt.plot(w, v, label = 'stim', color = 'blue')
    plt.plot(w, v2, label = 'test', color = 'red', linestyle = 'dotted')
    plt.xlabel('Time')
    plt.ylabel('OTOC')
    plt.title(f'OTOC for N = {N}')
    plt.legend()
    plt.show()  
    


 
if (__name__ == '__main__'):
    main()     
 
           


                      
