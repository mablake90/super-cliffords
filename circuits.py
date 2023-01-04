from typing import List, Collection
import stim
import numpy as np
import gates
import matplotlib.pyplot as plt
import time
import math
import entropy



def runFS1(N, T, M, l, cut):
    """
    - Purpose: Run a fast scrambling circuit. Each time step this circuit applies C3 to half the qubits (randomly chosen and random couplings between qubits) and also applies Z.H to the remaining qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - k: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS1Step(N, s)    
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w
    
    

def runFS2(N, T, M, l, cut):
    """
    - Purpose: Run a fast scrambling circuit. Step 1, apply Z.H to all qubits. Step 2, apply C3 to all qubits (with completely randomized couplings between the qubits). Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - k: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        for i in range(0, T):


            if (i % 2) == 0:
                s = gates.FS2StepE(N, s)
            else:
                s = gates.FS2StepO(N, s)
                
                    
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M

                
                
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w


def runFS3(N, T, M, l, cut):
    """
    - Purpose: Run a fast scrambling circuit. Each time step this circuit applies C3 to 3/4 of  the qubits (randomly chosen and random couplings between qubits) and also applies Z.H to the remaining qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - k: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS3Step(N, s)    
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w




def runFS3_Np(N, T, M, l, p, cut):
    """
    - Purpose: Run a fast scrambling circuit.Each time step this circuit acts on N/p qubits. On the qubits it acts on, it applies Z.H on 1/4 of the qubits and C3 on the remaining ones.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - l: integer (resolution i.e. how many times to compute operator entanglement)
          - p: integer (controls how much to slow down the circuit, so that each timestep the circuit acts on N/p qubits).
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        for i in range(0, T):
            if i == 0:
                s = gates.Id_Step(N, s)
            else:    
                s = gates.FS3_NpStep(N, s, p)
                
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w


def runLocInt(N, T, M, l, slow, cut):
    """
    - Purpose: Run a circuit which acts on O(N) qubits in each time step but only has local interactions. Step 1, apply Z.H to N/slow qubits. Step 2, apply SWAP to N/slow qubits (with completely randomized couplings between the qubits). Step 3, apply C3 to N/slow qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - l: integer (resolution i.e. how many times to compute operator entanglement)
          - slow: integer (how much to slow down the action of the circuit, so only acting on N/slow qubits each timestep.
          - cut: integer (where to make the cut for the entropy computation).
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        for i in range(0, T):
            j = i-1

            if i == 0:
                s = gates.Id_Step(N, s)
            else:
                s = gates.LocInt_Step1(N, s, slow)
                s = gates.LocInt_Step2(N, s, slow)

            if i % 2 == 0:
                s = gates.LocInt_Step3(N, s, slow)
            else:
                s = gates.LocInt_Step4(N, s, slow)


                
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M

                  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w


def runLocal_Scrambling(N, T, M, l, cut):
    """
    - Purpose: Run a circuit which acts on O(N) qubits in each time step but only has local interactions. Step 1, apply Z.H to N/slow qubits. Step 2, apply SWAP to N/slow qubits (with completely randomized couplings between the qubits). Step 3, apply C3 to N/slow qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - M: integer (Number of repetitions to average over).
          - l: integer (resolution i.e. how many times to compute operator entanglement)
          - cut: integer (where to make the cut for the entropy computation).
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/l)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, M):
        s = stim.TableauSimulator()
        s = gates.Local_Scrambling1_Init(N, s)

 
        
        for i in range(0, T):
            s = gates.Local_Scrambling1(N,s)
                
                     
            if (i % l) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/l)
            	v[j] += SA/M

                  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*l
        
    return v, w


    
