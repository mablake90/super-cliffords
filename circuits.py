from typing import List, Collection
import stim
import numpy as np
import gates
import matplotlib.pyplot as plt
import time
import math
import entropy


def runFS1(N, T, rep, res):
    """
    - Purpose: Run a fast scrambling circuit. Each time step this circuit applies C3 to half the qubits (randomly chosen and random couplings between qubits) and also applies Z.H to the remaining qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - rep: integer (Number of repetitions to average over).
          - res: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/res)
    v = np.zeros(B)
    w = np.zeros(B)

    cut = int(N/2)

    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS1Step(N, s)    
            if (i % res) == 0:
            	zs2, signs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/res)
            	v[j] += SA/rep
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*res
        
    return v, w
    
    

def runFS2(N, T, rep, res):
    """
    - Purpose: Run a fast scrambling circuit. Step 1, apply Z.H to all qubits. Step 2, apply C3 to all qubits (with completely randomized couplings between the qubits). Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - rep: integer (Number of repetitions to average over).
          - res: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/res)
    v = np.zeros(B)
    w = np.zeros(B)

    cut = int(N/2)

    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):


            if (i % 2) == 0:
                s = gates.FS2StepE(N, s)
            else:
                s = gates.FS2StepO(N, s)
                
                    
            if (i % res) == 0:
            	zs2, signs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/res)
            	v[j] += SA/rep

                
                
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*res
        
    return v, w


def runFS3(N, T, rep, res):
    """
    - Purpose: Run a fast scrambling circuit. Each time step this circuit applies C3 to 3/4 of  the qubits (randomly chosen and random couplings between qubits) and also applies Z.H to the remaining qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - rep: integer (Number of repetitions to average over).
          - res: integer (resolution i.e. how many times to compute operator entanglement)
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/res)
    v = np.zeros(B)
    w = np.zeros(B)

    cut = int(N/2)

    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS3Step(N, s)    
            if (i % res) == 0:
            	zs2, signs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/res)
            	v[j] += SA/rep
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*res
        
    return v, w




def runFS3_Np(N, T, rep, res, slow):
    """
    - Purpose: Run a fast scrambling circuit.Each time step this circuit acts on N/p qubits. On the qubits it acts on, it applies Z.H on 1/4 of the qubits and C3 on the remaining ones.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - rep: integer (Number of repetitions to average over).
          - res: integer (resolution i.e. how many times to compute operator entanglement)
          - slow: integer (controls how much to slow down the circuit, so that each timestep the circuit acts on N/slow qubits).
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/res)
    v = np.zeros(B)
    w = np.zeros(B)

    cut = int(N/2)

    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            s = gates.FS3_NpStep(N, s, slow)    
            if (i % res) == 0:
            	zs2, signs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/res)
            	v[j] += SA/rep
  
            else:
                pass
                
    for i in range(0, B):
        w[i] = i*res
        
    return v, w
    
def runLocInt(N, T, rep, res, slow, cut):
    """
    - Purpose: Run a circuit which acts on O(N) qubits in each time step but only has local interactions. Step 1, apply Z.H to N/slow qubits. Step 2, apply SWAP to N/slow qubits (with completely randomized couplings between the qubits). Step 3, apply C3 to N/slow qubits. Compute operator entanglement of this circuit.
    - Inputs:
          - N: integer (Number of qubits).
          - T: integer (Number of timesteps).
          - rep: integer (Number of repetitions to average over).
          - res: integer (resolution i.e. how many times to compute operator entanglement)
          - slow: integer (how much to slow down the action of the circuit, so only acting on N/slow qubits each timestep.
          - cut: integer (where to make the cut for the entropy computation).
    - Outputs:
          - v: array (Operator entanglement).
          - w: array (Timesteps at which the operator entanglement was computed).
    """

    B = int(T/res)
    v = np.zeros(B)
    w = np.zeros(B)


    for k in range(0, rep):
        s = stim.TableauSimulator()
        for i in range(0, T):
            j = i-1

            if (i % 2) == 0:
                s = gates.LocInt_Step2(N, s, slow)
                
            else:
                s = gates.LocInt_Step1(N,s, slow)
                
                     
            if (i % res) == 0:
            	zs2 = entropy.sample_stabilizers(s)
            	mat = entropy.binaryMatrix(zs2)
            	b2 = entropy.getCutStabilizers(mat, cut)
            	b3 = entropy.rows(b2)
            	SA = entropy.gf2_rank(b3.copy()) - cut
            	j = int(i/res)
            	v[j] += SA/rep

                  
            else:
                pass               
    for i in range(0, B):
        w[i] = i*res
        
    return v, w


