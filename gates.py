from typing import List, Collection
import stim
import numpy as np
import random


def C3(i, j, k):
     """
     - Purpose: Create C3 gate.
     - Inputs: 
           - i (integer): first qubit to act on with C3 (i.e. control qubit).
           - j (integer): second qubit to act on with C3
           - k (integer): third qubit to act on with C3
     - Outputs:
           - c (stim.circuit): a stim circuit that applies C3 to the chosen qubits.
     """

     c = stim.Circuit()
     c.append_operation("CY", [i, j])
     c.append_operation("CY", [i, k])
 
 
     return c
    
   
    

def ZH(k):
     """
     - Purpose: Create Z.H gate.
     - Inputs: 
           - k (integer): qubit to act on with Z.H

     - Outputs:
           - c (stim.circuit): a stim circuit that applies Z.H to the chosen qubit.
     """

     c = stim.Circuit()
     c.append_operation("H", [k])
     c.append_operation("Z" , [k])
      
     return c


def SWP(i, j):
    """
     - Purpose: Create SWAP gate.
     - Inputs: 
           - i (integer): first qubit to act on with SWAP
           - j (integer): second qubit to act on with SWAP
  
     - Outputs:
           - c (stim.circuit): a stim circuit that applies SWAP to the chosen qubits.
    """
     
    c = stim.Circuit()
    c.append_operation("SWAP", [i, j])      

    return c
    
    
def Rand1Step(N, s):
    """
    - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit applies Z.H to a random qubit and applies C3 to three randomly chosen adjacent qubits.
    
    - Inputs:
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
    - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step.     
    
    """
    r1 = random.randint(0, N)
    s.do(ZH(r1))
    
    r2 = random.randint(0, N-2)
    r3 = random.randint(0,5)
    if r3 == 0:
       s.do(C3(r2, r2+1, r2+2))
    elif r3 ==1:
       s.do(C3(r2+1, r2, r2+2))
    elif r3 == 2:
       s.do(C3(r2+1, r2+2, r2))
    elif r3 == 3:
       s.do(C3(r2, r2+2, r2+1))
    elif r3 ==4:
       s.do(C3(r2+2, r2+1, r2))
    else:
       s.do(C3(r2+2, r2, r2+1))
                    
    return s
    
    
    
def FS1Step(N, s):
    """
    - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit applies C3 to half the qubits (randomly chosen and random couplings between them) and applies Z.H to the other half.
    
    - Inputs:
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
    - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step.     
    
    """
    r = [i for i in range(N)]
    random.shuffle(r)   
    
    M2 = int(N/2)
    M6 = int(N/6)
    for i in range(0, M6):
        if i == 0:
            s.do(C3(r[i], r[M6], r[2*M6]))
        else:    
            s.do(C3(r[i], r[2*i], r[3*i]))
        
    for i in range(M2, N):
        s.do(ZH(r[i]))
        
    return s




def FS2StepE(N, s):
    """
    - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit applies Z.H to all qubits on the even timesteps and then applies C3 to all the qubits (randomly coupled) on the odd timesteps.
    
    - Inputs:
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
    - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step.     
    
    """

    for i in range(0, N):
        s.do(ZH(i))
        
    
    return s
    
    
def FS2StepO(N, s):
    """
    - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit applies Z.H to all qubits on the even timesteps and then applies C3 to all the qubits (randomly coupled) on the odd timesteps.
    
    - Inputs:
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
    - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step.     
    
    """

    r = [i for i in range(N)]
    random.shuffle(r)

    M3 = int(N/3)

    s.do(C3(r[0], r[M3], r[2*M3]))

    for i in range(1, M3):
        s.do(C3(r[i], r[2*i], r[3*i]))

    return s
    
    
    
def FS3Step(N, s):
     """
     - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit applies Z.H to N/4 qubits and C3 to the remaining qubits.
    
     - Inputs:
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """


     r = [i for i in range(N)]
     random.shuffle(r)

     M4 = int(N/4)

     for i in range(3*M4, N):
          s.do(ZH(r[i]))

     for i in range(1, M4):
          s.do(C3(r[i], r[2*i], r[3*i]))

     s.do(C3(r[0], r[M4], r[2*M4]))
     return s



     return s

def FS3_NpStep(N, s, p):
     """
     - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit acts on N/p qubits each step. On those qubits which it acts on, it acts as Z.H on 1/4 of them and then as C3 on the rest.
    
     - Inputs:
         - N (number of qubits): should be divisible by 12.
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - p (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """


     r = [i for i in range(N)]
     random.shuffle(r)

     Mp = int(N/ p)
     Mq = int(N/(p*4))

     for i in range(3*Mq, Mp):
          s.do(ZH(r[i]))

     for i in range(1, Mq):
          s.do(C3(r[i], r[2*i], r[3*i]))

     s.do(C3(r[0], r[Mq], r[2*Mq]))
     return s
     


     
