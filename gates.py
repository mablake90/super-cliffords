from typing import List, Collection
import stim
import numpy as np
import random

"""
This file contains gates which form the building blocks of the circuits that are used to compute the entropy and OTOC. The gate set we use are ZH, SWAP, C3. These gates are then organized into single timesteps of various circuits, which can then be easily run.

"""


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
    r1 = random.randint(0, N-1)
    s.do(ZH(r1))
    c = stim.Circuit()
    c.append_operation("I", [N-1])
    s.do(c) 
    
    r2 = random.randint(0, N-3)
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



def FS3_NpStep(N, s, slow):
     """
     - Purpose: Apply a circuit corresponding to one single step of a random circuit. This random circuit acts on N/p qubits each step. On those qubits which it acts on, it acts as Z.H on 1/4 of them and then as C3 on the rest.
    
     - Inputs:
         - N (number of qubits): should be divisible by 12.
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - slow (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """


     r = [i for i in range(N-1)] #Randomly chooses which qubits to act on with the gates.
     random.shuffle(r)

     Mp = int(N/slow - 1)
     
     Mq = int(N/(slow*4) -1)

     for i in range(3*Mq, Mp): #Apply ZH
          s.do(ZH(r[i]))

     for i in range(1, Mq):  #Apply C3
          s.do(C3(r[i], r[2*i], r[3*i]))

     s.do(C3(r[0], r[Mq], r[2*Mq]))


     c = stim.Circuit()
     c.append_operation("I", [N-1])
     s.do(c)
     
     return s
     

def LocInt_Step1(N, s, slow):
     """
     - Purpose: Apply a step of a random circuit with nearest neighbour interactions acting on O(N) qubits at each timestep. This random circuit acts on N/slow qubits. Step 1, ZH on N/slow qubits. Step 2, C3 on N/slow qubits. Steps 3, 4, SWAP gates on N/slow qubits. 
     - Inputs:
         - N (number of qubits).
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - slow (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """

     Choose = int(N/(slow))
     qubits = np.random.choice(N, size = Choose, replace=False) #Randomly chooses which qubits to act on with the gates.
     qubits = qubits.tolist()
     

     for i in qubits: #Act with ZH gates.
          s.do(ZH(i))
     c = stim.Circuit()
     c.append_operation("I", [N-1])
     s.do(c)
     
     return s     
  

     

def LocInt_Step2(N, s, slow):
     """
     - Purpose: Apply a step of a random circuit with nearest neighbour interactions acting on O(N) qubits at each timestep. This random circuit acts on N/slow qubits. Step 1, ZH on N/slow qubits. Step 2, C3 on N/slow qubits. Steps 3, 4, SWAP gates on N/slow qubits. 
     - Inputs:
         - N (number of qubits).
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - slow (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """



     Choose = int(N/(3*slow))

     div3 =[]
     for i in range(N):
          if i % 3 == 0:
               div3.append(i)
     div3a = np.array(div3)

     qubits = np.random.choice(div3a, size = Choose, replace = False)#Randomly chooses which qubits to act on with the gates.
     qubits = qubits.tolist()

     for i in qubits: #Act with C3 gates on the chosen qubits, randomizing which one is control.
          k = random.randint(0,6)
          if k == 0:
               s.do(C3(i, i+1, i+2))
          if k == 1:
               s.do(C3(i, i+2, i+1))
          if k == 2:
               s.do(C3(i+1, i, i+2))
          if k == 3:
               s.do(C3(i+1, i+2, i))
          if k == 4:
               s.do(C3(i+2, i+1, i))
          if k == 5:
               s.do(C3(i+2, i, i+1))
     c = stim.Circuit()
     c.append_operation("I", [N-1])
     s.do(c)
     
     return s

def LocInt_Step3(N, s, slow):
     """
     - Purpose: Apply a step of a random circuit with nearest neighbour interactions acting on O(N) qubits at each timestep. This random circuit acts on N/slow qubits. Step 1, ZH on N/slow qubits. Step 2, C3 on N/slow qubits. Steps 3, 4, SWAP gates on N/slow qubits. 
     - Inputs:
         - N (number of qubits).
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - slow (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """


     Choose = int(N/(2*slow))

     div2 = []
     for i in range(N):
          if i % 2 == 0:
               div2.append(i)
     div2a = np.array(div2)

     qubits = np.random.choice(div2a, size = Choose, replace = False) #Randomly chooses which qubits to act on with the gates.
     qubits = qubits.tolist()

     for i in qubits: #Act with SWAP gates.
          s.do(SWP(i, i+1))
     c = stim.Circuit()
     c.append_operation("I", [N-1])

     s.do(c)

     return s

def LocInt_Step4(N, s, slow):
     """
     - Purpose: Apply a step of a random circuit with nearest neighbour interactions acting on O(N) qubits at each timestep. This random circuit acts on N/slow qubits. Step 1, ZH on N/slow qubits. Step 2, C3 on N/slow qubits. Steps 3, 4, SWAP gates on N/slow qubits. 
     - Inputs:
         - N (number of qubits).
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
         - slow (integer): controls how much to slow down the circuit. (e.g. so we only act on N/p qubits at each timestep).
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step. 
     """

     

     Choose = int((N-1)/(2*slow) -1)

     ndiv2 = []
     for i in range(N-1):
          if i % 2 == 1:
               ndiv2.append(i)
     ndiv2a = np.array(ndiv2)

     qubits = np.random.choice(ndiv2a, size = Choose, replace = False) #Randomly chooses which qubits to act on with gates.
     qubits = qubits.tolist()

     for i in qubits: #Act with SWAP gates.
          s.do(SWP(i, i+1))

     c = stim.Circuit()
     c.append_operation("I", [N-1])

     s.do(c)

     return s

     
     
     
     

def Id_Step(N, s):
    """
     - Purpose: Create a Tableau that acts with identity on all qubits. Useful for ensuring tha          t nothing happens in the first step of the circuit.
     - Inputs:
         - N (number of qubits).
         - s (stim.circuit): a stim circuit to evolve forward by one time step.
     - Outputs:
         - s (stim.circuit): a stim circuit which has been evolved forward by one time step, by           acting with identity.
    """
    
     
    
    c = stim.Circuit()
    c.append_operation("I", [N-1])
    s.do(c)
    return s


