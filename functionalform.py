import numpy as np
from numpy import load
import scipy
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt




def ansatz(t, a, c):
    """
      - Purpose: Define the fitting ansatz which we will use to fit the entropy.
      - Inputs:
            - t (array): number of timesteps over which the entropy was computed.
            - a (real number): this is the parameter which we will fit.
      - Outputs:
            - S (array): Function to fit to the entropy.
    """
    cut = 60
    size = np.size(t)
    S = np.zeros(size)
    for i in range(size):
        S[i] = -math.log(2**(-cut+c) + 2**(-a*t[i]), 2)
        
    return S

 


def main():
    """
    This can be used to fit the entropy using the functional form of the entropy. The fitting ansatz for the entropy is given by the equation:

2^(-S(t)) = 2^(-cut) + 2^(-a*t)
     
Where cut is the number of spins in the subsystem we are considering. The inputs needed for this is some file which has suitable data for the entropy stored, in .npz format. 
    """
    

    data1 = load('/home/at16718/Dropbox/Dropbox_Documents/Maths/Scrambling/FS3_Entropy_M50_N360.npz')

    
    lst = data1.files
    list1 = []


    t = np.zeros(1000)
    for i in range(1000):
        t[i] = 5*i

    N = 120
    M = 50
    cut = 60
    
    for item in lst:
        list1.append(data1[item])
        array1 = data1[item]
 

    popt, pcov = curve_fit(ansatz, t,  array1)

    ansatz1 = ansatz(t, popt[0], popt[1])

    plt.plot(t, array1, color = 'blue', alpha = 0.7) 
    plt.plot(t, ansatz1, linestyle = '--', label = 'fit', color = 'red')


    plt.xlabel('Time')
    plt.ylabel('S')
    plt.title(f'Entropy, N = {N}, averaged over {M} runs, cut = {cut}.')
    plt.legend()
    plt.show()

if (__name__ == '__main__'):
    main()    

        
