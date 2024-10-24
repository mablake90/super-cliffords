### Entropy Data
The files in this directory contain all the data needed to reproduce Figures 1-3. See `/plotters/fig1.ipynb`, `/plotters/fig2.ipynb` , `/plotters/fig3.ipynb` for Jupyter notebooks with detailed explanations and code showing how to reproduce these figures from the data.


#### Explanation of data

This directory contains 24 `.npz`files, labelled by system size $N= 240, 360, ..., 3000$. Each of these files contains the results from calculating the entropy on a system with $N$ qubits, using the ThreeQuarterCircuit. In all cases the entropy has been averaged over $500$ realisations of the circuit, to smooth out fluctuations and in all cases the fraction of the total number of qubits that we use to calculate the entropy is $m=1/4$. The number of timesteps over which we compute the entropy is $200$.

Each `.npz` file contains a single array of shape $(1, 200)$ with the average entropy over all $500$ realisations.
