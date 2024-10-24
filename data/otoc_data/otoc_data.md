### OTOC Data
The files in this directory contain all the data needed to reproduce Figures 4-5. See `/plotters/fig4.ipynb`, `/plotters/fig5.ipynb` for Jupyter notebooks with detailed explanations and code showing how to reproduce these figures from the data.

#### Explanation of data

##### Figure 4:
The data to reproduce Figure 4 can be found in the directory `/otoc_data/fig4/`.

The left panel of Figure 4, requires the files `/data/otoc_data/fig4/fig4_left1.npz` and `/data/otoc_data/fig4/fig4_left2.npz`. These contain data from simulations of the out-of-time ordered correlator, on a system of $120$ qubits, averaged over $1000$ realisations of the circuit. The operator $V(0) = C3$ and the initial operators $W(0)$ were of the form:

 $$
D_1 XX...X
 $$

 Where $D_1 = X$ for the data in `/data/otoc_data/fig4/fig4_left1.npz` and $D_1 = Y$ for the data in `/data/otoc_data/fig4/fig4_left2.npz`.

 Each of these files store a single array of shape $(1, 200)$ containing the averaged otoc values for $200$ timesteps.

The data for reproducing the right panel of figure 4, can be found in `/data/otoc_data/fig4/fig4_right.csv`, it contains the data for OTOC calculations for $2000$ realisations of a circuit with $120$ qubits, these have not been averaged and instead all of the individual results for the calculations are included. The data found here computed the OTOC where $V(0) = C3$ and $W(0) = XX...X$. The OTOC was computed over $200$ timesteps.

This file contains $2000$ rows and $200$ columns. Each row corresponds to the results from a single realisation of the circuit, and each column corresponds to the value of the OTOC at that timestep.

##### Figure 5
The data for reproducing the left panel of figure 5, can be found in `/data/otoc_data/fig5/fig5_left`, this contains two files: `fig5_leftX.csv` and `fig5_leftY.csv`. These include data for the OTOC calculated on a system of $1000$ qubits over $250$ individual realisations over $200$ timesteps, with $V(0) = T_3 C3$ and $W(0)$ of the form:


 $$
D_1 XX...X
 $$

 Where $D_1 = X$ for the data in `/data/otoc_data/fig5/fig5_left/fig5_leftX.csv` and $D_1 = Y$ for the data in `/data/otoc_data/fig5/fig5_left/fig5_leftY.csv`.

 Each of these files contains a single array of shape $(1, 200)$ containing the averaged value of the OTOC over $250$ realisations.

 The data for reproducing the right panel of Figure 5 can be found in `/data/otoc_data/fig5/fig5_right`, this contains 5 different `.npz`files labelled by system-size and each of these contains the data for the out-of-time-ordered correlator calculated for $N = 120, 240, 360, 480$, each of which was calculated and averaged over $250$ realisations. In all cases the operator $V(0) = C3$ and $W(0) = XXX...X$, and the OTOC was computed over $200$ timesteps.

 Each of these $4$ `.npz`files contains a single array of shape $(1, 200)$ containing the averaged values of the OTOC over $250$ realisations.
