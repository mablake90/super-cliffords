### Appendix Data
The files in this directory contain all the data needed to reproduce Figure 6, from the Appendix. See `/plotters/fig6.ipynb` for a Jupyter notebook with detailed explanation and code showing how to reproduce this figure from the data.

#### Explanation of Data
This directory contains $19$ `.csv` files labelled by system sizes ranging from $N = 240, 360, ..., 2400$. Each of these files contains all the data from calculating the entropy in $100$ individual realisations of the ThreeQuarterChain circuit. In all cases, the fraction of the chain on which we compute the entanglement entropy is $m = 1/4$. In all cases the number of timesteps over which we compute the entropy is $200$.

Each individual file contains $100$ rows (one for each realisation) and $200$ columns with the value of the entropy at that time-step.
