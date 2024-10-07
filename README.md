title: "super-cliffords"
header-includes:
- \usepackage{braket}

# super-cliffords

This repository contains a module that can be used to compute operator entanglement entropy and out-of-time ordered correlators of super-stabilizer circuits.  This was developed as part of a research project carried out at the University of Bristol, by Anthony Thompson, Mike Blake, and Noah Linden.

The simulation techniques used are those introduced in [arXiv:2002.12824[quant-ph]](https://arxiv.org/abs/2002.12824). All code is written in python, and the simulation of stabilizer circuits is done using stim by Craig Gidney, a Clifford circuit simulator with python interface. This was introduced in [arXiv:2103.02202 [quant-ph]](https://arxiv.org/abs/2103.02202).


## Installation

### Clone repository

`pip install git+https://github.com/anthonypetertc/super-cliffords`

Alternatively, [download](https://github.com/anthonypetertc/super-cliffords/archive/refs/heads/main.zip) the package and unzip it.

### environment

`pip install numpy, pandas, stim, pytest, scipy, matplotlib`

### Install

Navigate to the super-cliffords directory and install the package:

`pip install .`

### Test

To check that this has been installed correctly, run the unit tests:

`pytest`

## Background

Super-clifford circuits are circuits that act as Clifford circuits in Operator space. Concretely, we consider the a system of qubits and let $S$ denote the subspace of operator space spanned by all strings of Pauli gates $\{X, Y\}$: this is a $2^N$ dimensional (Real) Hilbert space. It has previously been [demonstrated](https://arxiv.org/abs/2002.12824) that the Unitary dynamics obtained from three gates (and aribtrary compositions of them) generates a dynamics in the subspace $S$ of operator space that can be efficiently simulated.

Two probes of operator scrambling that can be efficiently simulated using these techniques are operator entanglement entropy, and out-of-time ordered correlators. Recall, the entanglement entropy of a state $\ket{\psi} \in \mathcal{H}$ on some region $A$ is defined as:

$$
S_A(\ket{\psi}) = -\mathrm{Tr}(\rho_A \rm log (\rho_A))
$$

Where $\rho_A$ denotes the reduced density matrix of the state $\ket{\psi}$ on region $A$.

The operator entanglement entropy of an operator $O$ is simply defined by associating each operator to a state in the doubled Hilbert space $\mathcal{H} \otimes \mathcal{H^*}$. The operator entanglement entropy, is simply the state entanglement entropy of this state.

An out-of-time ordered correlator is defined by:

$$
F(t) = \frac{1}{2^N} \mathrm{Tr}(W(t)^\dagger V(0)^\dagger W(t) V(0))
$$

This module includes functions that allow these quantities to be computed for super-Clifford circuits. Further details on how and why it is possible to simulate these quantities using super-Clifford techniques can be found in [arXiv:2103.02202 [quant-ph]](https://arxiv.org/abs/2103.02202).
