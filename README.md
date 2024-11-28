
# super-cliffords


[![arXiv](https://img.shields.io/badge/arXiv-2410.19614-b31b1b.svg)](https://arxiv.org/abs/2410.19614)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![Test passing](https://github.com/mablake90/super-cliffords/actions/workflows/python-package.yml/badge.svg)

- [super-cliffords](#super-cliffords)
  - [Introduction](#introduction)
    - [Overview](#overview)
  - [Installation](#installation)
    - [Clone repository](#clone-repository)
    - [Environment](#environment)
    - [Install](#install)
    - [Test](#test)
  - [Usage](#usage)
  - [Background](#background)
  - [Data](#data)
  - [Acknowledgements](#acknowledgements)
  - [Citation](#citation)

## Introduction

This repository contains a module that can be used to compute operator entanglement entropy and out-of-time ordered correlators of super-stabilizer circuits.  This was developed as part of a research project carried out at the University of Bristol, by Anthony Thompson, Mike Blake, and Noah Linden.

The simulation techniques used are those introduced in [arXiv:2002.12824[quant-ph]](https://arxiv.org/abs/2002.12824). All code is written in python, and the simulation of stabilizer circuits is done using stim by Craig Gidney, a Clifford circuit simulator with python interface. This was introduced in [arXiv:2103.02202 [quant-ph]](https://arxiv.org/abs/2103.02202).

### Overview
The directory `/docs/examples`contains Jupyter notebooks with clear explanations on how to use this module. The directory `/data/` contains all the data needed to reproduce the plots from "Fast Scrambling in Classically Simulable Quantum Circuits.", and the directory `/docs/plotters/` contains Jupyter notebooks with code to reproduce those plots from the data. `/src/` contains the source code, and `/test/` contains unit tests.


## Installation

### Clone repository

`git clone git@github.com:anthonypetertc/super-cliffords.git`

Alternatively, [download](https://github.com/anthonypetertc/super-cliffords/archive/refs/heads/main.zip) the package and unzip it.

### Environment

`pip install numpy pandas stim pytest scipy matplotlib`

### Install

Navigate to the super-cliffords directory and install the package:

`pip install .`

### Test

To check that this has been installed correctly, run the unit tests:

`pytest`

## Usage
This module comes with two examples of circuits built from the allowed super-Clifford gates. See the directory `/docs/examples` for Jupyter notebooks with multiple examples and explanations on how to use these pre-built circuits and how to compute the operator entanglement entropy and out-of-time ordered correlators for these circuits. It is also possible to build new circuits from the same gates and to simulate these using the functions provided here.

## Background

Super-clifford circuits are circuits that act as Clifford circuits in Operator space. Concretely, we consider a system of $N$ qubits and let $S$ denote the subspace of operator space spanned by all strings of Pauli gates $\{X, Y\}$: this is a $2^N$ dimensional (Real) Hilbert space. It has previously been demonstrated that the Unitary dynamics obtained from three gates (and aribtrary compositions of them) generates a dynamics in the subspace $S$ of operator space that can be efficiently simulated. For more details on the specific gates and how they are simulated please see [this](https://arxiv.org/abs/2002.12824) reference.

Two probes of operator scrambling that can be efficiently simulated using these techniques are operator entanglement entropy, and out-of-time ordered correlators. Recall, the entanglement entropy of a state $\ket{\psi} \in \mathcal{H}$ on some region $A$ is defined as:

$$
S_A(\ket{\psi}) = -\mathrm{Tr}(\rho_A \rm log (\rho_A))
$$

Where $\rho_A$ denotes the reduced density matrix of the state $\ket{\psi}$ on region $A$.

The operator entanglement entropy of an operator $O$ is defined by associating each operator to a state in the doubled Hilbert space $\mathcal{H} \otimes \mathcal{H^*}$. The operator entanglement entropy, is simply the state entanglement entropy of this state.

An out-of-time ordered correlator is defined by:

$$
F(t) = \frac{1}{2^N} \mathrm{Tr}(W(t)^\dagger V(0)^\dagger W(t) V(0))
$$

This module includes functions that allow these quantities to be computed for super-Clifford circuits.


## Data
The directory `/data`  contains all the data that is needed to reproduce the results that appeared in "Fast Scrambling in Classically Simulable Quantum Circuits". The directory `/docs/plotters` includes Jupyter notebooks with explanations on how to plot all the figures from that paper as well as detailed explanation on how each of these sets of data was obtained. All data presented here was computed using the circuit which was explained in the paper and which is referred to here as `ThreeQuarterCircuit`.

## Acknowledgements
Anthony Thompson acknowledges support from UK Engineering and Physical Sciences Research Council  (EP/SO23607/1). Mike Blake acknowledges support from UK Research and Innovation (UKRI) under the UK government's Horizon Europe guarantee (EP/Y00468X/1).  Noah Linden gratefully acknowledges support from the UK Engineering and Physical Sciences Research Council through grants EP/R043957/1, EP/S005021/1, EP/T001062/1.

## Citation

[Fast Scrambling in Classically Simulable Quantum Circuits](https://arxiv.org/abs/2410.19614)

Bibtex:

```
@misc{blake2024fastscramblingclassicallysimulable,
      title={Fast Scrambling in Classically Simulable Quantum Circuits},
      author={Mike Blake and Noah Linden and Anthony P. Thompson},
      year={2024},
      eprint={2410.19614},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2410.19614},
}
```
