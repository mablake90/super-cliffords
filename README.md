# super-cliffords

This repository contains code used to compute operator entanglement entropy and out-of-time ordered correlators of super-stabilizer circuits. This formed part of a project carried out at the University of Bristol, with Mike Blake and Noah Linden. The repository is primarily for my own personal use.

The simulation techniques used are those introduced in [arXiv:2002.12824[quant-ph]]. All code is written in python, and the simulation of stabilizer circuits is done using stim by Craig Gidney, a Clifford circuit simulator with python interface. This was introduced in [arXiv:2103.02202 [quant-ph]]. It can be installed by running the command:


`pip install stim`


The files in this repository have the following functions:

1. OTOC.py is a file which can be used to compute the OTOC calculation for super-stabilizer circuits.

2. OTOC_Test.py includes code used to test the function for computing the OTOC in the OTOC.py file.

3. entropy.py is a file which can be used to compute the operator entanglement entropy for super-stabilizer circuits.

4. circuits.py and gates.py are files which include functions for super-stabilizer gates and functions that compute the entropy over time for these gates. They are needed for the OTOC.py and entropy.py files to work.
