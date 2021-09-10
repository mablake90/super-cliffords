# ProjectB

Repository for keeping track of code created during project on simulations of operator scrambling. 

OTOC.py is a file which can be used to compute the OTOC calculation for super-stabilizer circuits. 
OTOC_Test.py includes code used to test the function for computing the OTOC in the OTOC.py file. 
entropy.py is a file which can be used to compute the operator entanglement entropy for super-stabilizer circuits. 
circuits.py and gates.py are files which include functions for super-stabilizer gates and functions that compute the entropy over time for these gates. They are needed for the OTOC.py and entropy.py files to work. 

All code is written in python, and the simulation of super-Clifford circuits is done using stim, a Clifford circuit simulator with python interface. This can be installed using the following command:

`pip install stim`
