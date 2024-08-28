"""Module for defining super-clifford circuits."""

import numpy as np
import stim
from supercliffords.steps import (
    IdStep,
    ThreeQuarterStep,
    AlternatingEven,
    AlternatingOdd,
    StepSequence,
)
from supercliffords.entropy import compute_entropy
from supercliffords.otoc import compute_otoc


class Circuit:
    """
    A super-clifford circuit.
    params:
        N (int): The number of qubits in the circuit.
        steps (supercliffords.StepSequence): The steps of the circuit.
    """

    def __init__(self, N, steps):
        """
        Initialize the circuit.
        """
        self.N = N
        self.steps = steps

    def compute_entropy(self, t, cut, res, rep):
        """
        Compute the entropy of the circuit.
        params:
            t (int): number of timesteps.
            cut (int): The cut across which to compute the entropy.
            res (int): resolution (i.e. how often to compute the operator
            entanglement).
            rep (int): number of times to repeat the simulation and
            average over.
        returns:
            S (np.array): Operator entanglement.
            ts (np.array): Timesteps at which the operator entanglement was
            computed.
        """
        ts = np.zeros(t // res)
        S = np.zeros(t // res)
        for i in range(t // res):
            ts[i] = i * res

        for _ in range(rep):
            s = stim.TableauSimulator()
            for stepcount in range(0, t):
                s = self.steps.apply(s, stepcount)
                if stepcount % res == 0:
                    S[stepcount // res] += compute_entropy(s, cut) / rep
        return S, ts

    def compute_otoc(self, t, res, rep, op):
        """
        Compute the out-of-time-ordered correlator of the circuit.
        params:
            t (int): number of timesteps.
            cut (int): The cut across which to compute the entropy.
            res (int): resolution (i.e. how often to compute the operator
            entanglement).
            rep (int): number of times to repeat the simulation and average
              over.
            op (stim.TableauSimulator): The perturbation operator V0.
        returns:
            f (np.array): Out-of-time-ordered correlator.
            ts (np.array): Timesteps at which the operator entanglement was
            computed.
        """
        ts = np.zeros(t // res)
        f = np.zeros(t // res)
        for i in range(t // res):
            ts[i] = i * res

        for _ in range(rep):
            s = stim.TableauSimulator()
            for stepcount in range(0, t):
                s = self.steps.apply(s, stepcount)
                if stepcount % res == 0:
                    f[stepcount // res] += compute_otoc(s, self.N, op) / rep
        return f, ts


class ThreeQuarterCircuit(Circuit):
    """
    A super-clifford circuit that acts with C3 on three quarters of the
    qubits that it acts on, and with T on the remaining quarter.

    params:
        N (int): The number of qubits in the circuit.
        slow (int): The proportion of qubits to act on at each timestep.
    """

    def __init__(self, N, slow):
        """
        Initialize the circuit.
        """
        steps = StepSequence(
            N,
            [
                IdStep(N),
                ThreeQuarterStep(N, slow),
            ],
        )
        super().__init__(N, steps)


class AlternatingCircuit(Circuit):
    """
    A super-clifford circui that acts with T on all qubits it acts on on even
      steps,
    and with C3 on all qubits it acts on on odd steps.

    params:
        N (int): The number of qubits in the circuit.
        slow (int): The proportion of qubits to act on at each timestep.
    """

    def __init__(self, N, slow):
        """
        Initialize the circuit.
        """
        steps = StepSequence(
            N,
            [
                IdStep(N),
                AlternatingEven(N, slow),
                AlternatingOdd(N, slow),
            ],
        )
        super().__init__(N, steps)
