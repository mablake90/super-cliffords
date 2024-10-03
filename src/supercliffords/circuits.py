"""Module for defining super-clifford circuits."""

import os
import numpy as np
import stim
from supercliffords.steps import (
    IdStep,
    Initialize,
    ThreeQuarterStep,
    AlternatingEven,
    AlternatingOdd,
    StepSequence,
)
from supercliffords.entropy import compute_entropy
from supercliffords.otoc import compute_otoc
from multiprocessing import Pool


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
        np.random.seed(int.from_bytes(os.urandom(4), "big"))
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

    def compute_entropy_parallel(self, t, cut, res, rep, n_jobs):
        """
        Distribute the calculation of entropy over multiple cores.

        params:
            t (int): number of timesteps.
            cut (int): The cut across which to compute the entropy.
            res (int): resolution (i.e. how often to compute the operator
            entanglement).
            rep (int): number of times to repeat the simulation and average
              over.
            n_jobs (int): number of cores to use.
        returns:
            S (np.array): Operator entanglement.
            ts (np.array): Timesteps at which the operator entanglement was
            computed.
        """
        if rep < n_jobs:
            n_jobs = rep
        ts = np.zeros(t // res)
        S = np.zeros(t // res)
        for i in range(t // res):
            ts[i] = i * res
        reps_per_job = rep // n_jobs
        args = [(t, cut, res, reps_per_job) for _ in range(n_jobs)]
        with Pool(n_jobs) as p:
            results = p.starmap(self.compute_entropy, args)
            for result in results:
                S += result[0] / n_jobs
        return S, ts

    def compute_otoc(self, t, res, rep, op):
        """
        Compute the out-of-time-ordered correlator of the circuit.
        params:
            t (int): number of timesteps.
            cut (int): The cut across which to compute the otoc.
            res (int): resolution (i.e. how often to compute the otoc).
            rep (int): number of times to repeat the simulation and average
              over.
            op (stim.TableauSimulator): The perturbation operator V0.
        returns:
            f (np.array): Out-of-time-ordered correlator.
            ts (np.array): Timesteps at which the otoc was
            computed.
        """
        np.random.seed(int.from_bytes(os.urandom(4), "big"))
        if isinstance(op, stim.TableauSimulator):
            op = op.current_inverse_tableau() ** -1
        elif not isinstance(op, stim.Tableau):
            raise ValueError(
                "op must be a stim.TableauSimulator or stim.Tableau"
            )

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

    def compute_otoc_parallel(self, t, res, rep, op, n_jobs):
        """
        Distribute the calculation of the out-of-time-ordered correlator over
        multiple cores.

        params:
            t (int): number of timesteps.
            cut (int): The cut across which to compute the entropy.
            res (int): resolution (i.e. how often to compute the otoc).
            rep (int): number of times to repeat the simulation and average
            over.
            op (stim.TableauSimulator): The perturbation operator V0.
            n_jobs (int): number of cores to use.

        returns:
            f (np.array): Out-of-time-ordered correlator.
            ts (np.array): Timesteps at which the otoc was
            computed.
        """
        if rep < n_jobs:
            n_jobs = rep
        ts = np.zeros(t // res)
        f = np.zeros(t // res)
        for i in range(t // res):
            ts[i] = i * res
        reps_per_job = rep // n_jobs

        op_tableau: stim.Tableau = op.current_inverse_tableau() ** -1
        args = [(t, res, reps_per_job, op_tableau) for _ in range(n_jobs)]
        with Pool(n_jobs) as p:
            results = p.starmap(self.compute_otoc, args)
            for result in results:
                f += result[0] / n_jobs
        return f, ts


class ThreeQuarterCircuit(Circuit):
    """
    A super-clifford circuit that acts with C3 on three quarters of the
    qubits that it acts on, and with T on the remaining quarter.

    params:
        N (int): The number of qubits in the circuit.
        slow (int): The proportion of qubits to act on at each timestep.
    """

    def __init__(self, N, slow, op_string=None):
        """
        Initialize the circuit.
        """
        steps = StepSequence(
            N,
            [
                Initialize(N, op_string),
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
