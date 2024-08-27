"""Module for defining the steps of a super-cliffor circuit."""
from abc import ABC, abstractmethod
import numpy as np
import stim
from supercliffords.gates import C3, ZH

class Step(ABC):
    """
    A step in a super-clifford circuit.
    params:
        N (int): The number of qubits in the circuit.
        when (str): Condition for when step should be applied.
    """
    def __init__(self, N, when):
        """
        Initialize the step.
        """
        self.N = N
        if when is None:
            when = "always"
        if isinstance(when, str) and when not in ["first", "always", "even", "odd"]:
            raise ValueError("Invalid value for 'when'. Currently accepted values are 'first', 'always', 'even','odd'.")

        self.when = when

    @abstractmethod
    def apply(self):
        pass

    def validate(self, step_count):
        when = self.when
        if when  == "always" and step_count > 0:
            return True
        elif when == "first" and step_count == 0:
            return True
        elif when == "even" and step_count % 2 == 0 and step_count > 0:
            return True
        elif when == "odd" and step_count % 2 == 1 and step_count > 0:
            return True
        else:
            return False


class IdStep(Step):
    """
    Identity step.
    """
    def __init__(self, N):
        """
        Initialize the step.
        """
        super().__init__(N, when='first')

    def apply(self, s, step_count):
        """
        Apply the step.
        """
        if self.validate(step_count):
            c = stim.Circuit()
            c.append_operation("I", [self.N-1])
            s.do(c)
        return s


class ThreeQuarterCircuit(Step):
    """
    Circuit that splits the qubits into 4 and applies C3 to 3/4 of the qubits
    being acted on, and T to the remaining qubits.
    """
    def __init__(self, N):
        """
        Initialize the step.
        """
        super().__init__(N, when='always')
    
    def apply(self, s, slow, step_count):
        """
        Apply the step.
        """
        if self.validate(step_count):
            r = [i for i in range(self.N-1)] #Randomly chooses which qubits to act on with the gates.
            np.random.shuffle(r)
            acted_on = int(self.N/slow - 1)
            quarter = int(self.N/(slow*4) -1)
            for i in range(3*quarter, acted_on): #Apply ZH
                s.do(ZH(r[i]))

            for i in range(1, quarter):  #Apply C3
                s.do(C3(r[i], r[2*i], r[3*i]))
            s.do(C3(r[0], r[quarter], r[2*quarter]))
            c = stim.Circuit()
            c.append_operation("I", [self.N-1])
            s.do(c)
        return s
    

class AlternatingEven(Step):
    """
    Circuit that acts with T on all qubits it acts on on even steps,
    and C3 on all qubits it acts on on odd steps.
    """
    def __init__(self, N):
        """
        Initialize the step.
        """
        super().__init__(N, when='even')
    
    def apply(self, s, slow, step_count):
        """
        Apply the step.
        """
        if self.validate(step_count):
            r = [i for i in range(self.N-1)] #Randomly chooses which qubits to act on with the gates.
            np.random.shuffle(r)
            acted_on = int(self.N/slow - 1)
            for i in range(acted_on):
                s.do(ZH(r[i]))
        return s

class AlternatingOdd(Step):
    """
    Circuit that acts with T on all qubits it acts on on odd steps,
    and C3 on all qubits it acts on on even steps.
    """
    def __init__(self, N):
        """
        Initialize the step.
        """
        super().__init__(N, when='odd')
    
    def apply(self, s, slow, step_count):
        """
        Apply the step.
        """
        if self.validate(step_count):
            r = [i for i in range(self.N-1)] #Randomly chooses which qubits to act on with the gates.
            np.random.shuffle(r)
            acted_on = int(self.N/slow - 1)
            third = acted_on/3
            s.do(C3(r[0], r[third], r[2*third]))
            for i in range(third):
                s.do(C3(r[i], r[2*i], r[3*third]))
        return s
