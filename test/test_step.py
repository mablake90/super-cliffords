import numpy as np
import stim
from supercliffords.gates import C3, ZH
from supercliffords.steps import Step, IdStep

class StepT(Step):
    def __init__(self, N, when):
        super().__init__(N, when)
    
    def apply(self, s, step_count):
        return s


def test_validate():
    i = np.random.randint(1, 1000)
    assert StepT(3, "always").validate(i)
    assert StepT(3, "first").validate(0)
    assert StepT(3, "even").validate(2*i)
    assert StepT(3, "odd").validate(2*i+1)
    assert not StepT(3, "first").validate(i)
    assert not StepT(3, "even").validate(2*i+1)
    assert not StepT(3, "odd").validate(2*i)
    assert not StepT(3, "always").validate(0)
    assert not StepT(3, "even").validate(0)
    assert not StepT(3, "odd").validate(0)


def test_IdStep():
    step = IdStep(3)
    s = stim.TableauSimulator()
    s3 = step.apply(s, 0)
    assert s3.peek_observable_expectation(stim.PauliString("ZII")) == 1
    assert s3.peek_observable_expectation(stim.PauliString("IZI")) == 1
    assert s3.peek_observable_expectation(stim.PauliString("IIZ")) == 1
    assert s3.peek_observable_expectation(stim.PauliString("XII")) == 0
    assert s3.peek_observable_expectation(stim.PauliString("IXI")) == 0
    assert s3.peek_observable_expectation(stim.PauliString("IIX")) == 0

    
