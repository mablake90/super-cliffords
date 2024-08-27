import numpy as np
import stim
from supercliffords.gates import C3, ZH
from utils import T, C3_123, X, I, HXY, XXX

def test_C3():
    s = stim.TableauSimulator()
    s.do(C3(0,1,2))
    z1 = stim.PauliString("ZII")
    z2 = stim.PauliString("IZI")
    z3 = stim.PauliString("IIZ")
    x1 = stim.PauliString("XII")
    x2 = stim.PauliString("IXI")
    x3 = stim.PauliString("IIX")
    obs = [z1, z2, z3, x1, x2, x3]
    z1u = np.kron(np.kron(X, I), I)
    z2u = np.kron(np.kron(I, X), I)
    z3u = np.kron(np.kron(I, I), X)
    x1u = np.kron(np.kron(HXY, I), I)
    x2u = np.kron(np.kron(I, HXY), I)
    x3u = np.kron(np.kron(I, I), HXY)
    obsu = [z1u, z2u, z3u, x1u, x2u, x3u]
    ket0 = C3_123.conjugate().transpose() @ XXX @ C3_123
    for ob, obu in zip(obs, obsu):
        res = s.peek_observable_expectation(ob)
        resu = 1/8 * np.trace(ket0.conjugate().transpose() @ obu.conjugate().transpose() @ ket0 @ obu).item()
        assert np.isclose(res, resu)


def test_ZH():
    s = stim.TableauSimulator()
    s.do(ZH(0))
    z1 = stim.PauliString("Z")
    x1 = stim.PauliString("X")
    obs = [z1, x1]
    z1u = X
    x1u = HXY
    obsu = [z1u, x1u]
    ket0 = T.conjugate().transpose() @ X @ T
    for ob, obu in zip(obs, obsu):
        print(ob)
        print(obu)
        res = s.peek_observable_expectation(ob)
        resu = 1/2 * np.trace(ket0.conjugate().transpose() @ obu.conjugate().transpose() @ ket0 @ obu).item()
        print(res, resu)
        assert np.isclose(res, resu)


    

    