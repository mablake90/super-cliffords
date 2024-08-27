import numpy as np
import stim
from supercliffords.entropy import sample_stabilisers, binary_matrix
from supercliffords.gates import ZH

def test_sample_stabilisers():
    s = stim.TableauSimulator()
    s.do(ZH(0))
    s.do(ZH(1))
    s.do(ZH(2))
    zs = sample_stabilisers(s)
    bin_mat = binary_matrix(zs)
    assert bin_mat.shape == (3, 6)
    assert np.allclose(bin_mat[0:3, 0:3], np.eye(3))
    assert np.allclose(bin_mat[0:3, 3:6], np.zeros((3, 3)))
