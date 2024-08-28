import numpy as np
import stim
from supercliffords.entropy import (
    sample_stabilisers,
    binary_matrix,
    convert_signs,
    get_cut_stabilizers,
    rows,
    gf2_rank,
    compute_entropy,
)
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


def test_convert_signs():
    signs = np.array([1, 1, 0, 1, 0])
    converted_signs = np.array([0, 0, 1, 0, 1])
    assert np.allclose(convert_signs(signs), converted_signs)


def test_get_cut_stabilizers():
    matrix = np.array(
        [[1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0], [0, 0, 1, 0, 0, 1]]
    )
    cut = 2
    cut_matrix = get_cut_stabilizers(matrix, cut)
    expected_cut_matrix = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 0, 0]])
    assert np.allclose(cut_matrix, expected_cut_matrix)


def test_rows():
    matrix = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 0, 0]])
    print(rows(matrix))
    assert np.allclose(rows(matrix), [5, 10, 0])

    matrix = np.array([[1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]])
    assert np.allclose(rows(matrix), [5, 10, 4])

    matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    assert np.allclose(rows(matrix), [1, 2, 4])


def test_gf2_rank():
    assert gf2_rank([5, 10, 15]) == 2
    assert gf2_rank([5, 5, 5]) == 1
    assert gf2_rank([0, 0, 0]) == 0
    assert gf2_rank([1, 2, 4]) == 3


def test_compute_entropy():
    # bell pair test.
    s = stim.TableauSimulator()
    c = stim.Circuit()
    c.append_operation("H", [0])
    c.append_operation("CX", [0, 1])
    s.do(c)
    assert compute_entropy(s, 1) == 1

    # product state test.
    sprod = stim.TableauSimulator()
    c = stim.Circuit()
    c.append_operation("H", [0])
    c.append_operation("H", [1])
    sprod.do(c)
    assert compute_entropy(sprod, 1) == 0

    # GHZ test.
    sghz = stim.TableauSimulator()
    c = stim.Circuit()
    c.append_operation("H", [0])
    c.append_operation("CX", [0, 1])
    c.append_operation("CX", [1, 2])
    sghz.do(c)
    assert compute_entropy(sghz, 1) == 1
    assert compute_entropy(sghz, 2) == 1
