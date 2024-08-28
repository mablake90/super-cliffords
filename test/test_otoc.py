import numpy as np
import stim
from supercliffords.otoc import (
    ref_binary,
    g,
    row_sum,
    xs,
    small_zs,
    compute_otoc,
)
from utils import (
    F,
    op,
    circuit_one_step_test,
    get_x,
)


def test_ref_binary():
    M = np.array([[1, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0], [1, 0, 1, 0, 0, 0]])
    signs = np.array([1, 1, 0])
    N = 3
    ref, signs = ref_binary(M, signs, N)
    assert np.allclose(
        ref, np.array([[1, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0]])
    )
    assert np.allclose(signs, np.array([1, 1, 1]))

    M = np.array([[1, 1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1]])
    signs = np.array([1, 1, 0])


def test_g():
    assert g(0, 0, 0, 1) == 0
    assert g(0, 0, 1, 0) == 0
    assert g(1, 1, 0, 1) == 1
    assert g(1, 1, 1, 0) == -1
    assert g(1, 0, 0, 1) == -1
    assert g(0, 1, 0, 1) == 0
    assert g(0, 1, 1, 1) == -1


def test_row_sum():
    h = np.array([1, 1, 0, 0])
    i = np.array([0, 0, 1, 1])
    rh = 1
    ri = 0
    assert row_sum(h, i, rh, ri, 2) == 0
    assert row_sum(i, h, ri, rh, 2) == 0

    h = np.array([1, 0, 0, 0])
    i = np.array([0, 0, 0, 1])
    rh = 0
    ri = 0
    assert row_sum(h, i, rh, ri, 1) == 0
    assert row_sum(i, h, ri, rh, 1) == 0


def test_xs():
    M = np.array([[1, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0], [1, 0, 1, 0, 0, 0]])
    expected_result = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]])
    assert np.allclose(xs(M), expected_result)

    M = np.array([[1, 1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1]])
    expected_result = np.array([[1, 1, 0], [0, 1, 0], [1, 0, 1]])
    assert np.allclose(xs(M), expected_result)


def test_small_zs():
    M = np.array([[1, 1, 0, 1, 0, 1], [0, 1, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1]])
    expected_result = np.array([[0, 0, 1], [1, 0, 1]])
    assert np.allclose(small_zs(M, 2, 3), expected_result)

    M = np.array([[1, 0, 1, 0, 0, 0], [0, 1, 0, 1, 0, 0], [1, 0, 1, 0, 0, 0]])
    expected_result = np.array([[0, 0, 0]])
    assert np.allclose(small_zs(M, 1, 3), expected_result)


def test_compute_otoc():
    N = 6
    V0_stim, V0 = op(6)
    t = 20
    U = np.eye(2**N)
    otoc_stim, otoc_exact = [], []
    W0 = get_x(N)
    s = stim.TableauSimulator()
    for _ in range(t):
        print(t)
        s, U, _ = circuit_one_step_test(s, U, N)
        otoc_stim.append(compute_otoc(s, N, V0_stim))
        otoc_exact.append(F(U, W0, V0, N))

    assert np.allclose(otoc_stim, otoc_exact)
