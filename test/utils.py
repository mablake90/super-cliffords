import math
import random
from functools import reduce
import numpy as np
import stim
from supercliffords import gates


X = np.array([[0, 1], [1, 0]])

Id = np.array([[1, 0], [0, 1]])

Y = np.array([[0, complex(0, -1)], [complex(0, 1), 0]])

Z = np.array([[1, 0], [0, -1]])

XXX = np.kron(X, np.kron(X, X))

YXX = np.kron(Y, np.kron(X, X))

YYX = np.kron(Y, np.kron(Y, X))

YXY = np.kron(Y, np.kron(X, Y))

XYX = np.kron(X, np.kron(Y, X))

XXY = np.kron(X, np.kron(X, Y))

XYY = np.kron(X, np.kron(Y, Y))

YYY = np.kron(Y, np.kron(Y, Y))

CX12 = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]), Id)

CZ12 = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]), Id)

CX21 = np.kron(np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]), Id)

CX31 = np.array(
    [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
    ]
)


T = np.array([[1, 0], [0, (complex(1, 1)) * (1 / math.sqrt(2))]])

T_trip = np.kron(T, np.kron(T, T))

T1 = np.kron(T, np.kron(Id, Id))

T1_6 = np.linalg.matrix_power(T1, 6)

T2 = np.kron(Id, np.kron(T, Id))

T3 = np.kron(Id, np.kron(Id, T))

T2_6 = np.linalg.matrix_power(T2, 6)

SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

SWAP12 = np.kron(SWAP, Id)

SWAP23 = np.kron(Id, SWAP)

SWAP13 = np.matmul(SWAP12, np.matmul(SWAP23, SWAP12))

C3_123 = np.matmul(CX21, np.matmul(CX31, np.matmul(CZ12, np.matmul(T1_6, T2_6))))

C3_123_dag = np.conjugate(np.transpose(C3_123))

C3_213 = np.matmul(SWAP12, np.matmul(C3_123, SWAP12))

C3_132 = np.matmul(SWAP23, np.matmul(C3_123, SWAP23))

C3_312 = np.matmul(SWAP13, np.matmul(C3_132, SWAP13))

C3_321 = np.matmul(SWAP13, np.matmul(C3_123, SWAP13))

C3_231 = np.matmul(SWAP12, np.matmul(C3_132, SWAP12))

H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
S = np.array([[1, 0], [0, 1j]])
HXY = S @ X


def get_x(N):
    return reduce(lambda a, b: np.kron(a, b), [X] * N)


def dagger(U):
    """
    Takes the conjugate transpose of an operator.
    """
    return np.conjugate(np.transpose(U))


def sqrt(n):
    return math.sqrt(n)


def circuit_one_step_test(s, U, N):
    """
    Purpose: Move forward by one step of a random circuit.

    Inputs: - s: a stim circuit.
            - U: a matrix on the underlying Hilbert space, should be the same operator as that found in s.
            - N (int): number of qubits.

    Outputs: -s: stim circuit, updated by one step of random circuit.
             - U: a matrix on underlying Hilbert space, updated by one step of random circuit.
    """
    assert N > 2, "Must be three qubits or more."
    gates_dict = {}

    c = stim.Circuit()
    c.append_operation("I", [N - 1])
    c.append_operation("I", [0])
    s.do(c)

    r1 = random.randint(0, N - 1)
    s.do(gates.ZH(r1))
    gates_dict["T"] = r1
    left_dims = 2 ** (r1)
    right_dims = 2 ** (N - r1 - 1)
    left_I = np.eye(left_dims)
    right_I = np.eye(right_dims)

    U = np.matmul(U, np.kron(left_I, np.kron(T, right_I)))

    r3 = random.randint(0, 5)
    r4 = random.randint(0, N - 3)
    left_dims = 2**r4
    right_dims = 2 ** (N - r4 - 3)
    left_I = np.eye(left_dims)
    right_I = np.eye(right_dims)

    if r3 == 0:
        s.do(gates.C3(r4, r4 + 1, r4 + 2))
        gates_dict["C3"] = [r4, r4 + 1, r4 + 2]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_123, right_I)))
    elif r3 == 1:
        s.do(gates.C3(r4 + 1, r4, r4 + 2))
        gates_dict["C3"] = [r4 + 1, r4, r4 + 2]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_213, right_I)))
    elif r3 == 2:
        s.do(gates.C3(r4 + 1, r4 + 2, r4))
        gates_dict["C3"] = [r4 + 1, r4 + 2, r4]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_231, right_I)))
    elif r3 == 3:
        s.do(gates.C3(r4, r4 + 2, r4 + 1))
        gates_dict["C3"] = [r4, r4 + 2, r4 + 1]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_132, right_I)))
    elif r3 == 4:
        s.do(gates.C3(r4 + 2, r4 + 1, r4))
        gates_dict["C3"] = [r4 + 2, r4 + 1, r4]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_321, right_I)))
    else:
        s.do(gates.C3(r4 + 2, r4, r4 + 1))
        gates_dict["C3"] = [r4 + 2, r4, r4 + 1]
        U = np.matmul(U, np.kron(left_I, np.kron(C3_312, right_I)))

    return s, U, gates_dict


def op(N):
    """
    Purpose: Create an operator which will be the perturbation V(0) in the OTOC calculation.

    Inputs: - None.
    Outputs: - s: a stim circuit with operator V(0) encoded in it.
             - V0: a matrix on underlying Hilbert space representing the same operator.

    """

    s = stim.TableauSimulator()
    c = stim.Circuit()
    c.append_operation("I", [N - 1])
    c.append_operation("I", [0])
    s.do(c)
    s.do(gates.C3(0, 1, 2))

    right_dims = 2 ** (N - 3)
    V0 = np.kron(C3_123, np.eye(right_dims))

    return s, V0


def otoc_operator(U, W0, V0):
    return np.matmul(
        W0,
        np.matmul(
            U,
            np.matmul(
                dagger(V0),
                np.matmul(
                    dagger(U), np.matmul(W0, np.matmul(U, np.matmul(V0, dagger(U))))
                ),
            ),
        ),
    )


def F(U, W0, V0, N):
    return abs(np.trace(otoc_operator(U, W0, V0))) / (2**N)
