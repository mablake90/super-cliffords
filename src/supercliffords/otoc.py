import stim
import numpy as np
import supercliffords.entropy as entropy

"""
The code in this file allows one to compute the OTOC of certain random
 circuits.

The function row_sum and g are taken directly from [arXiv:quant-ph/0406196].
"""


def ref_binary(A, signs, N):
    """
    - Purpose: Given a N x 2N matrix, A and an array of signs. This will
      convert the matrix to row echelon form (REF) and convert the signs using
        the rowsum operation.
    - Inputs:
            - A (binary N x 2N np.ndarray).
            - signs (np.ndarray of length N)
            - N (integer).
    - Outputs:
            - A (binary N x 2N np.ndarray - in REF).
            - signs (array of length N - updated using rowsum operation).

    """

    n_rows, n_cols = A.shape
    assert n_cols == 2 * n_rows, "Matrix must be of shape (N, 2N)"

    # Compute row echelon form (REF)
    current_row = 0
    for j in range(n_cols):  # For each column
        if current_row >= n_rows:
            break

        pivot_row = current_row

        # find the first row in this column with non-zero entry.
        # this becomes the pivot row
        while pivot_row < n_rows and A[pivot_row, j] == 0:
            pivot_row += 1

        # if we reach the end, this column cannot be eliminated.
        if pivot_row == n_rows:
            continue

        # otherwise, swap current row with the pivot row
        A[[current_row, pivot_row]] = A[[pivot_row, current_row]]
        a = signs[current_row]
        signs[current_row] = signs[pivot_row]
        signs[pivot_row] = a

        pivot_row = current_row
        current_row += 1

        # Eliminate rows below
        for i in range(current_row, n_rows):
            # subtract current row from any other rows beneath with
            # a non-zero entry in this current column
            if A[i, j] == 1:
                A[i] = (
                    A[i] + A[pivot_row]
                ) % 2  # subtracting is same as adding in GF(2)
                signs[i] = row_sum(
                    A[i], A[pivot_row], signs[i], signs[pivot_row], N
                )

    return A, signs


def g(x1, z1, x2, z2):
    """
    Purpose: Computes the function g needed for the rowsum operation.
    Inputs:
         - x1, z1, x2, z2  in {0, 1} (i.e. four bits).
    Outputs:
         - g in {-1, 0, 1}
    """
    if (x1 == 0) and (z1 == 0):
        g = 0
    if (x1 == 0) and (z1 == 1):
        g = x2 * (1 - 2 * z2)
    if (x1 == 1) and (z1 == 0):
        g = z2 * (2 * x2 - 1)
    if (x1 == 1) and (z1 == 1):
        g = z2 - x2

    return g


def row_sum(h, i, rh, ri, N):
    """
    Purpose: Compute the row_sum operation.
    Inputs:
         - h, i -  two arrays of bits, length 2N
           (two rows from a binary matrix).
         - rh, ri - two bits (the signs corresponding to the rows, h, i).
         - N - an integer. The number of qubits in the chain.
    Outputs:
         - rh - a bit.
         The sign of the new row obtained from adding h + i as bitstrings.
    """
    k = 0
    for j in range(N):
        k += g(i[j], i[j + N], h[j], h[j + N])
    f = 2 * rh + 2 * ri + k
    if (f % 4) == 0:
        rh = 0
    elif (f % 4) == 2:
        rh = 1
    else:
        raise ValueError("Error in row_sum operation")
    return rh


def xs(binary_array):
    """
    Purpose: Given a stabilizer tableau (N, 2N) extract the X's of
    the tableau (i.e. the first N columns).
    Inputs:
         - bin_array (N, 2N), np.ndarray, binary - a stabilizer tableau
            ( i.e. an (N, 2N) binary matrix).
    Outputs:
         - xs (np.ndarray) - a binary matrix of size (N, N).
    """
    N = len(binary_array)
    xs = np.zeros((N, N))
    xs[:, :] = binary_array[:, :N]
    return xs


def small_zs(bin_array, starting_row_index, N):
    """
    Purpose: Given a stabilizer tableau (N, 2N), extract
        a small portion of the Z's.
      I.e. the second N columns and only the rows with index larger
        than N-starting_row_index.
    Inputs:
          - bin_array (np.ndarray) - a (N, 2N) binary matrix.
          - starting_row_index (int)- an integer, less than N.
          - N (int) - an integer, number of qubits.
    Outputs:
          - small_zs (np.ndarray)- a (N-starting_row_index, N) binary matrix.

    """
    small_zs = np.zeros((starting_row_index, N))
    small_zs[:, :] = bin_array[N - starting_row_index :, N:]
    return small_zs


def prepare_op_string(op_string, N):
    """
    Purpose: Prepare the operator string for the OTOC calculation.
    Inputs:
         - op_string (str or None) - a string of length N, containing only "X" and "Y".
         - N (int) - an integer, number of qubits.
    Outputs:
         - op_tableau (stim.TableauSimulator) - the operator in tableau form.
    """
    if op_string is None:
        op_string = "X" * N
    if isinstance(op_string, str):
        raise ValueError("op_string must be a string")
    return op_string


def compute_otoc(s, N, op_tableau):
    """
    Purpose: Compute the OTOC of a given circuit.
    Inputs:
         - s (stim.TableauSimulator) - the circuit.
         - op_tableau (stim.Tableau) - the operator.
         - N (int) - the number of qubits.

    Outputs:
         - otoc (float) - the out-of-time-order correlator.
    """
    tableau1: stim.Tableau = s.current_inverse_tableau() ** -1
    tableau3: stim.Tableau = s.current_inverse_tableau()
    tableau_tot: stim.Tableau = (tableau3 * op_tableau) * tableau1
    n = len(tableau_tot)
    zs = [tableau_tot.z_output(k) for k in range(n)]
    zs_array = np.array(zs)
    signs = np.array([(zs[k].sign).real for k in range(n)])
    bin_arr = entropy.binary_matrix(zs_array)
    converted_signs = entropy.convert_signs(signs)

    ref, converted_signs = ref_binary(bin_arr, converted_signs, N)
    x = xs(ref)
    rows = entropy.rows(x)
    rank = entropy.gf2_rank(rows)
    starting_rows = N - rank
    reduced_signs = converted_signs[rank:]

    if any(reduced_signs[:starting_rows]) == 1:
        return 0
    else:
        return 2 ** (-rank / 2)
