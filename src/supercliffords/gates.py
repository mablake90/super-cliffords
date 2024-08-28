"""Module defining the super-clifford gates."""

import stim


def C3(i, j, k):
    """
    - Purpose: Create C3 gate.
    - Inputs:
          - i (integer): first qubit to act on with C3 (i.e. control qubit).
          - j (integer): second qubit to act on with C3
          - k (integer): third qubit to act on with C3
    - Outputs:
          - c (stim.circuit): a stim circuit that applies C3 to the chosen qubits.
    """
    c = stim.Circuit()
    c.append_operation("CY", [i, j])
    c.append_operation("CY", [i, k])
    return c


def ZH(k):
    """
    - Purpose: Create Z.H gate.
    - Inputs:
          - k (integer): qubit to act on with Z.H

    - Outputs:
          - c (stim.circuit): a stim circuit that applies Z.H to the chosen qubit.
    """
    c = stim.Circuit()
    c.append_operation("H", [k])
    c.append_operation("Z", [k])
    return c


def SWP(i, j):
    """
    - Purpose: Create SWAP gate.
    - Inputs:
          - i (integer): first qubit to act on with SWAP
          - j (integer): second qubit to act on with SWAP

    - Outputs:
          - c (stim.circuit): a stim circuit that applies SWAP to the chosen qubits.
    """
    c = stim.Circuit()
    c.append_operation("SWAP", [i, j])
    return c
