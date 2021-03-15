import numpy as np
import math
import cmath
from utils import PAULI_X, is_unitary, is_special_unitary

class Gate:
    """Represents gate acting on one qubit."""

    def __init__(self, name, arg=None):
        assert name in ['R', 'S', 'E']
        self.name = name
        self.arg = arg

    def to_matrix(self):
        if self.name == 'R':
            return np.array([[cmath.cos(self.arg / 2), -cmath.sin(self.arg / 2)],
                            [cmath.sin(self.arg / 2), cmath.cos(self.arg / 2)]])
        elif self.name == 'S':
            return np.diag([cmath.exp(-1j * self.arg), cmath.exp(1j * self.arg)])
        elif self.name == 'E':
            return np.diag([1.0, cmath.exp(1j * self.arg)])

    def is_identity(self):
        return np.allclose(self.to_matrix(), np.eye(2))

    def __repr__(self):
        if self.arg is not None:
            return self.name + "(" + str(self.arg) + ")"
        else:
            return self.name


class Gate2:
    """Represents gate acting on one qubit.
    Definitions:
    Ry(a) = exp(0.5*i*a*sigma_y)
    Rz(a) = exp(0.5*i*a*sigma_z)
    R1(a) = diag(1, exp(i*a))
    """

    def __init__(self, name, arg=None):
        assert name in ['Ry', 'Rz', 'R1', 'X']
        self.name = name
        self.arg = arg

    def to_matrix(self):
        if self.name == 'Ry':
            return np.array([[cmath.cos(self.arg / 2), cmath.sin(self.arg / 2)],
                             [-cmath.sin(self.arg / 2), cmath.cos(self.arg / 2)]])
        elif self.name == 'Rz':
            return np.diag([cmath.exp(0.5j * self.arg), cmath.exp(-0.5j * self.arg)])
        elif self.name == 'R1':
            return np.diag([1.0, cmath.exp(1j * self.arg)])
        elif self.name == 'X':
            return PAULI_X

    def is_identity(self):
        return np.allclose(self.to_matrix(), np.eye(2))

    def __repr__(self):
        if self.arg is not None:
            return self.name + "(" + str(self.arg) + ")"
        else:
            return self.name

def su_to_gates(A):
    """Decomposes 2x2 special unitary to gates Ry, Rz.
    R_k(x) = exp(0.5*i*x*sigma_k).
    """
    assert is_special_unitary(A)
    u00 = A[0, 0]
    u01 = A[0, 1]
    theta = np.arccos(np.abs(u00))
    lmbda = np.angle(u00)
    mu = np.angle(u01)

    result = []
    result.append(Gate2('Rz', lmbda - mu))
    result.append(Gate2('Ry', 2 * theta))
    result.append(Gate2('Rz', lmbda + mu))
    return result


def unitary2x2_to_gates(A):
    """Decomposes 2x2 unitary to gates Ry, Rz, R1.
    R1(x) = diag(1, exp(i*x)).
    """
    assert is_unitary(A)
    phi = np.angle(np.linalg.det(A))
    if np.abs(phi) < 1e-9:
        return su_to_gates(A)
    elif np.allclose(A, PAULI_X):
        return [Gate2('X')]
    else:
        A = np.diag([1.0, np.exp(-1j * phi)]) @ A
        return su_to_gates(A) + [Gate2('R1', phi)]
