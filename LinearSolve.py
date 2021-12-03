import numpy as np
from FiniteFields import *
import shelve

def LUDecompose(A):
    """
    Returns a lower-triangular matrix L and an upper-triangular matrix U
    such that A = LU.

    input: A = a matrix object of size n x n
    output: a tuple (L, U) where L is a lower-triangular matrix and U is a
        an upper-triangular matrix satisfying A = LU

    Note: Does not permute any rows, i.e., it will not work for all square
    matrices.
    """
    assert type(A) == np.matrixlib.defmatrix.matrix
    assert A.shape[0] == A.shape[1], "The matrix " + str(A) + " is not square."

    size = A.shape[0] # the number of rows (and columns) in A
    U = A.copy()

    # initialize the lower-triangular matrix
    # Note: could use the built-in function giving the identity
    # matrix
    L = np.zeros_like(A)
    for i in range(size):
        L[i, i] = 1

    # for each column, multiply by a matrix to get zeros below the diagonal
    for i in range(size):
        # Note: could use the built-in function giving the identity
        T = np.zeros_like(A)
        for j in range(size):
            T[j, j] = 1
        for j in range(i+1, size):
            T[j, i] = -U[j, i]/U[i, i]
            # make L the product of the inverses of the matrices
            L[j, i] = -T[j, i]
        U = T*U
    return (L, U)

def ForwardSubSolve(L, b):
    """
    Returns an n x 1 matrix x such that Lx = b.

    input: L = a lower-triangular matrix of size n x n
        b = a matrix of size n x 1

    Uses forward substitution to solve the equations.
    """
    assert type(L) == np.matrixlib.defmatrix.matrix
    assert L.shape[0] == L.shape[1], "The matrix \n" + str(L) + "\n is not square."
    assert L.shape[0] == b.shape[0] and b.shape[1] == 1, "The matrix \n" + \
           str(L) + "\n and \n" + str(b) + "\n are not of the correct sizes."

    size = L.shape[0] # the number of rows (and columns) in L

    # check that L is lower-triangular
    for i in range(size):
        for j in range(i+1, size):
            if int(L[i, j]) != 0:
                raise AssertionError("The matrix \n" + str(L) + "\n is not lower"\
                                     + " triangular.")
    x = np.zeros_like(b)

    # since I can't add field elements and None, calculate the first value
    x[0, 0] = b[0, 0]/L[0, 0]

    for i in range(1, size):
        x[i, 0] = (b[i, 0] - (L[i, :i]*x[:i, 0])[0, 0])/L[i, i]
    return x

def BackwardSubSolve(U, b):
    """
    Returns an n x 1 matrix x such that Ux = b.

    input: U = an upper-triangular matrix of size n x n
        b = a matrix of size n x 1

    Uses backward substitution to solve the equations.
    """
    assert type(U) == np.matrixlib.defmatrix.matrix
    assert U.shape[0] == U.shape[1], "The matrix \n" + str(U) + "\n is not square."
    assert U.shape[0] == b.shape[0] and b.shape[1] == 1, "The matrix \n" + \
           str(U) + "\n and \n" + str(b) + "\n are not of the correct sizes."

    size = U.shape[0] # the number of rows (and columns) in L

    # check that U is upper-triangular
    for i in range(size):
        for j in range(i+1, size):
            if int(U[j, i]) != 0:
                raise AssertionError("The matrix \n" + str(U) + "\n is not upper-"\
                                     + "triangular.")
    x = np.zeros_like(b)

    # since I can't add field elements and None, calculate the first value
    x[size-1, 0] = b[size-1, 0]/U[size-1, size-1]

    for i in range(size-2, -1, -1):
        x[i, 0] = (b[i, 0] - (U[i, i+1:]*x[i+1:, 0])[0, 0])/U[i, i]
    return x

def SolveFromLUDecomposition(L, U, b):
    """
    Returns an n x 1 matrix x such that LUx = b.

    input: L = a lower-triangular matrix of size n x n
        U = an upper-triangular matrix of size n x n
        b = a matrix of size n x 1
    """
    assert type(U) == np.matrixlib.defmatrix.matrix and \
           type(L) == np.matrixlib.defmatrix.matrix and \
           type(b) == np.matrixlib.defmatrix.matrix
    assert U.shape[0] == U.shape[1], "The matrix \n" + str(U) + "\n is not square."
    assert L.shape[0] == L.shape[1], "The matrix \n" + str(L) + "\n is not square."
    assert U.shape[0] == L.shape[0] and U.shape[0] == b.shape[0]

    size = U.shape[0]

    # solve Ly = b for y
    y = ForwardSubSolve(L, b)

    # solve Ux = y for x
    x = BackwardSubSolve(U, y)

    return x

def LinearSolve(A, b):
    """
    Returns an nx1 matrix x such that Ax = b.

    input: A = a matrix of size nxn
        b = a matrix of sixe nx1 or a matrix or size 1xn
    """

    # decompDict = shelve.open('LUDecompositions.p')
    # if decompDict.has_key(str(A)):
    #    (L, U) = decompDict[str(A)]
    #    decompDict.close()
    # else:
    (L, U) = LUDecompose(A)
    #    decompDict[str(A)] = (L, U)
    #    decompDict.close()
    if b.shape == (A.shape[0], 1):
        return SolveFromLUDecomposition(L, U, b)
    if b.shape == (1, A.shape[0]):
        return SolveFromLUDecomposition(L, U, b.T)

# matrices for testing
A = np.mat([[one, zero, zero], [one, one, one], [one, t, t2]])
b = np.mat([t2, 0, 0])
B = np.mat([[one, one, one], [one, t2, t4], [one, t3, t6]])
c = np.mat([zero, t4, t2])
