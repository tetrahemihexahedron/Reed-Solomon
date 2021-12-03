from FiniteFields import *
from itertools import combinations as combinations
from LinearSolve import *
import numpy as np

def MakeVandermondeMat(E):
    """
    Returns a matrix object that represents the Vandermonde matrix formed
    from the elements in the tuple E.

    input: E = a tuple
    """
    n = len(E) # M should be an n x n matrix
    M = []
    for e in E:
        row = []
        for i in range(n):
            row.append(e**i)
        M.append(row)
    return np.mat(M)


def FindPossSoln(T, l, F = FiniteField((1, 1, 0, 1))):
    """
    Returns a dictionary containing all the possible solutions with the number
    of times that they appear.

    input: T = a tuple representing the transmission (len(T) = the size of F)
        l = an int representing the length of the message
        F = a FiniteField object (default F_8)

    The transmission tuple T should be in the order (m(0), m(1), m(t), m(t^2),
    ..., m(t^k)) where t is the generator of the field and m(x) is the
    message polynomial. The message polynomial should be m(x) = a_0 + a_1 x +
    ... + a_(l-1)x^(l-1).

    """
    assert len(T) == F.size, "The transmission is not the right length."

    T = np.array(T)
    possSols = {}
    for c in combinations(range(len(T)), l):
        outputValues = np.mat(T[np.array(c)])
        inputValues = []
        for i in c:
            if i == 0:
                inputValues.append(0)
            else:
                inputValues.append(F.getEltFromPower(i-1))
        A = MakeVandermondeMat(inputValues)
        sol = tuple(np.array(LinearSolve(A, outputValues)).flatten())
        if sol in possSols.keys():
            possSols[sol] = possSols[sol] + 1
        else:
            possSols[sol] = 1
    return possSols

def RSDecode(T, l, F = FiniteField((1, 1, 0, 1)), showall  = False):
    possSols = FindPossSoln(T, l, F)
    if showall == True:
        solns = []
        for s in possSols.keys():
            solns += [(s, possSols[s])]
        return solns
    else:
        soln = ()
        count = 0
        for s in possSols.keys():
            if possSols[s] > count:
                soln = s
                count = possSols[soln]
        return soln


T = (t2, 0, 0, t4, t2, t, t4, t)
G = FiniteField((1, 1, 0, 0, 1)) #F_16
s = G.getGenerator()
s2 = s**2
s3 = s**3
s4 = s**4
s5 = s**5
s6 = s**6
s7 = s**7
s8 = s**8
s9 = s**9
s10 = s**10
s11 = s**11
s12 = s**12
s13 = s**13
s14 = s**14
S = (1, s5, s2, s12, s9, s5, s5, s4, 1, 0, 1, s9, s7, s14, s12, s)
c = np.array(range(10))
output = np.array(S)[c]
inpVal = [0]
for i in range(0, 9):
    inpVal = inpVal + [G.getEltFromPower(i)]
A = MakeVandermondeMat(inpVal)
