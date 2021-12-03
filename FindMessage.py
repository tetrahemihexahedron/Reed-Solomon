from FiniteFields import *
from LinearSolve import *
from RS import MakeVandermondeMat
import numpy as np
from numpy.polynomial import Polynomial as Poly
import binascii
import random

def ConvertSol(sol):
    """
    Finds the alphanumeric message whose message polynomial corresponds
    to sol.

    input: sol = tuple of coefficients of the message polynomial
    """
    binstr = ''
    for coeff in sol:
        binstr += str(coeff)
    return binascii.unhexlify(hex(int(binstr, 2))[2:])

def FindMess(I, O):
    """
    Finds the coefficients of the message polynomial with the required outputs.

    input: I = a tuple representing the values at which the message
        polynomial is evaluated
        O = a tuple representing the result of the evaluation
    """
    inputs = list(I)
    outputs = np.mat(np.array(O))
    A = MakeVandermondeMat(inputs)
    sol = tuple(np.array(LinearSolve(A, outputs)).flatten())
    message = ConvertSol(sol)
    return message

def FindGoodWord(I, partO, num = 10):
    """
    Pads partO with additional field elements and finds the corresponding
    messages.

    input: I = a tuple representing the values of the inputs
    partO: a tuple representing the required outputs
    """
    numInputs = len(I)
    currNumOutputs = len(partO)
    # pad the outputs with random field elements
    messages = []
    for i in range(num):
        extra = []
        for i in range(numInputs - currNumOutputs):
            extra.append(random.choice(range(16)))
        extraOutputs = []
        for i in extra:
            if i == 0:
                extraOutputs.append(0)
            else:
                extraOutputs.append(G.getEltFromPower(i - 1))
        outputs = partO + tuple(extraOutputs)
        messages.append((outputs, FindMess(I, outputs)))
    return messages

def getListAllFieldElts(F = FiniteField((1, 1, 0, 0, 1))):
    fieldElts = [0]
    for i in range(1, F.getSize()):
        fieldElts.append(F.getEltFromPower(i - 1))
    return fieldElts

def findTrans(m, F = FiniteField((1, 1, 0, 0, 1))):
    trans = []
    for i in getListAllFieldElts(F):
        trans.append(m(i))
    return trans

def findCoeff(s):
    """
    input: s = word to be sent
    """
    allCoeff = []
    for letter in s:
        letterHex = binascii.hexlify(letter)
        coeffOne = "{0:04b}".format(int(letterHex[0], 16))
        allCoeff.append(coeffOne)
        coeffTwo = "{0:04b}".format(int(letterHex[1], 16))
        allCoeff.append(coeffTwo)
    return allCoeff

G = FiniteField((1, 1, 0, 0, 1))
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
I = (0, 1, s, s2, s3, s4, s5, s6, s7, s8)
O = (s2, s5, s10, 1, s11, s13, s9, s9, 0, s12)
WantedO = (s2, s4, s2, s4, s8, s5)
WantedI = (s4, s7, s8, s11, s12, s14)
