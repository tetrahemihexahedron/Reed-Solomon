import numpy as np
from numpy.polynomial import Polynomial as Poly

class FiniteField(object):
    """
    Creates a FiniteField object representing a field of size 2^k.

    Usage: FiniteField((a0, a1, a2,..., ak))
        The tuple (a0, a1, a2,..., ak) represents an irreducible polynomial,
        a0 + a1*x + a2*x^2 + ... + ak*x^k, in Z_2[x].

    Does not check whether the polynomial is irreducible!
    """
    def __init__(self, irrPoly):
        # Note: Does not check whether the polynomial is irreducible
        self.irrPoly = Poly(np.array(irrPoly))
        self.size = 2**self.irrPoly.degree()
        # the generator is set as x
        self.generator = FFieldElt(self, [0, 1] + [0]*(self.irrPoly.degree() - 2))
        self.lookUp = {} # non-zero field elements as power of the generator
        for i in range(self.size - 1):
            self.lookUp[i] = self.generator**i
            self.lookUp[self.generator**i] = i

    def getSize(self):
        return self.size

    def getIrrPoly(self):
        return self.irrPoly

    def getGenerator(self):
        return self.generator

    # Could add a method to change the way the field is represented
    # and a way to change the generator

    def getPowerOfGenerator(self, elt):
        assert (type(elt) == FFieldElt and self == elt.getField())\
               or elt == 1 or elt == 0, str(elt) + " is not in the field " \
               + str(self) + "."
        if (type(elt) == FFieldElt and elt.iszero()) or elt == 0:
            raise AssertionError("The element 0 is not a power of the \
                                 generator.")

        if elt == 1 or elt.isone():
            return 0
        else:
            return self.lookUp[elt]

    def getEltFromPower(self, exp):
        assert type(exp) == int
        if exp == 1:
            return self.generator
        else:
            return self.lookUp[exp % (self.size - 1)]

    def __eq__(self, other):
        return self.irrPoly == other.irrPoly

    def __ne__(self, other):
        return self.irrPoly != other.irrPoly

class FFieldElt(object):
    """
    Creates a FFieldElt object representing an element of a finite field.

    Usage: FFieldElt(field, [a0, a1, a2,..., an])
        field: a FiniteField object (field.getSize() = 2*k)
        [a0, a1, a2,..., an]: the list of coefficients of the element
            a0 + a1*x + a2*x^2 + ... + an*x^n, n = k-1
    """
    def __init__(self, field, poly):
        assert type(field) == FiniteField
        if 2**len(poly) != field.getSize():
            raise ValueError('The finite field element has the wrong number of\
                             coefficients.')
        self.poly = Poly(poly)
        self.field = field

    def getField(self):
        return self.field

    def __str__(self):
        """
        Prints the FFieldElt as string of bits starting with the coefficient
        of the highest degree term.
        """
        name = ''
        for i in reversed(self.poly.coef):
            name += str(int(i))
        return name

    def __repr__(self):
        name = ''
        for i in reversed(self.poly.coef):
            name += str(int(i))
        return name

    def iszero(self):
        return (self.poly.coef != 0).nonzero()[0].size == 0

    def isone(self):
        return (self.poly.coef[0] == 1 and
                (self.poly.coef != 0).nonzero()[0].size == 1)

    def __eq__(self, other):
        return (type(other) == FFieldElt and self.field == other.field\
                and self.poly == other.poly) or (self.iszero() and other == 0)\
                or (self.isone() and other == 1)

    def __ne__(self, other):
        return not(self == other)

    def __hash__(self):
        # TODO: create a better hash function--the value of the binary string
        # the number of 1's
        return (self.poly.coef != 0).nonzero()[0].size

    def __add__(self, other):
        assert (type(other) == FFieldElt and self.field == other.field) \
               or other == 0 or abs(other) == 1 ,"You can't add " + str(self)\
               + " and " + str(other) + "."
        if type(other) == FFieldElt:
            newPoly = []
            for i in range(len(self.poly.coef)):
                newPoly += [int(self.poly.coef[i]) ^ int(other.poly.coef[i])]
            return FFieldElt(self.getField(), newPoly)
        if other == 0:
            return self
        if abs(other) == 1:
            newPoly = self.poly.coef.tolist()
            newPoly[0] = int(newPoly[0]) ^ 1
            return FFieldElt(self.field, newPoly)

    def __radd__(self, other):
        assert other == 0 or abs(other) == 1, "You can't add " + \
               str(other) + " and " + str(self) + "."
        if other == 0:
            return self
        if abs(other) == 1:
            return self + other

    def __sub__(self, other):
        """
        For finite fields with characteristic 2, subtraction is the same as
        addition.
        """
        return self + other

    def __rsub__(self, other):
        return other + self

    def __mul__(self, other):
        assert (type(other) == FFieldElt and self.field == other.field) \
               or other == 0 or abs(other) == 1, "You can't multiply " + \
               str(self) + " and " + str(other) + "."
        if type(other) == FFieldElt:
            prodPoly = (self.poly * other.poly) % self.field.getIrrPoly()
            newPoly = []
            for i in prodPoly.coef:
                newPoly += [i % 2]
            # pad with 0's if needed
            for i in range(len(self.poly.coef) - len(newPoly)):
                newPoly += [0]
            return FFieldElt(self.field, newPoly)
        if abs(other) == 1:
            return self
        if other == 0:
            return 0

    def __rmul__(self, other):
        assert other == 0 or abs(other) == 1, "You can't multiply " + \
               str(other) + " and " + str(self) + "."
        if other == 0:
            return 0
        if abs(other) == 1:
            return self

    def __pow__(self, exp):
        assert type(exp) == int
        newFieldElt = self
        if exp==0:
            newFieldElt = 1
        for i in range(exp - 1):
            newFieldElt = newFieldElt * self
        return newFieldElt

    def inv(self):
        if self.iszero():
            raise ZeroDivisionError
        powerOfGen = self.field.getPowerOfGenerator(self)
        invPowerOfGen = (-powerOfGen) % (self.field.getSize() - 1)
        return self.field.getEltFromPower(invPowerOfGen)

    def __div__(self, other):
        assert (type(other) == FFieldElt and self.field == other.field)\
               or other == 1 or other == 0
        if type(other) == FFieldElt:
            return self*other.inv()
        if other == 1:
            return self
        if other == 0:
            raise ZeroDivisionError

    def __rdiv__(self, other):
        assert other == 1 or other == 0, "You can't divide " + str(other) + \
               " and " + str(self) + "."
        if self.iszero():
            raise ZeroDivisionError
        if other == 1:
            return self.inv()
        if other == 0:
            return 0

    def __neg__(self):
        return self

    def __pos__(self):
        return self

    def __int__(self):
        if self.iszero():
            return 0
        if self.isone():
            return 1

F = FiniteField((1, 1, 0, 1))
zero = FFieldElt(F, [0, 0, 0])
one = FFieldElt(F, [1, 0, 0])
t = FFieldElt(F, [0, 1, 0])
t2 = t**2
t3 = t**3
t4 = t**4
t5 = t**5
t6 = t**6
