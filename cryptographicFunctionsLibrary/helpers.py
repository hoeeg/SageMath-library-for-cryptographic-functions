from sage.all import *


def is_primitive_element(F, e):
    """
    Check if an element e is a primitive element of the finite field F
    """
    return e.multiplicative_order() == F.order() - 1

def reduce_canonical(F, poly):
    """
    Reduce polynomial to its canonical representative in F_{2^n}[x]/(x^{2^n}+x).
    """
    R = PolynomialRing(F, 'x')
    x = R.gen()
    return R(poly).mod(x**(2**F.degree()) + x)

def get_terms(poly):
    """
    Return {exponent: coefficient} for every nonzero term.
    """
    if poly == 0:
        return {}
    return dict(zip(poly.exponents(), poly.coefficients()))
