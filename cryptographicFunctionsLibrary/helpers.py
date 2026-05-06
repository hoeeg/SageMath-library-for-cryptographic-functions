from sage.all import *
from sage.crypto.sbox import SBox
from sage.all import Polynomial


def construct_truth_table(F, polynomial):
    """
    Construct the truth table of a polynomial function over GF(2^n) by evaluating it at all field elements.
    """
    return [polynomial(F.from_integer(i)).to_integer() for i in range(F.order())]

def is_apn(F, function):
    """
    Check if a given univariate polynomial or truth table over GF(2^n) is APN by evaluating the differential uniformity.
    """
    if isinstance(function, list):
        tt = function
    elif isinstance(function, Polynomial):
        tt = construct_truth_table(F, function)
    else:
        raise ValueError("Input must be either a univariate polynomial or a truth table")
   
    return SBox(tt).differential_uniformity() == 2

def is_primitive_element(F, e):
    """
    Check if an element e is a primitive element of the finite field F.
    """
    return e.multiplicative_order() == F.order() - 1

def get_terms(n, poly):
    """
    Return {exponent: coefficient} for every nonzero term in the polynomial reduced modulo x^(2^n) - x.
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    poly_reduced = R(poly).mod(x**(2**n) - x)
    return dict(zip(poly_reduced.exponents(), poly_reduced.coefficients())) if poly_reduced != 0 else {}

def is_cube(F, x):
    """
    Check if an element is a cube in GF(2^n).
    """
    if x == F(0): 
        return True
    order = F.order() - 1
    return x**(order // 3) == F(1)

def family12_conditions(F, i_val, b_val, c_val):
    """
    Conditions for s in family 12, given F, i, b and c. 
    Returns a list of s values that satisfy the conditions.
    """
    n = F.degree()
    m = n // 2
    q = 2**m
    s = []

    if not is_cube(F, b_val):
        term = b_val**(2**(2 * i_val) - 2**i_val + 1) * c_val**(-1)
        if term**q == term:
            s.append(3 * i_val)
        
        term = b_val**(2**i_val - 1) * c_val**(-(2**(2 * i_val)))
        if term**q == term and (m - 2 * i_val) > 0:
            s.append(m - 2 * i_val)

        term = c_val * b_val**(2**i_val - 1)
        if term**q == term:
            s.append(m + 2 * i_val)

        if i_val == 1:
            for s_can in (1, n):
                if (s_can * (m - 2)) % n == 1:
                    term = b_val**(2**(2 * s_can)) * c_val**(-(2**s_can - 1))
                    if term**q == term:
                        s.append(s_can)
        
        if c_val**q != c_val:
            s.append(m)
        
    term = b_val * c_val**(-(2**i_val))
    if term**q == term:
        s.append(n - i_val)

    return list(set(s))
