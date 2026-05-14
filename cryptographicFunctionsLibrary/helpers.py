from sage.all import *
from sage.crypto.sbox import SBox
from sage.all import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from collections import defaultdict


def construct_truth_table(F, polynomial):
    """
    Construct the truth table of a polynomial function over GF(2^n) by evaluating it at all field elements.
    """
    return [polynomial(F.from_integer(i)).to_integer() for i in range(F.order())]

def check_apn(F, function):
    """
    Check if a given univariate polynomial or truth table over GF(2^n) is APN by evaluating the differential uniformity.
    """
    if isinstance(function, list):
        tt = function
    elif isinstance(function, Polynomial):
        tt = construct_truth_table(F, function)
    else:
        raise ValueError("Input must be either a univariate polynomial or a truth table")
   
    return SBox(tt).is_apn()

def is_primitive_element(F, e):
    """
    Check if an element e is a primitive element of the finite field F.
    """
    return e != 0 and e.multiplicative_order() == F.order() - 1

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
    x = F(x)
    if x == F(0): 
        return True
    order = F.order() - 1
    if order % 3 != 0:
        return True
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

    term = b_val * c_val**(-(2**i_val))
    if term**q == term:
        s.append(n - i_val)
    
    if is_cube(F, b_val):
        return list(set(s))

    term = b_val**(2**(2 * i_val) - 2**i_val + 1) * c_val**(-1)
    if term**q == term and 0 < 3 * i_val < n:
        s.append(3 * i_val)
    
    term = b_val**(2**i_val - 1) * c_val**(-(2**(2 * i_val)))
    if term**q == term and 0 < (m - 2 * i_val) < n:
        s.append(m - 2 * i_val)

    term = c_val * b_val**(2**i_val - 1)
    if term**q == term and 0 < (m + 2 * i_val) < n:
        s.append(m + 2 * i_val)

    if c_val**q != c_val:
        s.append(m)
        
    if i_val == 1 and gcd(m - 2, n) == 1:
        s_can = pow(m - 2, -1, n)
        term = b_val**(2**(2 * s_can)) * c_val**(-(2**s_can - 1))
        if term**q == term:
            s.append(s_can)
    
    return list(set(s))

def family12_check_s(F, i_val, b_val, c_val, s_val):
    """
    True if s satisfies the conditions for the given (i, b, c).
    """
    n = F.degree()
    m = n // 2
    q = 2**m

    if s_val == n - i_val:
        term = b_val * c_val**(-(2**i_val))
        if term**q == term: 
            return True

    # All remaining cases requires that b is not a cube
    if is_cube(F, b_val):
        return False

    if s_val == 3 * i_val and 0 < s_val < n:
        term = b_val**(2**(2*i_val) - 2**i_val + 1) * c_val**(-1)
        if term**q == term: 
            return True
    
    if s_val == m - 2*i_val and 0 < s_val < n:
        term = b_val**(2**i_val - 1) * c_val**(-(2**(2*i_val)))
        if term**q == term:
            return True
    
    if s_val == m + 2*i_val and 0 < s_val < n:
        term = c_val * b_val**(2**i_val - 1)
        if term**q == term:
            return True
    
    if s_val == m and c_val**q != c_val:
        return True
    
    if i_val == 1 and gcd(m - 2, n) == 1 and s_val == pow(m - 2, -1, n):
        term = b_val**(2**(2*s_val)) * c_val**(-(2**s_val - 1))
        if term**q == term:
            return True
    
    return False

def build_table(pairs):
    """
    Build the result table from a list of (polynomial, parameters) pairs.    
    """
    table = defaultdict(list)
    for poly, params in pairs:
        table[poly].append(params)

    if not table:
        raise TypeError("No valid polynomials found")
    
    if len(table) == 1:
        return next(iter(table))
 
    return dict(table)
