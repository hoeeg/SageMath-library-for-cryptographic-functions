from math import gcd

from sage.all import *
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import Matrix
from sage.crypto.sbox import SBox
from sage.all import Polynomial, SR
from sage.crypto.sboxes import monomial_function

# from sage.all import GF, PolynomialRing, Matrix, Polynomial, SBox, SR

def is_apn(F, function):
    """
    Check if a given univariate polynomial over GF(2^n) is APN by evaluating the differential uniformity.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``function`` -- either
        - a univariate polynomial over GF(2^n)
        - a truth table (look up table) represented as a list of integers
    """
    if isinstance(function, list):
        tt = function
    elif isinstance(function, Polynomial):
        tt = construct_truth_table(F, function)
    else:
        raise ValueError("Input must be either a univariate polynomial or a truth table")
   
    return SBox(tt).differential_uniformity() == 2


def polynomial_to_matrix(F, polynomial, basis, output_format='univariate'):
    """
    Compute the quadratic APN matrix (QAM) of a univariate polynomial over GF(2^n) with respect to a normal basis.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    - ``basis`` -- a normal basis of GF(2^n) over GF(2)
    - ``output_format`` -- default as 'univariate', if 'power', output matrix entries as powers of the field generator a
    
    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import polynomial_to_matrix
        sage: F.<a> = GF(2^4)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: basis = [(a^3)**(2^i) for i in range(4)]
        sage: M = polynomial_to_matrix(F, polynomial, basis, output_format='power'); M
        [   0 a^11 a^10 a^13]
        [a^11    0  a^7  a^5]
        [a^10  a^7    0 a^14]
        [a^13  a^5 a^14    0]

        sage: from cryptographicFunctionsLibrary import polynomial_to_matrix
        sage: F.<a> = GF(2^4)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: basis = [(a^3)**(2^i) for i in range(4)]
        sage: M = polynomial_to_matrix(F, polynomial, basis); M
        [            0 a^3 + a^2 + a   a^2 + a + 1 a^3 + a^2 + 1]
        [a^3 + a^2 + a             0   a^3 + a + 1       a^2 + a]
        [  a^2 + a + 1   a^3 + a + 1             0       a^3 + 1]
        [a^3 + a^2 + 1       a^2 + a       a^3 + 1             0]
    """
    if not is_apn(F, polynomial):
        raise ValueError("The provided polynomial is not a quadratic APN function")

    if output_format == 'power':
        a = F.gen()
        a_sym = SR.var('a')

        def to_power(val):
            if val == 0: 
                return 0
            elif val == 1:
                return 1
            return a_sym**val.log(a)
        
        return Matrix(SR, [[to_power(polynomial(bi + bj) + polynomial(bi) + polynomial(bj) + polynomial(0)) for bj in basis] for bi in basis])


    return Matrix(F, [[polynomial(bi + bj) + polynomial(bi) + polynomial(bj) + polynomial(0) for bj in basis] for bi in basis])



def matrix_to_polynomial(F, M, basis):
    """
    Compute the univariate polynomial representation from a matrix over GF(2^n) with respect to a normal basis.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``M`` -- a matrix over GF(2^n)
    - ``basis`` -- a normal basis of GF(2^n) over GF(2)
    
    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import matrix_to_polynomial
        sage: F.<a> = GF(2^4)
        sage: basis = [(a^3)**(2^i) for i in range(4)]
        sage: M = Matrix(F, [[ 0, a^3 + a^2 + a, a^2 + a + 1, a^3 + a^2 + 1], [a^3 + a^2 + a, 0, a^3 + a + 1, a^2 + a], [ a^2 + a + 1, a^3 + a + 1, 0, a^3 + 1], [a^3 + a^2 + 1, a^2 + a, a^3 + 1, 0]])
        sage: polynomial = matrix_to_polynomial(F, M, basis); polynomial
        x^3
    """
    if M.base_ring() != F:
        raise ValueError("Matrix base ring must match the field F")
        
    n = F.degree()
    R = PolynomialRing(F, 'x')
    x = R.gen()

    M_alpha = Matrix(F, [[b**(2**j) for b in basis] for j in range(n)])
    CF = M_alpha.inverse().transpose() * M * M_alpha.inverse()
    
    polynomial = R.zero()
    for i in range(n):
        for j in range(i):
            if CF[i,j] != 0:
                polynomial += CF[i,j] * x**(2**i + 2**j)

    return polynomial


# DO NOT KNOW IF NEEDED
# ONLY IF WE WANT CUSTOM BASIS
def field_element_to_binary_integer(a, basis):
    """
    Encode a field element of GF(2^n) as an integer via its coordinate vector.
    """
    vec = a.vector(basis=basis)
    return sum(int(vec[i]) << i for i in range(len(vec)))


# DO NOT KNOW IF NEEDED
# ONLY IF WE WANT CUSTOM BASIS
def binary_integer_to_field_element(F, basis, integer):
    """
    Decode a binary integer into a field element of GF(2^n) using the given GF(2)-basis.
    """
    bits = Integer(integer).digits(base=2, padto=F.degree())
    return sum(bit * b for bit, b in zip(bits, basis))



def apn_function_to_algebra_sequence(F, polynomial):
    """
    Compute the APN algebra sequence associated to a quadratic APN function F.

    INPUT:
    - ``F`` -- a finite field GF(2^n) 
    - ``polynomial`` -- a polynomial over F
    
    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import apn_function_to_algebra_sequence
        sage: F.<a> = GF(2^7)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: seq = apn_function_to_algebra_sequence(F, polynomial)
        [6, 20, 72, 22, 56, 32, 48, 35, 76, 51, 69, 5, 30, 108, 29, 40, 115, 106, 70, 17, 60]
    """
    n = F.degree()
    basis = [F.gen()**i for i in range(n)]

    if not is_apn(F, polynomial):
        raise ValueError("The provided polynomial is not a quadratic APN function")

    M = polynomial_to_matrix(F, polynomial, basis)

    sequence = [M[i, j].to_integer() for i in range(n) for j in range(i + 1, n)]
    print(sequence)

    return sequence


def algebra_sequence_to_apn_function(F, sequence):
    """
    Reconstruct a quadratic APN polynomial from its APN algebraic sequence.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``sequence`` -- an APN algebraic sequence
    
    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import algebra_sequence_to_apn_function
        sage: F.<a> = GF(2^7)
        sage: seq = [6, 20, 72, 22, 56, 32, 48, 35, 76, 51, 69, 5, 30, 108, 29, 40, 115, 106, 70, 17, 60]
        sage: polynomial = algebra_sequence_to_apn_function(F, seq); polynomial
        x^3
    """
    n = F.degree()
    basis = [F.gen()**i for i in range(n)]

    seq_iter = (F.from_integer(s) for s in sequence)

    M = Matrix(F, n)
    for i in range(n):
        for j in range(i + 1, n):
            M[i,j] =  M[j,i] = next(seq_iter)

    return matrix_to_polynomial(F, M, basis)


def construct_truth_table(F, polynomial):
    """
    Construct the truth table of a polynomial function over GF(2^n) by evaluating it at all field elements.
    """
    return [polynomial(F.from_integer(i)).to_integer() for i in range(F.order())]


def polynomial_to_truth_table(F, polynomial):
    """
    Evaluates a polynomial for all elements in the field F to create a truth table.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import polynomial_to_truth_table
        sage: F.<a> = GF(2^3)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: tt = polynomial_to_truth_table(F, polynomial); tt
        [0, 1, 3, 4, 5, 6, 7, 2]
    """
    tt = construct_truth_table(F, polynomial)

    if not is_apn(F, tt):
        raise ValueError("The provided polynomial is not a quadratic APN function")
    
    return tt


def truth_table_to_polynomial(F, tt: list):
    """
    Convert a truth table (look up table) into a univariate polynomial representation using Lagrange interpolation.

    INPUT:
    - ``F`` -- a finite field GF(2^n)
    - ``tt`` -- a truth table (look up table) represented as a list of integers

    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import truth_table_to_polynomial
        sage: F.<a> = GF(2^3)
        sage: tt = [0, 1, 3, 4, 5, 6, 7, 2]
        sage: polynomial = truth_table_to_polynomial(F, tt); polynomial
        x^3

        sage: F.<a> = GF(2^8)
        sage: tt = [0, 1, 0, 6, 0, 12, 2, 9, 0, 24, 50, 45, 4, 17, 52, 38, 0, 48, 47, 24, 100, 89, 73, 115, 8, 33, 21, 59, 104, 76, 119, 84, 0, 96, 229, 130, 93, 48, 186, 208, 200,
        ....: 177, 31, 97, 145, 229, 68, 55, 16, 65, 218, 140, 41, 117, 225, 186, 208, 152, 40, 103, 237, 168, 23, 85, 0, 192, 76, 139, 203, 6, 133, 79, 185, 96, 199, 25, 118, 162,
        ....: 10, 217, 146, 99, 241, 7, 61, 193, 92, 167, 35, 203, 114, 157, 136, 109, 219, 57, 32, 129, 137, 47, 182, 26, 29, 182, 81, 233, 202, 117, 195, 118, 90, 232, 162, 50, 36
        ....: , 179, 80, 205, 212, 78, 219, 82, 111, 225, 45, 169, 155, 24, 0, 130, 128, 5, 152, 23, 26, 146, 151, 12, 37, 185, 11, 157, 187, 42, 115, 192, 220, 104, 143, 49, 34, 15
        ....: 5, 236, 70, 113, 220, 20, 179, 139, 43, 38, 197, 67, 167, 227, 13, 132, 109, 121, 131, 46, 211, 184, 79, 237, 29, 69, 151, 15, 218, 228, 59, 172, 116, 18, 217, 106, 16
        ....: 6, 183, 113, 205, 12, 64, 3, 140, 200, 19, 93, 221, 148, 110, 52, 144, 205, 57, 110, 197, 149, 161, 211, 66, 55, 150, 233, 119, 15, 135, 236, 86, 58, 180, 210, 103, 6,
        ....:  70, 100, 111, 74, 72, 103, 99, 75, 160, 155, 187, 135, 170, 156, 179, 130, 183, 164, 177, 165, 221, 195, 217, 192, 89, 83, 109, 96, 55, 48, 1, 1]
        sage: polynomial = truth_table_to_polynomial(F, tt); polynomial
        (a^7 + a^6 + a^4 + a^3 + a + 1)*x^192 + (a^7 + a^6 + a^2 + a + 1)*x^160 + (a^7 + a^6 + a^5)*x^144 + (a^6 + a^5 + a^4 + a^3 + a^2)*x^136 + (a^6 + a^5 + a^3 + a)*x^132 + 
        (a^7 + a^5 + a^3 + a)*x^130 + (a^5 + a^4 + a + 1)*x^129 + (a^7 + a^5 + a^4 + a^3 + a)*x^128 + (a^6 + a^5 + a^4 + a^3 + 1)*x^96 + (a^5 + a^4 + a^2 + 1)*x^80 + 
        (a^6 + a^4 + a^3 + a^2 + a)*x^72 + (a^5 + a^4 + a)*x^68 + (a^7 + a^4 + a + 1)*x^66 + (a^7 + a^5 + 1)*x^65 + (a^4 + a^3 + a^2 + a)*x^64 + (a^3 + a^2)*x^48 + (a^4 + a^3 + a^2 + a + 1)*x^40 + 
        (a^7 + a^4 + a + 1)*x^36 + (a^5 + a^4 + a^3 + a^2)*x^34 + (a^7 + a^3 + a^2)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a + 1)*x^32 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^24 + (a^6 + a^5 + 1)*x^20 + 
        (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^18 + (a^6 + a^5 + a^4)*x^17 + (a^7 + 1)*x^16 + (a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^12 + (a^7 + a^5 + a^4 + 1)*x^10 + (a^7 + a^6 + a^5 + a^4 + a + 1)*x^9 + 
        (a^7 + a^3)*x^8 + (a^6 + a^3 + a^2 + 1)*x^6 + (a^7 + a^6)*x^5 + (a^6 + a^5 + a^2 + a)*x^4 + (a^7 + a^6 + a^5 + a^3 + a)*x^3 + (a^7 + a^6 + a^4 + a + 1)*x^2 + (a^6 + a^5 + a^3 + a^2 + a + 1)*x
    """
    R = PolynomialRing(F, 'x')
    
    points = [(F.from_integer(i), F.from_integer(val)) for i, val in enumerate(tt)]

    return R.lagrange_polynomial(points)
