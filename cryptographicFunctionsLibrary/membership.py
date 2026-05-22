from sage.all import *
from helpers import is_primitive_element, get_terms, family12_validates_s


def _membership_family1_2(n, p, poly):
    """
    Shared implementation for Family 1 (p=3) and Family 2 (p=4).
    """
    if n < 12:
        return False, {}
    
    if n % p != 0:
        return False, {}
    
    k = n // p
    if gcd(k, 3) != 1:
        return False, {}
    
    F = GF(2**n, 'a')
  
    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    for s in range(1, n):
        if gcd(s, 3 * k) != 1:
            continue
  
        i = (s * k) % p
        m = p - i

        e_x = 2**s + 1                                                  # exponent of x^(2^s + 1)
        e_ux = (2**(i * k) + 2**(m * k + s)) % (2**n - 1)               # exponent of x^(2^(ik) + 2^(mk + s)) mod 2^n - 1

        # Single term case
        if e_x == e_ux:
            if set(terms) != {e_ux}:
                continue
            # Add 1 to the coefficient of x^(2^s + 1) to get u^(2^k - 1)
            u_power = terms.get(e_ux, F(0)) + F(1) 
        
        # Two term case
        else:
            if set(terms) != {e_x, e_ux}:
                continue
            # Extract the coefficients for each term, defaulting to 0 if the term is not present
            if terms.get(e_x, F(0)) != F(1):
                continue
            u_power = terms.get(e_ux)
        
        if u_power == F(0):
            continue

        try:
            all_roots = u_power.nth_root(2**k - 1, all=True)
        except ValueError:
            continue

        primitive_roots = [r for r in all_roots if r != 0 and is_primitive_element(F, r)]
        if not primitive_roots:
            continue

        return True, {'s': s, 'k': k, 'u': primitive_roots}
    
    return False, {}


def membership_family_1(n, poly):
    r"""
    Check if the polynomial is a member of Family 1, the Budaghyan-Carlet-Leander construction from 2008 for `p = 3`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_1
        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33
        sage: membership_family_1(12, poly)
        (True,
        {'s': 5,
        'k': 4,
        'u': [a^11 + a^8 + a^6 + a^3 + a + 1,
        a^11 + a^8 + a^5 + a^3 + a^2 + a,
        a^11 + a^8 + a^6 + a^3 + a^2 + a + 1,
        a^10 + a^8 + a^7 + a^5 + a^4 + a^2 + 1,
        a^6 + a^5 + 1,
        a^11 + a^10 + a^7 + a^6 + a^5 + a^4 + a^3 + a,
        a^10 + a^8 + a^7 + a^6 + a^4,
        a^2,
        a^11 + a^8 + a^5 + a^3 + a,
        a^10 + a^8 + a^7 + a^6 + a^4 + a^2,
        a^11 + a^10 + a^7 + a^4 + a^3 + a^2 + a + 1,
        a^10 + a^8 + a^7 + a^5 + a^4 + 1]})
        
        sage: poly = x^129 + (a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^4)*x^24
        sage: membership_family_1(12, poly)
        (True,
        {'s': 7,
        'k': 4,
        'u': [a^11 + a^6 + a^4 + a^3 + a^2,
        a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1,
        a^11 + a^10 + a^9 + a^8 + a^4 + a^2 + a,
        a^11 + a^10 + a^8 + a^6 + a^5 + a^2,
        a^10 + a^7 + a^6 + a^5 + a^4 + a^3 + a + 1,
        a^9 + a^6 + a^5 + a^4 + a,
        a^11 + a^10 + a^9 + a^7 + a^6 + a^4 + a^2 + 1,
        a^10 + a^9 + a^8 + a^6 + a^3 + a,
        a^11 + a^10 + a^7 + a^5 + a^2 + a + 1,
        a^11 + a^8 + a^7 + a^4 + a^3 + a^2 + a + 1,
        a^8 + a^7 + a^6 + a + 1,
        a^11 + a^9 + a^5 + a^3 + a^2 + a]})
    """
    return _membership_family1_2(n, 3, poly)


def membership_family_2(n, poly):
    r"""
    Check if the polynomial is a member of Family 2, the Budaghyan-Carlet-Leander construction from 2008 for `p = 4`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(i*k) + 2^(mk + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    
        sage: from cryptographicFunctionsLibrary import membership_family_2
        sage: F.<a> = GF(2^16)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3
        sage: membership_family_2(16, poly)
        (True,
        {'s': 1,
        'k': 4,
        'u': [a^12 + a^8 + a^6 + a^5 + a^4 + a^3 + a^2,
        a^15 + a^14 + a^7 + a^6 + a^3 + a,
        a^11 + a^10 + a^9 + a^8 + a^7 + a^3 + a + 1,
        a^13 + a^11 + a^10 + a^7 + a^4 + a + 1,
        a^15 + a^14 + a^12 + a^11 + a^10 + a^9 + a^5 + a^4 + a^3 + a^2 + 1,
        a^15 + a^14 + a^13 +	a^12 +	a^9 +	a^7 +	a^5 +	a^3 +	a^2 +	a,
       	a^15 +	a^14 +	a^13 +	a^12 +	a^11 +	a^10 +	a^8 +
        a^13 + a^12 + a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a + 1]})

        sage: poly = (a^15 + a^14 + a^12 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2)*x^2049
        sage: membership_family_2(16, poly)
        (True,
        {'s': 11,
        'k': 4,
        'u': [a^14 + a^9 + a + 1,
        a^13 + a^12 + a^11 + a^10 + a^4 + a + 1,
        a^14 + a^11 + a^10 + a^9 + a^5 + a^4 + a^3 + a^2 + 1,
        a^11 + a^10 + a^9 + a^4 + a^2 + a,
        a^13 + a^12 + a^5 + a^3 + a^2 + 1,
        a^13 + a^12 + a^11 + a^10 + a^9 + a^5 + a^4 + a^3 + a + 1,
        a^14 + a^13 + a^12 + a^2 + a,
        a^14 + a^11 + a^10 + a^4 + a^2 + 1]})
    """
    return _membership_family1_2(n, 4, poly)


def membership_family_3(n, poly):
    r"""
    Check if the polynomial is a member of Family 3, the Budaghyan-Carlet construction from 2008.
    Defined by `f(x) = sx^(q + 1) + x^(2^i + 1) + x^(q * (2^i + 1)) + cx^(2^i * q + 1) + c^q * x^(2^i + q))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    
        sage: from cryptographicFunctionsLibrary import membership_family_3
        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3
        sage: membership_family_3(6, poly)
        (True, {'q': 8, 'i': 1, 's': a, 'c': a})  

        sage: poly = x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5
        sage: membership_family_3(6, poly)
        (True, {'q': 8, 'i': 2, 's': a^2 + a + 1, 'c': a^5 + a^4 + a^2 + a + 1})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')

    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    for i in range(1, m): 
        if gcd(i, m) != 1:
            continue

        e_s   = q + 1                                               # exponent of sx^(q + 1)
        e_x  = 2**i + 1                                             # exponent of x^(2^i + 1)  
        e_xq  = (q * (2**i + 1)) % (2**n - 1)                       # exponent of x^(q * (2^i + 1))
        e_c   = (2**i * q + 1) % (2**n - 1)                         # exponent of cx^(2^i * q + 1)
        e_cq  = 2**i + q                                            # exponent of x^(2^i + q)

        if not set(terms).issubset({e_s, e_x, e_xq, e_c, e_cq}):
            continue

        if terms.get(e_x) != F(1):
            continue
        if terms.get(e_xq) != F(1):
            continue

        # Extract the coefficients for each term, defaulting to 0 if the term is not present
        s = terms.get(e_s, F(0))
        c = terms.get(e_c, F(0))
        cq = terms.get(e_cq, F(0))

        if c**q != cq:
            continue

        if s not in F or s in GF(q):
            continue

        return True, {'q': q, 'i': i, 's': s, 'c': c}

    return False, {}


def membership_family_4(n, poly):
    r"""
    Check if the polynomial is a member of Family 4, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr_n(a^3 * x^9)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_4
        sage: F.<a> = GF(2^9)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^8 + a^3 + a)*x^288 + (a^8 + a^6 + a^5 + a^4 + a)*x^260 + (a^7 + a^4 + a^3 + a^2 + a)*x^144 + (a^8 + a^6 + a^2 + a)*x^130 + (a^6 + a^5 + a^4)*x^72 + (a^7 + a^6 + a^4 + a^3)*x^65 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^36 + (a^8 + a^4 + a^2 + 1)*x^18 + (a^5 + a^3 + a)*x^9 + x^3
        sage: membership_family_4(9, poly)
        (True, {'a': a^8 + a^6 + a^5 + a^3 + a})

        poly = x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3
        sage: membership_family_4(9, poly)
        (True, {'a': 1})

        sage: F.<a> = GF(2^7)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^6 + a^5 + a^3 + a^2 + a + 1)*x^72 + (a^5 + a^4 + a^2 + 1)*x^68 + (a^5 + a^4 + a^3 + 1)*x^36 + (a^3 + a^2 + 1)*x^34 + (a^6 + a^5 + a^4 + a)*x^18 + (a^5 + a^3 + a^2 + a)*x^17 + (a^6 + a^4 + a^2 + 1)*x^9 + x^3
        sage: membership_family_4(7, poly)
        (True, {'a': a^3 + a^2 + a + 1})
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    # Get the terms of the polynomial reduced modulo x^(2^n) - x
    terms = get_terms(n, poly)
    if not terms:
        return False, {}    

    # All exponents must be 3 or of the form 9*2^i mod (2^n - 1)
    e_9 = {(9 * 2**i) % (2**n - 1) for i in range(n)}
    if not set(terms).issubset({3} | e_9):
        return False, {}

    # Recover a by square root, x^9 = a^(3-1) = a^2
    a_square = terms.get(9 % (2**n - 1), F(0))
    if a_square == F(0):
        return False, {}
    a = a_square**(2**(n - 1))
    
    # Verify the reduced polynomial matches the expected form
    trace = sum(a**(3 * 2**i) * x**((9 * 2**i) % (2**n - 1)) for i in range(n))
    expected = (x**3 + (F(1) / a) * trace).mod(x**(2**n) - x)
    reduced = R(poly).mod(x**(2**n) - x)

    if reduced == expected:
        return True, {'a': a}
    return False, {}


def membership_family_5(n, poly):
    r"""
    Check if the polynomial is a member of Family 5, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^3 * x^9 + a^6 * x^18)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
        sage: from cryptographicFunctionsLibrary import membership_family_5
        sage: F.<a> = GF(2^9)
        sage: R.<x> = PolynomialRing(F)
        sage: sage: poly = (a^8 + a^4 + a)*x^144 + (a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^130 + (a^6 + a^4 + a^3 + a^2 + a + 1)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^65 + (a^7 + a^3 + a^2)*x^18 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3)*x^9 + x^3
        sage: membership_family_5(9, poly)
        (True, {'a': a^7 + a^6 + a^4 + a^3 + 1})

        sage: poly = x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3
        sage: membership_family_5(9, poly)
        (True, {'a': 1})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^11 + a^10 + a^9 + a^8 + a^7 + a^3 + a^2)*x^1152 + (a^10 + a^9 + a^7 + a^3 + a^2)*x^1026 + (a^11 + a^8 + a^6 + a^3 + 1)*x^576 + (a^11 + a^10 + a^6 + a^5 + a^4 + a)*x^513 + (a^10 + a^8 + a^5 + a^3 + a^2 + 1)*x^144 + (a^9 + a^8 + a^5 + a^4 + a^3 + 1)*x^72 + (a^9 + a^8 + a^7 + a^5 + a^3 + 1)*x^18 + (a^11 + a^5 + a^4 + a^2 + 1)*x^9 + x^3
        sage: membership_family_5(12, poly)
        (True, {'a': a^10 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2})
    """
    if n % 3 != 0:
        return False, {}
    
    k = n // 3
    F = GF(2**n, 'a')

    # Get the terms of the polynomial reduced modulo x^(2^n) - x
    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    # x^3 must be present
    if 3 not in terms:
        return False, {}
    
    # All exponents must be 3 or of the form 9*2^3i mod (2^n - 1) or 18*2^3i mod (2^n - 1)
    e_9 = {(9 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    e_18 = {(18 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    if not set(terms).issubset({3} | e_9 | e_18):
        return False, {}
    
    # Recover a by square root, x^9 = a^(3-1) = a^2
    a_square = terms.get(9 % (2**n - 1), F(0))
    if a_square == F(0):
        return False, {}
    a = a_square**(2**(n - 1))
    
    # Verify x^(9 * 2^(3i)) terms: coeff = a^(3 * 2^(3i) - 1) and x^(18 * 2^(3i)) terms: coeff = a^(6 * 2^(3i) - 1)        
    if all(
        terms.get((9 * 2**(3*i)) % (2**n - 1), F(0)) == a ** (3  * 2**(3*i) - 1) and
        terms.get((18 * 2**(3*i)) % (2**n - 1), F(0)) == a ** (6 * 2**(3*i) - 1)
        for i in range(k)
    ):
        return True, {'a': a}

    return False, {}


def membership_family_6(n, poly):
    r"""
    Check if the polynomial is a member of Family 6, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^6 * x^18 + a^12 * x^36)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_6
        sage: F.<a> = GF(2^9)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^7 + a^5 + a^3)*x^288 + (a^8 + a^4 + a^3)*x^260 + (a^6 + a^5 + a^3 + a^2)*x^144 + (a^7 + a^2 + 1)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^36 + (a^8 + a^4 + a^3 + a^2 + 1)*x^18 + x^3
        sage: membership_family_6(9,  poly)
        (True, {'a': a^8 + a^7 + a^2})

        sage: poly = x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3
        sage: membership_family_6(9,  poly)
        (True, {'a': 1})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^2304 + (a^11 + a^7 + a^6 + a^5 + a^4 + a^3)*x^2052 + (a^10 + a^7 + a^6 + a^4 + a + 1)*x^1152 + (a^10 + a^8 + a^7 + a^6 + a^3 + 1)*x^1026 + (a^11 + a^8 + a^5 + a^2 + a)*x^288 + (a^10 + a^8 + a^6 + a^5 + a^4 + a^3 + a^2)*x^144 + (a^9 + a^8 + a^5 + a^4 + a^3 + a^2 + a)*x^36 + (a^9 + a^7 + a^6 + a^5 + a^3 + a^2)*x^18 + x^3
        sage: membership_family_6(12, poly)
        (True, {'a': a^11 + a^10 + a^9 + a^8 + a^7 + a^5 + a^4 + a})
    """
    if n % 3 != 0:
        return False, {}
    
    k = n // 3
    F = GF(2**n, 'a')

    # Get the terms of the polynomial reduced modulo x^(2^n) - x
    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    # x^3 must be present with coefficient 1
    if terms.get(3, F(0)) != F(1):
        return False, {}
    
    # All exponents must be 3 or of the form 18*2^3i mod (2^n - 1) or 36*2^3i mod (2^n - 1)
    e_18 = {(18 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    e_36 = {(36 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    if not set(terms).issubset({3} | e_18 | e_36):
        return False, {}
    
    # Recover a by fifth root, x^18 = a^(6-1) = a^5
    a_fifth = terms.get(18 % (2**n - 1), F(0))
    if a_fifth == F(0):
        return False, {}
    
    try:
        all_roots = a_fifth.nth_root(5, all=True)
    except ValueError:
        return False, {}
    
    # Verify all fifth roots to find a valid a such that all x^(18 * 2^(3i)) terms: coeff = a^(6 * 2^(3i) - 1) and x^(36 * 2^(3i)) terms: coeff = a^(12 * 2**(3i) - 1)
    for a in all_roots:
        if a == F(0):
            continue
        
        if all(
            terms.get((18 * 2**(3*i)) % (2**n - 1), F(0)) == a ** (6  * 2**(3*i) - 1) and
            terms.get((36 * 2**(3*i)) % (2**n - 1), F(0)) == a ** (12 * 2**(3*i) - 1)
            for i in range(k)
        ):
            return True, {'a': a}

    return False, {}


def membership_family_7_9(n, poly):
    r"""
    Check if the polynomial is a member of Family 7 - Family 9, the Bracken-Byrne-Markin-McGuire construction from 2011.
    Defined by `f(x) = ux^(2^s + 1) + u^(2^k) * x^(2^-k + 2^(k + s)) + vx^(2^-k + 1) + wu^(2^k + 1) * x^(2^s + 2^(k + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_7_9
        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^10 + a^6 + a^5 + a^3 + a)*x^768 + (a^9 + a^8 + a^6 + a^5 + a^2)*x^544 + (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^257 + (a^7 + a^5 + a^3 + a^2 + 1)*x^33
        sage: membership_family_7_9(12, poly)
        (True, {'s': 5, 'u': a^7 + a^5 + a^3 + a^2 + 1, 'v': a4^3 + a4, 'w': a4^3 + a4})

        sage: poly = a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + x^257
        sage: membership_family_7_9(12, poly)
        (True, {'s': 11, 'u': a, 'v': 1, 'w': 0})

        sage: poly =  (a^10 + a^3 + a^2 + 1)*x^2056 + a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257
        sage: membership_family_7_9(12, poly)
        (True, {'s': 11, 'u': a, 'v': a4^2, 'w': a4^2})
    """ 
    if n % 3 != 0:
        return False, {}
    
    k = n // 3
    if gcd(k, 3) != 1:
        return False, {}
    
    F = GF(2**n, 'a')
    K = F.subfield(k)
  
    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    for s in range(1, n):
        if gcd(s, 3*k) != 1 or (k + s) % 3 != 0:
            continue

        e1 = (2**s + 1)
        e2 = (2**(n - k) + 2**(k + s)) % (2**n - 1)
        e3 = (2**(n - k) + 1) % (2**n - 1)
        e4 = (2**s + 2**(k + s)) % (2**n - 1)

        if not set(terms).issubset({e1, e2, e3, e4}):
            continue

        u = terms.get(e1,  F(0))
        uk = terms.get(e2, F(0))
        v = terms.get(e3,  F(0))
        wu = terms.get(e4,  F(0))

        if u == F(0) or not is_primitive_element(F, u):
            continue
 
        
        w = wu / (u**(2**k + 1))
        if v * w == K(1) or v not in K or w not in K:
            continue

        return True, {'s': s, 'u': u, 'v': K(v), 'w': K(w)}

    return False, {}


def membership_family_11(n, poly):
    r"""
    Check if the polynomial is a member of Family 11, the Budaghyan-Helleseth-Kaleyski construction from 2020.
    Defined by `f(x) = x^3 + a(x^2^i + 1) + bx^(3 * 2^m) + c(x^(2^(i + m) + 2^m))^2^k`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_11
        sage: F.<a> = GF(2^10)
        sage: R.<x> = PolynomialRing(F)
        sage: poly =  (a^5 + a^3 + a)*x^576 + (a^5 + a^3 + a + 1)*x^96 + x^18 + x^3
        sage: membership_family_11(10, poly)
        (True, {'k': 6, 'i': 3, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: poly = x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3
        sage: membership_family_11(10, poly)
        (True, {'k': 2, 'i': 9, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: F.<a> = GF(2^14)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^768 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a + 1)*x^384 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^6 + x^3
        sage: membership_family_11(14, poly)
        (True,
        {'k': 2,
        'i': 13,
        'a': a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a,
        'b': a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a + 1,
        'c': 1})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    if m % 2 == 0 or m % 3 == 0:
        return False, {}
    
    F = GF(2**n, 'a')
    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    if terms.get(3, F(0)) != F(1):
        print("test 1 failed")
        return False, {}
    
    def _make_i_set(ref, base):
        inv_set = {pow(ref, -1, n)} if gcd(ref, n) == 1 else set()
        return base | inv_set

    i_list_even = _make_i_set(m - 2, {m - 2, m, n - 1})
    i_list_odd  = _make_i_set(m + 2, {m + 2, m})
    
    for k in range(0, n):
        i_list = i_list_even if k % 2 == 0 else i_list_odd
        for i in i_list:
            e_a = ((2**i + 1) * 2**k) % (2**n - 1)                  # exponent of a(x^2^i + 1)^2^k
            e_b = (3 * 2**m)                                        # exponent of bx^(3 * 2^m)
            e_c = ((2**(i + m) + 2**m) * 2**k) % (2**n - 1)         # exponent of c(x^(2^(i + m) + 2^m))^2^k
    
            if not set(terms).issubset({3, e_a, e_b, e_c}):
                continue
            
            # Single term case
            if e_a == e_c:
                a = terms.get(e_a, F(0)) + F(1)
                c = F(1)
            # Two term case
            else:
                a = terms.get(e_a, F(0))
                c = terms.get(e_c, F(0))
            
            b = terms.get(e_b, F(0))
            if b != a**2 or c != F(1):
                continue
       
            # a must be a primitive element of GF(2^2)
            if a == F(0) or not is_primitive_element(F.subfield(2), a):
                continue
            
            return True, {'k': k, 'i': i, 'a': a, 'b': b, 'c': c}
        
    return False, {}


def membership_family_12(n, poly):
    r"""
    Check if the polynomial is a member of Family 12, the Zheng-Kan-Li-Peng-Tang from 2022.
    Defined by `f(x) = a * Tr^n_m(bx^(2^i + 1)) + a^q * Tr^n_m(cx^(2^s + 1))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    
    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_12
        sage: F.<a> = GF(2^10)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^8 + a^6 + a^4 + a^3 + 1)*x^96 + (a^9 + a^7 + a^6 + a^4)*x^33 + (a^8 + a^6 + a^5 + a^3)*x^3
        sage: membership_family_12(10, poly)
        (True,
        {'i': 1,
        's': 5,
        'a': a^8 + a^5 + a^4 + a^3 + a^2,
        'b': a^8 + a^6 + a^5 + a^4 + a^3 + a + 1,
        'c': a})

        sage: poly = a^3*x^513 + (a^8 + a^7 + a^2 + a + 1)*x^48 + (a^3 + a^2 + 1)*x^33
        sage: membership_family_12(10, poly)
        (True,
        {'i': 9,
        's': 5,
        'a': a^9 + a^6 + a^2 + a + 1,
        'b': a^6 + a^5 + a^3,
        'c': a^9 + a^7 + a^6 + a^4 + a^3 + a})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    q = 2**m
    if m % 2 == 0:
        return False, {}
    
    F = GF(2**n, 'a')
    Q = GF(q)

    # Get the terms of the polynomial reduced modulo x^(2^n) - x
    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    for i in range(1, n):
        if gcd(i, n) != 1:
            continue

        for s in range(1, n):
            if i == s:
                continue

            e_x1 = (2**i + 1)
            e_b = (q * e_x1) % (2**n - 1)
            e_x2 = (2**s + 1)
            e_c = (q * e_x2) % (2**n - 1)

            if not set(terms).issubset({e_x1, e_b, e_x2, e_c}):
                continue

            c1 = terms.get(e_x1, F(0))
            c2 = terms.get(e_b, F(0))
            c3 = terms.get(e_x2, F(0))

            if c1 == F(0) or c2 == F(0) or c3 == F(0):
                continue

            # Recover a^(q-1) = c1^q / c2
            try:
                all_roots = (c1**q / c2).nth_root(q - 1, all=True)
            except ValueError:
                continue

            for a in all_roots:
                if a == F(0) or a in Q:
                    continue
            
                b = c1 / a
                if b == F(0) or a*b**q != c2:
                    continue

                c_can = c3 / a**q
                if e_x2 == e_c:
                    c_list = [c for c in F if c != F(0) and c**q + c == c_can]
                elif c_can == F(0):
                    continue
                else:
                    c_list = [c_can]

                for c in c_list:
                    if family12_validates_s(F, i, b, c, s):
                        return True, {'i': i, 's': s, 'a': a, 'b': b, 'c': c}
    
    return False, {}


def membership_family_13(n, poly):
    r"""
    Check if the polynomial is a member of Family 13, the Li-Zhou-Li-Qu construction from 2022.
    Defined by `f(x) = L(z)^(2^m + 1) + cz^(2^m + 1)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_family_13
        sage: F.<a> = GF(2^9)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^144 + (a^8 + a^7 + a^3 + a^2 + a)*x^130 + x^129 + (a^8 + a^2)*x^32 + x^24 + (a^7 + a^6 + a^4 + a^3 + a)*x^18 + (a^8 + a^2)*x^17 + (a^8 + a^7 + a^3 + a^2 + a)*x^10 + (a^4 + a^3 + a^2 + 1)*x^9
        sage: membership_family_13(9, poly)
        (True, {'s': 1, 'v': a^4 + a^3 + a^2, 'mu': a^8 + a^7 + a^3 + a^2 + a})

        sage: poly = x^288 + (a^8 + a^7 + a^6 + a^3 + a^2)*x^260 + x^257 + (a^8 + a^7 + a^5 + a^3 + a^2 + a + 1)*x^64 + x^40 + (a^7 + a^3 + a^2)*x^36 + (a^8 + a^7 + a^5 + a^3 + a^2 + a + 1)*x^33 + (a^8 + a^7 + a^6 + a^3 + a^2)*x^12
        sage: membership_family_13(9, poly)
        (True, {'s': 2, 'v': 1, 'mu': a^8 + a^7 + a^6 + a^3 + a^2})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^544 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^514 + x^513 + (a^7 + a^6 + a^4 + a^3 + a)*x^64 + x^48 + (a^11 + a^9 + a^5)*x^34 + (a^7 + a^6 + a^4 + a^3 + a)*x^33 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^18 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^17
        sage: membership_family_13(12, poly)
        (True,
        {'s': 1,
        'v': a^11 + a^9 + a^8 + a^6 + a^3 + a,
        'mu': a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1})
    """
    if n % 3 != 0:
        return False, {}
    
    m = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    Fm = F.subfield(m)

    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    for s_val in range(1, m + 1):
        if gcd(s_val, m) != 1:
            continue

        P = 2**(2*m + s_val)
        Q = 2**(m + s_val)
        Rp = 2**m
        B = 2**s_val

        e_PQ = (P + Q) % (2**n - 1)
        e_PB = (P + B) % (2**n - 1)
        e_P1 = (P + 1) % (2**n - 1)
        e_2Q = (2 * Q) % (2**n - 1)
        e_QB = (Q + B) % (2**n - 1)
        e_Q1 = (Q + 1) % (2**n - 1)
        e_RQ = (Rp + Q) % (2**n - 1)
        e_RB = (Rp + B) % (2**n - 1)
        e_R1 = (Rp + 1) % (2**n - 1)

        expected = {e_PQ, e_PB, e_P1, e_2Q, e_QB, e_Q1, e_RQ, e_RB, e_R1}
        if not set(terms).issubset(expected):
            continue

        # Recover mu from x^(P+B)
        mu = terms.get(e_PB, F(0))
        if mu == F(0):
            continue

        if terms.get(e_PQ, F(0)) != F(1) or terms.get(e_P1, F(0)) != F(1) or terms.get(e_RQ, F(0)) != F(1):
            continue
        if terms.get(e_RB, F(0)) != mu:
            continue
        mu_q = mu**(2**m)
        if terms.get(e_2Q, F(0)) != mu_q or terms.get(e_Q1, F(0)) != mu_q:
            continue
        if terms.get(e_QB, F(0)) != mu_q * mu:
            continue

        # v from x^(R+1) coefficient = 1 + v
        v = terms.get(e_R1, F(0)) + F(1)
        if v == F(0) or v not in Fm:
            continue

        # mu^(q^2 + q + 1) != 1
        if mu**(2**(2*m) + 2**m + 1) == F(1):
            continue

        # L permutes F
        L = x**(2**(m + s_val)) + mu*x**(2**s_val) + x
        if L.gcd(x**(2**n) - x) != x:
            continue

        return True, {'s': s_val, 'v': v, 'mu': mu}

    return False, {}


FAMILIES = {
    "Family 1": membership_family_1,
    "Family 2": membership_family_2,
    "Family 3": membership_family_3,
    "Family 4": membership_family_4,
    "Family 5": membership_family_5,
    "Family 6": membership_family_6,
    "Family 7-9": membership_family_7_9,
    "Family 11": membership_family_11,
    "Family 12": membership_family_12,
    "Family 13": membership_family_13
}

def membership_all(n, polynomial):
    r"""
    Check if a given function is a member of any of the known infinite families of quadratic APN polynomials. 
    
    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import membership_all
        sage: F.<a> = GF(2^16)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^14 + a^12 + a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^2049
        sage: membership_all(16, polynomial)
        Belongs to Family 2, with parameters: 
            {'s': 11, 'k': 4, 'u': 
                [a^14 + a^13 + a^12 + a^11 + a^7 + a^6 + a^3 + a^2 + a, 
                ...
                a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^4]}
        
        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3
        sage: membership_all(6, polynomial)
        Belongs to Family 3, with parameters: {'q': 8, 'i': 1, 's': a, 'c': a}
    """
    for name, function in FAMILIES.items():
        try:
            found, params = function(n, polynomial)
        except Exception:
            continue
        if found:
            print(f"Belongs to {name}, with parameters: {params}")
