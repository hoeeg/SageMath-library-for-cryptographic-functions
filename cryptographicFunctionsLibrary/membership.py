from sage.all import *
from helpers import is_primitive_element


def get_terms(poly):
    """
    Return {exponent: coefficient} for every nonzero term.
    """
    if poly == 0:
        return {}
    return dict(zip(poly.exponents(), poly.coefficients()))


def _belong_family1_2(n, p, poly):
    """
    Shared implementation for family1 (p=3) and family2 (p=4).
    """
    if n < 12:
        raise ValueError("n must be at least 12")
    
    if n % p != 0:
        return False, {}
    
    k = n // p
    if gcd(k, 3) != 1:
        return False, {}
    
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()    
    terms = get_terms((poly.mod(x**(2**n) - x)))

    for s in range(1, n):
        if gcd(s, 3 * k) != 1:
            continue
  
        i = (s * k) % p
        m = p - i

        e_x = 2**s + 1                                                  # exponent of x^(2^s + 1)
        e_ux = (2**(i * k) + 2**(m * k + s)) % (2**n - 1)               # exponent of x^(2^(ik) + 2^(mk + s)) mod 2^n - 1

        # Single term case
        if e_x == e_ux:
            if set(terms) != {e_x}:
                continue
            # Add 1 to the coefficient of x^(2^s + 1) to get u^(2^k - 1)
            u_power = terms.get(e_x, F(0)) + F(1) 
        
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

        primitive_roots = [r for r in all_roots if is_primitive_element(F, r)]
        if not primitive_roots:
            continue

        return True, {'s': s, 'k': k, 'u': primitive_roots}
    
    return False, {}


def belong_family1(n, poly):
    r"""
    Check if the polynomial belongs to Family1, the Budaghyan-Carlet-Leander construction from 2008 for `p = 3`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import belong_family1
        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33
        sage: belong_family1(12, poly)
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
        sage: belong_family1(12, poly)
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
    return _belong_family1_2(n, 3, poly)


def belong_family2(n, poly):
    r"""
    Check if the polynomial belongs to Family2, the Budaghyan-Carlet-Leander construction from 2008 for `p = 4`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(i*k) + 2^(mk + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    
        sage: from cryptographicFunctionsLibrary import belong_family2
        sage: F.<a> = GF(2^16)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3
        sage: belong_family2(16, poly)
        (True,
        {'s': 1,
        'k': 4,
        'u': [a^12 + a^8 + a^6 + a^5 + a^4 + a^3 + a^2,
        a^15 + a^14 + a^7 + a^6 + a^3 + a,
        a^11 + a^10 + a^9 + a^8 + a^7 + a^3 + a + 1,
        a^13 + a^11 + a^10 + a^7 + a^4 + a + 1,
        a^15 + a^14 + a^12 + a^11 + a^10 + a^9 + a^5 + a^4 + a^3 + a^2 + 1,
        a^15 + a^14 + a^13 + a^12 + a^9 + a^7 + a^5 + a^3 + a^2 + a,
        a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^8 + a^5 + a^2 + 1,
        a^13 + a^12 + a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a + 1]})

        sage: poly = (a^15 + a^14 + a^12 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2)*x^2049
        sage: belong_family2(16, poly)
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
    return _belong_family1_2(n, 4, poly)


def belong_family3(n, poly):
    r"""
    Check if the polynomial belongs to Family3, the Budaghyan-Carlet construction from 2008.
    Defined by `f(x) = sx^(q + 1) + x^(2^i + 1) + x^(q * (2^i + 1)) + cx^(2^i * q + 1) + c^q * x^(2^i + q))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    
        sage: from cryptographicFunctionsLibrary import belong_family3
        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3
        sage: belong_family3(6, poly)
        (True, {'q': 8, 'i': 1, 's': a, 'c': a})  

        sage: poly = x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5
        sage: belong_family3(6, poly)
        (True, {'q': 8, 'i': 2, 's': a^2 + a + 1, 'c': a^5 + a^4 + a^2 + a + 1})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()    
    terms = get_terms((poly.mod(x**(2**n) - x)))

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

def belong_family4(n, poly):
    r"""
    Check if the polynomial belongs to Family4, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr_n(a^3 * x^9)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    """
    # Placeholder for future implementation
    return False, {}


def belong_family5(n, poly):
    r"""
    Check if the polynomial belongs to Family5, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^3 * x^9 + a^6 * x^18)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    """
    # Placeholder for future implementation
    return False, {}


def belong_family6(n, poly):
    r"""
    Check if the polynomial belongs to Family6, the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^6 * x^18 + a^12 * x^36)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    """
    # Placeholder for future implementation
    return False, {}


def belong_family7_9(n, poly):
    r"""
    Check if the polynomial belongs to Family7 - Family9, the Bracken-Byrne-Markin-McGuire construction from 2011.
    Defined by `f(x) = ux^(2^s + 1) + u^(2^k) * x^(2^-k + 2^(k + s)) + vx^(2^-k + 1) + wu^(2^k + 1) * x^(2^s + 2^(k + s))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import belong_family7_9
        sage: F.<a> = GF(2^3)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = (a^10 + a^6 + a^5 + a^3 + a)*x^768 + (a^9 + a^8 + a^6 + a^5 + a^2)*x^544 + (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^257 + (a^7 + a^5 + a^3 + a^2 + 1)*x^33
        sage: belong_family7_9(12, poly)
        (True,
        {'s': 5,
        'k': 4,
        'u': a^7 + a^5 + a^3 + a^2 + 1,
        'v': a^10 + a^9 + a^6 + a^5 + a^3 + 1,
        'w': a^10 + a^9 + a^6 + a^5 + a^3 + 1})

        sage: poly = a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + x^257
        sage: belong_family7_9(12, poly)
        (True, {'s': 11, 'k': 4, 'u': a, 'v': 1, 'w': 0})

        sage: poly =  (a^10 + a^3 + a^2 + 1)*x^2056 + a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257
        sage: belong_family7_9(12, poly)
        (True,
        {'s': 11,
        'k': 4,
        'u': a,
        'v': a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1,
        'w': a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1})
    """ 
    # If the degree is less than 12, the polynomial can collapse to a simpler form that does not fit the family structure
    if n < 12:
        return False, {}
    
    if n % 3 != 0:
        return False, {}
    
    k = n // 3
    if gcd(k, 3) != 1:
        return False, {}
    
    F = GF(2**n, 'a')
    K = F.subfield(k)
    R = PolynomialRing(F, 'x')
    x = R.gen()    

    terms = get_terms((poly.mod(x**(2**n) - x)))
    if not terms:
        return False, {}
    
    for s in range(1, n):
        if gcd(s, 3*k) != 1 or (k + s) % 3 != 0:
            continue

        e_u = 2**s + 1                                          # exponent of ux^(2^s + 1)
        e_ux = (2**(n - k) + 2**(k + s)) % (2**n - 1)           # exponent of u^(2^k) * x^(2^-k + 2^(k + s))
        e_v = (2**(n - k) + 1) % (2**n - 1)                     # exponent of vx^(2^-k + 1)
        e_wu = (2**s + 2**(k + s)) % (2**n - 1)                 # exponent of wu^(2^k + 1)

        if not set(terms).issubset({e_u, e_ux, e_v, e_wu}):
            continue

        u = terms.get(e_u,  F(0))
        uk = terms.get(e_ux, F(0))
        v = terms.get(e_v,  F(0))
        wu = terms.get(e_wu,  F(0))

        if u == F(0) or not is_primitive_element(F, u):
            continue
        if uk != u**(2**k):
            continue
        
        w = wu / (u**(2**k + 1))
        if v * w == K(1):
            continue

        return True, {'s': s, 'k': k, 'u': u, 'v': v, 'w': w}

    return False, {}


def belong_family10(n, poly):
    pass


def belong_family11(n, poly):
    r"""
    Check if the polynomial belongs to Family11, the Budaghyan-Helleseth-Kaleyski construction from 2020.
    Defined by `f(x) = x^3 + a(x^2^i + 1) + bx^(3 * 2^m) + c(x^(2^(i + m) + 2^m))^2^k`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import belong_family11
        sage: F.<a> = GF(2^10)
        sage: R.<x> = PolynomialRing(F)
        sage: poly =  (a^5 + a^3 + a)*x^576 + (a^5 + a^3 + a + 1)*x^96 + x^18 + x^3
        sage: belong_family11(10, poly)
        (True, {'k': 6, 'i': 3, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: poly = x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3
        sage: belong_family11(10, poly)
        (True, {'k': 2, 'i': 9, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: F.<a> = GF(2^14)
        sage: R.<x> = PolynomialRing(F)
        sage: poly = x^768 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a + 1)*x^384 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^6 + x^3
        sage: belong_family11(14, poly)
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
    K = F.subfield(2)
    R = PolynomialRing(F, 'x')
    x = R.gen()    

    terms = get_terms((poly.mod(x**(2**n) - x)))
    if not terms:
        return False, {}

    if terms.get(3, F(0)) != F(1):
        return False, {}
    
    for k in range(0, n):
        if k % 2 == 0:
            i_list = {m-2, m, n-1} | {i for i in range(1, n) if (i * (m - 2)) % n == 1}
        else:
            i_list = {m+2, m} | {i for i in range(1, n) if (i * (m + 2)) % n == 1}
        
        for i in i_list:
            e_3 = 3                                                 # exponent of x^3
            e_a = ((2**i + 1) * 2**k) % (2**n - 1)                  # exponent of a(x^2^i + 1)^2^k
            e_b = (3 * 2**m)                                        # exponent of bx^(3 * 2^m)
            e_c = ((2**(i + m) + 2**m) * 2**k) % (2**n - 1)         # exponent of c(x^(2^(i + m) + 2^m))^2^k
    
            if not set(terms).issubset({e_3, e_a, e_b, e_c}):
                continue

            if e_a == e_c:
                c = K(1)
                a = terms.get(e_a, K(0)) + K(1)
            else:
                a = terms.get(e_a, K(0))
                c = terms.get(e_c, K(0))
            
            b = terms.get(e_b, K(0))
            if b != a**2:
                continue
            if c != K(1):
                continue
       
            # a must be a primitive element of GF(2^2)
            if a == K(0) or not is_primitive_element(K, a):
                continue
            
            return True, {'k': k, 'i': i, 'a': a, 'b': b, 'c': c}
        
    return False, {}


def belong_family12(n, poly):
    r"""
    Check if the polynomial belongs to Family12, the Zheng-Kan-Li-Peng-Tang from 2022.
    Defined by `f(x) = a * Tr^n_m(bx^(2^i + 1)) + a^q * Tr^n_m(cx^(2^s + 1))`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    

    EXAMPLES::
    """
    # Placeholder for future implementation
    return False, {}


def belong_family13(n, poly):
    r"""
    Check if the polynomial belongs to Family13, the Li-Zhou-Li-Qu construction from 2022.
    Defined by `f(x) = L(z)^(2^m + 1) + cz^(2^m + 1)`.

    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::
    """
    # Placeholder for future implementation
    return False, {}


FAMILIES = {
    "Family1": belong_family1,
    "Family2": belong_family2,
    "Family3": belong_family3,
    "Family7_9": belong_family7_9,
    "Family11": belong_family11,
}

def belong(n, polynomial):
    r"""
    Check a function for membership in known infinite families of quadratic APN polynomials
    
    INPUT:

    - ``n`` -- the degree of the field GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: belong(6, x^3)
        Belong to Family1: False
        Belong to Family2: False
        Belong to Family3: False
        Belong to Family7_9: False
        Belong to Family11: False
    """
    for family_name, family_function in FAMILIES.items():
        try:
            found, params = family_function(n, polynomial)
        except Exception as e:
            print(f"Belong to {family_name}: {False}")
            continue
        print(f"Belong to {family_name}: {found}")
        if found:
            print(f"With params: {params}")
