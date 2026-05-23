from sage.all import *
from helpers import is_primitive_element, get_terms


def _membership_family1_2(n, p, poly):
    """
    Shared implementation for Family 1 (p=3) and Family 2 (p=4).
    """
    if n < 12 or n % p != 0 :
        return False, {}
      
    k = n // p
    if gcd(k, 3) != 1:
        return False, {}
    
    F = GF(2**n, 'a')
  
    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    for s_val in range(1, n):
        if gcd(s_val, 3 * k) != 1:
            continue
        i = (s_val * k) % p
        m = p - i

        e1 = 2**s_val + 1
        e2 = (2**(i * k) + 2**(m * k + s_val)) % (2**n - 1)

        # Single-term case: both exponents coincide
        if e1 == e2:
            if set(terms) != {e2}:
                continue
            # Combined coeff = 1 + u^(2^k-1), so u^(2^k-1) = coeff + 1
            u_power = terms.get(e2, F(0)) + F(1) 
        
        # Two term case: exponents are distinct
        else:
            if set(terms) != {e1, e2}:
                continue
            # Extract the coefficients for each term, defaulting to 0 if the term is not present
            if terms.get(e1, F(0)) != F(1):
                continue
            u_power = terms.get(e2, F(0))
        
        if u_power == F(0):
            continue
        
        # Recover u from u^(2^k - 1)
        try:
            all_roots = u_power.nth_root(2**k - 1, all=True)
        except ValueError:
            continue
        
        # u must be a primitive element of GF(2^n)
        u_vals = [u_val for u_val in all_roots if u_val != 0 and is_primitive_element(F, u_val)]
        if not u_vals:
            continue

        return True, {'s': s_val, 'u': u_vals}
    
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
        sage: polynomial = (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33
        sage: membership_family_1(12, polynomial)
        (True, {'s': 5, 'u': 
            [a^11 + a^8 + a^6 + a^3 + a + 1,
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
        
        sage: polynomial = x^129 + (a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^4)*x^24
        sage: membership_family_1(12, polynomial)
        (True, {'s': 5, 'u': 
            [a^11 + a^8 + a^6 + a^3 + a + 1,
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
        sage: polynomial = (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3
        sage: membership_family_2(16, polynomial)
        (True, {'s': 1, 'u': 
            [a^12 + a^8 + a^6 + a^5 + a^4 + a^3 + a^2,
            a^15 + a^14 + a^7 + a^6 + a^3 + a,
            a^11 + a^10 + a^9 + a^8 + a^7 + a^3 + a + 1,
            a^13 + a^11 + a^10 + a^7 + a^4 + a + 1,
            a^15 + a^14 + a^12 + a^11 + a^10 + a^9 + a^5 + a^4 + a^3 + a^2 + 1,
            a^15 + a^14 + a^13 + a^12 + a^9 + a^7 + a^5 + a^3 + a^2 + a,
            a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^8 + a^5 + a^2 + 1,
            a^13 + a^12 + a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a + 1]})

        sage: polynomial = (a^15 + a^14 + a^12 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2)*x^2049
        sage: membership_family_2(16, polynomial)
        (True, {'s': 11, 'u': 
            [a^14 + a^9 + a + 1,
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
        sage: polynomial = x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3
        sage: membership_family_3(6, polynomial)
        (True, {'i': 1, 's': a, 'c': a})  

        sage: polynomial = x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5
        sage: membership_family_3(6, polynomial)
        (True, {'i': 2, 's': a^2 + a + 1, 'c': a^5 + a^4 + a^2 + a + 1})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')
    K = GF(q)

    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    for i_val in range(1, m): 
        if gcd(i_val, m) != 1:
            continue

        e1 = q + 1
        e2 = 2**i_val + 1
        e3 = (q * (2**i_val + 1)) % (2**n - 1)
        e4 = (2**i_val * q + 1) % (2**n - 1)
        e5 = 2**i_val + q

        if not set(terms).issubset({e1, e2, e3, e4, e5}):
            continue

        if terms.get(e2) != F(1) or terms.get(e3) != F(1):
            continue

        s_val = terms.get(e1, F(0))
        c_val = terms.get(e4, F(0))
        cq = terms.get(e5, F(0))

        # Coeff at e_cq must equal c^q
        if c_val**q != cq:
            continue
        
        # s must be in GF(2^n) but not in the subfield GF(q)
        if s_val in K:
            continue

        return True, {'i': i_val, 's': s_val, 'c': c_val}

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
        sage: polynomial = (a^8 + a^3 + a)*x^288 + (a^8 + a^6 + a^5 + a^4 + a)*x^260 + (a^7 + a^4 + a^3 + a^2 + a)*x^144 + (a^8 + a^6 + a^2 + a)*x^130 + (a^6 + a^5 + a^4)*x^72 + (a^7 + a^6 + a^4 + a^3)*x^65 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^36 + (a^8 + a^4 + a^2 + 1)*x^18 + (a^5 + a^3 + a)*x^9 + x^3
        sage: membership_family_4(9, polynomial)
        (True, {'a': a^8 + a^6 + a^5 + a^3 + a})

        sage: polynomial = x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3
        sage: membership_family_4(9, polynomial)
        (True, {'a': 1})

        sage: F.<a> = GF(2^7)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^6 + a^5 + a^3 + a^2 + a + 1)*x^72 + (a^5 + a^4 + a^2 + 1)*x^68 + (a^5 + a^4 + a^3 + 1)*x^36 + (a^3 + a^2 + 1)*x^34 + (a^6 + a^5 + a^4 + a)*x^18 + (a^5 + a^3 + a^2 + a)*x^17 + (a^6 + a^4 + a^2 + 1)*x^9 + x^3
        sage: membership_family_4(7, polynomial)
        (True, {'a': a^3 + a^2 + a + 1})
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    terms = get_terms(n, poly)
    if not terms:
        return False, {}    

    # All exponents must lie in the orbit {9 * 2^i mod (2^n-1)} ∪ {3}
    e_9 = {(9 * 2**i) % (2**n - 1) for i in range(n)}
    if not set(terms).issubset({3} | e_9):
        return False, {}

    # Coeff of x^9 is a^2; recover a via square root
    a_square = terms.get(9 % (2**n - 1), F(0))
    if a_square == F(0):
        return False, {}
    a_val = a_square**(2**(n - 1))
    
    # Verify the reduced polynomial matches the expected form
    trace = sum(a_val**(3 * 2**i) * x**((9 * 2**i) % (2**n - 1)) for i in range(n))
    expected = (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)
    reduced = R(poly).mod(x**(2**n) - x)

    if reduced == expected:
        return True, {'a': a_val}
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
        sage: polynomial = (a^8 + a^4 + a)*x^144 + (a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^130 + (a^6 + a^4 + a^3 + a^2 + a + 1)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^65 + (a^7 + a^3 + a^2)*x^18 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3)*x^9 + x^3
        sage: membership_family_5(9, polynomial)
        (True, {'a': a^7 + a^6 + a^4 + a^3 + 1})

        sage: polynomial = x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3
        sage: membership_family_5(9, polynomial)
        (True, {'a': 1})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^11 + a^10 + a^9 + a^8 + a^7 + a^3 + a^2)*x^1152 + (a^10 + a^9 + a^7 + a^3 + a^2)*x^1026 + (a^11 + a^8 + a^6 + a^3 + 1)*x^576 + (a^11 + a^10 + a^6 + a^5 + a^4 + a)*x^513 + (a^10 + a^8 + a^5 + a^3 + a^2 + 1)*x^144 + (a^9 + a^8 + a^5 + a^4 + a^3 + 1)*x^72 + (a^9 + a^8 + a^7 + a^5 + a^3 + 1)*x^18 + (a^11 + a^5 + a^4 + a^2 + 1)*x^9 + x^3
        sage: membership_family_5(12, polynomial)
        (True, {'a': a^10 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2})
    """
    if n % 3 != 0:
        return False, {}
    
    k = n // 3

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    # All exponents must lie in the orbit {9 * 2^(3*i) mod (2^n-1)} ∪ {18 * 2^(3*i) mod (2^n-1)} ∪ {3} 
    e_9 = {(9 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    e_18 = {(18 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    if not set(terms).issubset({3} | e_9 | e_18):
        return False, {}
    
    # Coeff of x^9 is a^2; recover a via square root
    a_square = terms.get(9 % (2**n - 1), F(0))
    if a_square == F(0):
        return False, {}
    a_val = a_square**(2**(n - 1))

    # Verify the reduced polynomial matches the expected form
    trace = sum(a_val**(3 * 2**(3*i)) * x**((9 * 2**(3*i)) % (2**n - 1)) + (a_val**(6 * 2**(3*i)) * x**((18 * 2**(3*i)) % (2**n - 1))) for i in range(k))
    expected = (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)
    reduced = R(poly).mod(x**(2**n) - x)       

    if reduced == expected:
        return True, {'a': a_val}
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
        sage: polynomial = (a^7 + a^5 + a^3)*x^288 + (a^8 + a^4 + a^3)*x^260 + (a^6 + a^5 + a^3 + a^2)*x^144 + (a^7 + a^2 + 1)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^36 + (a^8 + a^4 + a^3 + a^2 + 1)*x^18 + x^3
        sage: membership_family_6(9,  polynomial)
        (True, {'a': a^8 + a^7 + a^2})

        sage: polynomial = x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3
        sage: membership_family_6(9,  polynomial)
        (True, {'a': 1})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^2304 + (a^11 + a^7 + a^6 + a^5 + a^4 + a^3)*x^2052 + (a^10 + a^7 + a^6 + a^4 + a + 1)*x^1152 + (a^10 + a^8 + a^7 + a^6 + a^3 + 1)*x^1026 + (a^11 + a^8 + a^5 + a^2 + a)*x^288 + (a^10 + a^8 + a^6 + a^5 + a^4 + a^3 + a^2)*x^144 + (a^9 + a^8 + a^5 + a^4 + a^3 + a^2 + a)*x^36 + (a^9 + a^7 + a^6 + a^5 + a^3 + a^2)*x^18 + x^3
        sage: membership_family_6(12, polynomial)
        (True, {'a': a^11 + a^10 + a^9 + a^8 + a^7 + a^5 + a^4 + a})
    """
    if n % 3 != 0:
        return False, {}
    
    k = n // 3

    F = GF(2**n, 'a')

    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    # x^3 must be present
    if 3 not in terms:
        return False, {}
    
    # All exponents must lie in the orbit {18 * 2^(3*i) mod (2^n-1)} ∪ {36 * 2^(3*i) mod (2^n-1)} ∪ {3} 
    e_18 = {(18 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    e_36 = {(36 * 2**(3*i)) % (2**n - 1) for i in range(k)}
    if not set(terms).issubset({3} | e_18 | e_36):
        return False, {}
    
    # Coeff of x^18 is a^(6-1) = a^5; recover a via fifth root
    a_fifth = terms.get(18 % (2**n - 1), F(0))
    if a_fifth == F(0):
        return False, {}
    try:
        all_roots = a_fifth.nth_root(5, all=True)
    except ValueError:
        return False, {}
    
    # Verify all fifth roots to find a valid a
    for a_val in all_roots:
        if a_val == F(0):
            continue
        if all(
            terms.get((18 * 2**(3*i)) % (2**n - 1), F(0)) == a_val ** (6  * 2**(3*i) - 1) and
            terms.get((36 * 2**(3*i)) % (2**n - 1), F(0)) == a_val ** (12 * 2**(3*i) - 1)
            for i in range(k)
        ):
            return True, {'a': a_val}

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
        sage: polynomial = (a^10 + a^6 + a^5 + a^3 + a)*x^768 + (a^9 + a^8 + a^6 + a^5 + a^2)*x^544 + (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^257 + (a^7 + a^5 + a^3 + a^2 + 1)*x^33
        sage: membership_family_7_9(12, polynomial)
        (True, {'s': 5, 'u': a^7 + a^5 + a^3 + a^2 + 1, 'v': a4^3 + a4, 'w': a4^3 + a4})

        sage: polynomial = a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + x^257
        sage: membership_family_7_9(12, polynomial)
        (True, {'s': 11, 'u': a, 'v': 1, 'w': 0})

        sage: polynomial = (a^10 + a^3 + a^2 + 1)*x^2056 + a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257
        sage: membership_family_7_9(12, polynomial)
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
    
    for s_val in range(1, n):
        if gcd(s_val, 3*k) != 1 or (k + s_val) % 3 != 0:
            continue

        e1 = (2**s_val + 1)
        e2 = (2**(n - k) + 2**(k + s_val)) % (2**n - 1)
        e3 = (2**(n - k) + 1) % (2**n - 1)
        e4 = (2**s_val + 2**(k + s_val)) % (2**n - 1)

        if not set(terms).issubset({e1, e2, e3, e4}):
            continue

        u_val = terms.get(e1,  F(0))
        uk = terms.get(e2, F(0))
        v_val = terms.get(e3,  F(0))
        wu = terms.get(e4,  F(0))

        # u must be primitive in GF(2^n)
        if u_val == F(0) or not is_primitive_element(F, u_val):
            continue
        # Coeff at e2 must equal u^(2^k)
        if uk != u_val**(2**k):
            continue
        
        # Recover w from w * u^(2^k+1) and verify subfield membership and vw ≠ 1
        w_val = wu / (u_val**(2**k + 1))
        if v_val * w_val == K(1) or v_val not in K or w_val not in K:
            continue

        return True, {'s': s_val, 'u': u_val, 'v': K(v_val), 'w': K(w_val)}

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
        sage: polynomial =  (a^5 + a^3 + a)*x^576 + (a^5 + a^3 + a + 1)*x^96 + x^18 + x^3
        sage: membership_family_11(10, polynomial)
        (True, {'k': 6, 'i': 3, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: polynomial = x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3
        sage: membership_family_11(10, polynomial)
        (True, {'k': 2, 'i': 9, 'a': a^5 + a^3 + a, 'b': a^5 + a^3 + a + 1, 'c': 1})

        sage: F.<a> = GF(2^14)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^768 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a + 1)*x^384 + (a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^6 + x^3
        sage: membership_family_11(14, polynomial)
        (True, {'k': 2, 'i': 13, 'a': a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a, 'b': a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a + 1, 'c': 1})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    if m % 2 == 0 or m % 3 == 0:
        return False, {}
    
    F = GF(2**n, 'a')
    K = F.subfield(2)

    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    # x^3 must be present
    if 3 not in terms:
        return False, {}
    
    def _valid_i(k_val):
        p = (m - 2) if k_val % 2 == 0 else (m + 2)
        base = {m - 2, m, n - 1} if k_val % 2 == 0 else {m + 2, m}
        inv_set = {pow(p, -1, n)} if gcd(p, n) == 1 else set()
        return tuple(base | inv_set)
    
    pair_ki = [(k_val, i_val) for k_val in range(n) for i_val in _valid_i(k_val)]

    for k_val, i_val in pair_ki:
        e1 = ((2**i_val + 1) * 2**k_val) % (2**n - 1)
        e2 = (3 * 2**m)
        e3 = ((2**(i_val + m) + 2**m) * 2**k_val) % (2**n - 1)

        if not set(terms).issubset({3, e1, e2, e3}):
            continue
        
        # Single-term case: a and c exponents coincide
        if e1 == e3:
            a_val = terms.get(e1, F(0)) + F(1)
            c_val = F(1)
        # Two-term case: exponents are distinct
        else:
            a_val = terms.get(e1, F(0))
            c_val = terms.get(e3, F(0))
        
        b_val = terms.get(e2, F(0))
        # Require b = a^2 and c = 1
        if b_val != a_val**2 or c_val != F(1):
            continue
    
        # a must be a primitive element of GF(4)
        if a_val == F(0) or not is_primitive_element(K, a_val):
            continue
        
        return True, {'k': k_val, 'i': i_val, 'a': a_val, 'b': b_val, 'c': c_val}
        
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
        sage: polynomial = (a^8 + a^6 + a^4 + a^3 + 1)*x^96 + (a^9 + a^7 + a^6 + a^4)*x^33 + (a^8 + a^6 + a^5 + a^3)*x^3
        sage: membership_family_12(10, polynomial)
        (True, {'i': 1, 's': 5, 'a': a^8 + a^5 + a^4 + a^3 + a^2, 'b': a^8 + a^6 + a^5 + a^4 + a^3 + a + 1, 'c': a})

        sage: polynomial = a^3*x^513 + (a^8 + a^7 + a^2 + a + 1)*x^48 + (a^3 + a^2 + 1)*x^33
        sage: membership_family_12(10, polynomial)
        (True, {'i': 9, 's': 5, 'a': a^9 + a^6 + a^2 + a + 1, 'b': a^6 + a^5 + a^3, 'c': a^9 + a^7 + a^6 + a^4 + a^3 + a})
    """
    if n % 2 != 0:
        return False, {}
    
    m = n // 2
    q = 2**m
    if m % 2 == 0:
        return False, {}
    
    F = GF(2**n, 'a')
    Q = GF(q)

    terms = get_terms(n, poly)
    if not terms:
        return False, {}
    
    def _terms_match(i_val, s_val, a_val, b_val, c_val):
        # Build the expected coefficient dict and compare with terms
        e1 = (2**i_val + 1)
        e2 = (q * e1) % (2**n - 1)
        e3 = (2**s_val + 1)
        e4 = (q * e3) % (2**n - 1)
        expected = {}
        for e, coeff in [(e1, a_val*b_val), (e2, a_val*b_val**q), (e3, a_val**q*c_val), (e4, a_val**q*c_val**q)]:
            expected[e] = expected.get(e, F(0)) + coeff
        return {e: v for e, v in expected.items() if v != F(0)} == terms
        
    for i_val in range(1, n):
        if gcd(i_val, n) != 1:
            continue

        for s_val in range(1, n):
            e1 = (2**i_val + 1)
            e2 = (q * e1) % (2**n - 1)
            e3 = (2**s_val + 1)
            e4 = (q * e3) % (2**n - 1)

            if not set(terms).issubset({e1, e2, e3, e4}):
                continue

            c1 = terms.get(e1, F(0))
            c2 = terms.get(e2, F(0))

            # Single-term case: i and s exponents coincide
            if i_val == s_val:
                if c1 == F(0) and c2 == F(0):
                    continue
                for a_val in F:
                    if a_val == F(0) or a_val in Q or a_val + a_val**q == F(0):
                        continue
                    # c^q = (c2 + c1^q * a^(1-q)) / (a^(2-q) + a^q)
                    denom = a_val**(2 - q) + a_val**q
                    if denom == F(0):
                        continue
                    cq = (c2 + c1**q * a_val**(1 - q)) / denom
                    if cq== F(0):
                        continue
                    # c = (c^q)^q since (c^q)^q = c^(q^2) = c
                    c_val = cq**q
                    b_val = (c1 + a_val**q * c_val) / a_val
                    if b_val == F(0):
                        continue

                    if _terms_match(i_val, s_val, a_val, b_val, c_val):
                        return True, {'i': i_val, 's': s_val, 'a': a_val, 'b': b_val, 'c': c_val}
        
            # Two-term case: i and s exponents are distinct
            else:
                if c1 == F(0) or c2 == F(0):
                    continue
                c3 = terms.get(e3, F(0))
                try:
                    all_roots = (c1**q / c2).nth_root(q - 1, all=True)
                except ValueError:
                    continue
            
                for a_val in all_roots:
                    if a_val == F(0) or a_val in Q or a_val + a_val**q == F(0):
                        continue
                
                    b_val = c1 / a_val
                    if b_val == F(0) or a_val*b_val**q != c2:
                        continue

                    c_can = c3 / a_val**q
                    if e3 == e4:
                        c_vals = [c for c in F if c != F(0) and c**q + c == c_can]
                    elif c_can == F(0):
                        continue
                    else:
                        c_vals = [c_can]

                    for c_val in c_vals:
                        if c_val == F(0):
                            continue
                        if _terms_match(i_val, s_val, a_val, b_val, c_val):
                            return True, {'i': i_val, 's': s_val, 'a': a_val, 'b': b_val, 'c': c_val}
    
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
        sage: polynomial = x^144 + (a^8 + a^7 + a^3 + a^2 + a)*x^130 + x^129 + (a^8 + a^2)*x^32 + x^24 + (a^7 + a^6 + a^4 + a^3 + a)*x^18 + (a^8 + a^2)*x^17 + (a^8 + a^7 + a^3 + a^2 + a)*x^10 + (a^4 + a^3 + a^2 + 1)*x^9
        sage: membership_family_13(9, polynomial)
        (True, {'s': 1, 'v': a^4 + a^3 + a^2, 'mu': a^8 + a^7 + a^3 + a^2 + a})

        sage: polynomial = x^288 + (a^8 + a^7 + a^6 + a^3 + a^2)*x^260 + x^257 + (a^8 + a^7 + a^5 + a^3 + a^2 + a + 1)*x^64 + x^40 + (a^7 + a^3 + a^2)*x^36 + (a^8 + a^7 + a^5 + a^3 + a^2 + a + 1)*x^33 + (a^8 + a^7 + a^6 + a^3 + a^2)*x^12
        sage: membership_family_13(9, polynomial)
        (True, {'s': 2, 'v': 1, 'mu': a^8 + a^7 + a^6 + a^3 + a^2})

        sage: F.<a> = GF(2^12)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^544 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^514 + x^513 + (a^7 + a^6 + a^4 + a^3 + a)*x^64 + x^48 + (a^11 + a^9 + a^5)*x^34 + (a^7 + a^6 + a^4 + a^3 + a)*x^33 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^18 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^17
        sage: membership_family_13(12, polynomial)
        (True, {'s': 1, 'v': a^11 + a^9 + a^8 + a^6 + a^3 + a, 'mu': a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1})
    """
    if n % 3 != 0:
        return False, {}
    
    m = n // 3

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    M = F.subfield(m)

    terms = get_terms(n, poly)
    if not terms:
        return False, {}

    for s_val in range(1, m + 1):
        if gcd(s_val, m) != 1:
            continue

        P, Q, Rp, B = 2**(2*m + s_val), 2**(m + s_val), 2**m, 2**s_val

        ePQ = (P + Q) % (2**n - 1)
        ePB = (P + B) % (2**n - 1)
        eP1 = (P + 1) % (2**n - 1)
        e2Q = (2 * Q) % (2**n - 1)
        eQB = (Q + B) % (2**n - 1)
        eQ1 = (Q + 1) % (2**n - 1)
        eRQ = (Rp + Q) % (2**n - 1)
        eRB = (Rp + B) % (2**n - 1)
        eR1 = (Rp + 1) % (2**n - 1)

        expected = {ePQ, ePB, eP1, e2Q, eQB, eQ1, eRQ, eRB, eR1}
        if not set(terms).issubset(expected):
            continue

        # Recover mu from x^(P+B)
        mu_val = terms.get(ePB, F(0))
        if mu_val == F(0):
            continue

        if terms.get(ePQ, F(0)) != F(1) or terms.get(eP1, F(0)) != F(1) or terms.get(eRQ, F(0)) != F(1):
            continue
        if terms.get(eRB, F(0)) != mu_val:
            continue
        mu_q = mu_val**(2**m)
        if terms.get(e2Q, F(0)) != mu_q or terms.get(eQ1, F(0)) != mu_q:
            continue
        if terms.get(eQB, F(0)) != mu_q * mu_val:
            continue

        # Recover v from x^(R+1) coefficient = 1 + v
        v_val = terms.get(eR1, F(0)) + F(1)
        if v_val == F(0) or v_val not in M:
            continue

        if mu_val**(2**(2*m) + 2**m + 1) == F(1):
            continue

        # L permutes F
        L = x**(2**(m + s_val)) + mu_val*x**(2**s_val) + x
        if L.gcd(x**(2**n) - x) != x:
            continue

        return True, {'s': s_val, 'v': v_val, 'mu': mu_val}

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
        Belongs to Family 2, with parameters: {'s': 11, 'u': 
            [a^14 + a^13 + a^12 + a^11 + a^7 + a^6 + a^3 + a^2 + a, 
            ...
            a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^4]}
        
        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3
        sage: membership_all(6, polynomial)
        Belongs to Family 3, with parameters: {'i': 1, 's': a, 'c': a}
    """
    for name, function in FAMILIES.items():
        try:
            found, params = function(n, polynomial)
        except Exception:
            continue
        if found:
            print(f"Belongs to {name}, with parameters: {params}")
