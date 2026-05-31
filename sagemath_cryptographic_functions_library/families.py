from sage.all import *
from .helpers import aggregate_results, is_primitive_element, family12_s_candidates, family12_validates_s


def _family_1_2(n, p, s, u):
    """
    Shared implementation for Family 1 (p=3) and Family 2 (p=4).
    """
    if n < 12:
        raise TypeError("n must be at least 12")

    if n % p != 0:
        raise TypeError(f"n must be divisible by p = {p}")
    
    k = n // p
    if gcd(k, 3) != 1:
        raise TypeError("gcd(k, 3) must be 1")

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    if u is None: 
        g = F.multiplicative_generator()
        u_vals = (g**i for i in range(2**n - 1) if gcd(i, 2**n - 1) == 1)
    elif not is_primitive_element(F, u):
        raise TypeError("u must be a primitive element of GF(2^n)")
    else:
        u_vals = (u,)

    if s is None:
        s_vals = tuple(s_val for s_val in range(1, n) if gcd(s_val, 3 * k) == 1)
    elif gcd(s, 3*k) != 1:
        raise TypeError("gcd(s, 3*k) must be 1")
    else:
        s_vals = (s,)

    def _poly(s_val, u_val):
        i = (s_val * k) % p
        m = p - i
        e1 = (2**(i * k) + 2**(m * k + s_val)) % (2**n - 1)
        return (x**(2**s_val + 1) + u_val**(2**k - 1) * x**e1)
    
    pairs = (
        (_poly(s_val, u_val), {'s': s_val, 'u': u_val}) 
        for u_val in u_vals
        for s_val in s_vals
    )
    return aggregate_results(pairs)


def family_1(n, s=None, u=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2008 for `p = 3`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and gcd(k, 3) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_1
        sage: F.<a> = GF(2^12)
        sage: family_1(12, 5, a^2)
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33

        sage: family_1(12, None, a^6 + a^5 + 1)
        {(a^11 + a^9 + a^7 + a^6 + a^2)*x^528 + x^3: [{'s': 1, 'u': a^6 + a^5 + 1}],
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33: [{'s': 5, 'u': a^6 + a^5 + 1}],
        x^129 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^24: [{'s': 7, 'u': a^6 + a^5 + 1}],
        x^2049 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^264: [{'s': 11, 'u': a^6 + a^5 + 1}]}

        sage: family_1(12, 1)
        {(a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^528 + x^3: 
            [{'s': 1, 'u': a},
            {'s': 1, 'u': a^11 + a^10 + a^9 + a^5 + a^4 + a^3},
            ...
            {'s': 1, 'u': a^9 + a^7 + a^6 + a^5 + a^3}],
        ...
        (a^10 + a^9 + a^8 + a^5 + a^4 + a^3 + 1)*x^528 + x^3: 
            [{'s': 1, 'u': a^11 + a^10 + a^9 + a^7 + a^3 + a^2},
            ...
            {'s': 1, 'u': a^11 + a^10 + a^9 + a^8 + a^5 + a^3 + a^2 + a}]}

        sage: result = family_1(12); list(result.keys())
        [(a^9 + a^6 + a^5 + a^4 + a^2 + 1)*x^768 + x^33,
        (a^11 + a^10 + a^8 + a^7 + a^6 + a^3 + a^2 + a)*x^768 + x^33,
        ...
        x^129 + (a^11 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^24]

        sage: len(result)
        576
    """
    return _family_1_2(n, 3, s, u)


def family_2(n, s=None, u=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2008 for `p = 4`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 4k
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and gcd(k, 3) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_2
        sage: F.<a> = GF(2^16)
        sage: family_2(16, 5, a^15 + a^14 + a^13 + a^11 + a^3 + a^2 + a)
        (a^15 + a^14 + a^13 + a^12 + a^10 + a^9 + a^8 + a^4 + a^2 + a + 1)*x^33

        sage: family_2(16, None, a^15 + a^14 + a^7 + a^6 + a^3 + a)
        {(a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3: 
            [{'s': 1, 'u': a^15 + a^14 + a^7 + a^6 + a^3 + a}],
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^33: 
            [{'s': 5, 'u': a^15 + a^14 + a^7 + a^6 + a^3 + a}],
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^129: 
            [{'s': 7, 'u': a^15 + a^14 + a^7 + a^6 + a^3 + a}],
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^2049: 
            [{'s': 11, 'u': a^15 + a^14 + a^7 + a^6 + a^3 + a}],
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^8193: 
            [{'s': 13, 'u': a^15 + a^14 + a^7 + a^6 + a^3 + a}]}
 
        sage: family_2(16, 1)
        {(a^15 + 1)*x^3: 
            [{'s': 1, 'u': a},
            ...
            {'s': 1, 'u': a^13 + a^11 + a^10 + a^6 + a^5 + a^4 + a^3 + a^2}],
        ...
         (a^14 + a^12 + a^11 + a^10 + a^9 + a^7 + a^5 + a^4 + a^3 + a)*x^3:
            [{'s': 1, 'u': a^13 + a^12 + a^11 + a^8 + a^4 + a},
            ...
            {'s': 1, 'u': a^15 + a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a}]}
        
        sage: result = family_2(16); list(result.keys())
        [(a^15 + 1)*x^3,
        (a^15 + 1)*x^33,
        ...
        (a^14 + a^12 + a^11 + a^10 + a^9 + a^7 + a^5 + a^4 + a^3 + a)*x^8193]
        
        sage: len(result)
        20480
    """
    return _family_1_2(n, 4, s, u)


def family_3(n, i=None, s=None, c=None):
    r"""
    Return the Budaghyan-Carlet construction from 2008.
    Defined by `f(x) = sx^(q + 1) + x^(2^i + 1) + x^(q * (2^i + 1)) + cx^(2^i * q + 1) + c^q * x^(2^i + q))`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be even 
    - ``i`` -- (optional) positive integer with gcd(i, n/2) = 1; if None, returns a list over all valid i in {1, ... , m - 1}
    - ``s`` -- (optional) element of GF(2^n) not in GF(q); if None, returns a list over all valid s in GF(2^n) \ GF(q)
    - ``c`` -- (optional) element of GF(2^n); if None, returns a list over all valid c in GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_3
        sage: F.<a> = GF(2^6)
        sage: family_3(6, 1, a, a)
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3

        sage: family_3(6, None, a^2 + a + 1, a^5 + a^4 + a^2 + a + 1)
        {x^24 + (a^5 + a^4 + a^2 + a + 1)*x^17 + a*x^10 + (a^2 + a + 1)*x^9 + x^3: 
            [{'i': 1, 's': a^2 + a + 1, 'c': a^5 + a^4 + a^2 + a + 1}],
        x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5: 
            [{'i': 2, 's': a^2 + a + 1, 'c': a^5 + a^4 + a^2 + a + 1}]}
        
        sage: family_3(6, 1)
        {x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + x^3: 
            [{'i': 1, 's': 0, 'c': a}],
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3: 
            [{'i': 1, 's': a, 'c': a}],
        ...
        x^24 + (a^3 + a^2)*x^17 + (a^5 + a^4 + a^3)*x^10 + x^9 + x^3: 
            [{'i': 1, 's': 1, 'c': a^3 + a^2}]}

        sage: result = family_3(6); list(result.keys())
        [x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + x^3,
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3,
        ...
        x^40 + (a^5 + a^3 + a^2 + a + 1)*x^33 + (a^4 + a^3 + a + 1)*x^12 + x^9 + x^5]

        sage: len(result)
        2304
    """
    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    K = GF(q)
    g = F.gen()**(q - 1)

    if c is None:
        c_vals = iter(F)
    elif c not in F:
        raise TypeError("c must be an element of GF(2^n)")
    else:
        c_vals = (c,)
    
    if s is None:
        s_vals = tuple(s_val for s_val in F if s_val not in K)
    elif s in K or s not in F:
        raise TypeError("s must be in GF(2^n) but not in GF(q)")
    else:
        s_vals = (s,)

    if i is None:
        i_vals = (i_val for i_val in range(1, m) if gcd(i_val, m) == 1)
    elif gcd(i, m) != 1 or not (0 < i < m):
        raise TypeError("gcd(i, m) must be 1")
    else:
        i_vals = (i,)

    def _verify_root(i_val, c_val):
        """
        Check that X^(2^i + 1) + c * X^(2^i) + c^q * X + 1 = 0 has no solutions in GF(2^n) such that x^(q + 1) = 1.
        Iterate through all q+1 elements of U, generated by g = F.gen()^(q-1), and reject (i, c) if any is a root.
        """
        cq = c_val**q
        X = F(1)
        i2 = 2**i_val
    
        for _ in range(q + 1):
            Xi2 = X**(i2) # X^(2^i)
            if X * Xi2 + c_val * Xi2 + cq * X + 1 == 0:
                return False
            X *= g # next element of U
        return True
    
    ic_pair = ((i_val, c_val) for i_val in i_vals for c_val in c_vals if _verify_root(i_val, c_val))
    
    def _poly(i_val, s_val, c_val):
        e1 = (q * (2**i_val + 1)) % (2**n - 1)
        e2 = (2**i_val * q + 1) % (2**n - 1)
        return (s_val * x**(q + 1) + x**(2**i_val + 1) + x**e1 + c_val * x**e2 + c_val**q * x**(2**i_val + q))
    
    pairs = (
        (_poly(i_val, s_val, c_val), {'i': i_val, 's': s_val, 'c': c_val})
        for i_val, c_val in ic_pair
        for s_val in s_vals
    )
    return aggregate_results(pairs)


def family_4(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr_n(a^3 * x^9)`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_4
        sage: F.<a> = GF(2^9)
        sage: family_4(9, a^5 + 1)
        (a^2 + a)*x^288 + (a^7 + a^2)*x^260 + (a^7 + a^3 + a^2 + 1)*x^144 + a*x^130 + (a^8 + a^3 + a^2)*x^72 + (a^7 + a^2 + 1)*x^65 + (a^8 + a^6 + a)*x^36 + (a^7 + a^6 + a^5 + a^2 + 1)*x^18 + (a^5 + a + 1)*x^9 + x^3
        
        sage: family_4(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3

        sage: result = family_4(9); list(result.keys())
        [x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + x^3,
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3,
        ...
        x^24 + (a^5 + a^4 + a^3)*x^17 + (a^3 + a^2)*x^10 + (a^5 + a^3 + a^2 + 1)*x^9 + x^3]

        sage: len(result)
        511
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    exponents = tuple((3 * 2**i, (9 * 2**i) % (2**n - 1)) for i in range(n))
  
    def _poly(a_val):
        trace = sum(a_val**e_a * x**e_x for e_a, e_x in exponents)
        return x**3 + (F(1) / a_val) * trace
    
    if a is None:
        a_vals = (a_val for a_val in F if a_val != 0)
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    else:
        a_vals = (a,)
    
    pairs = (
        (_poly(a_val), {'a': a_val})
        for a_val in a_vals
    )
    return aggregate_results(pairs)


def family_5(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^3 * x^9 + a^6 * x^18)`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_5
        sage: F.<a> = GF(2^9)
        sage: family_5(9, a^6 + a^5)
        (a^3 + a^2 + 1)*x^144 + (a^8 + a^5 + a^3 + a^2)*x^130 + (a^4 + a^3 + a + 1)*x^72 + (a^3 + a + 1)*x^65 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a)*x^9 + x^3
        
        sage: family_5(9, F(1))
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3

        sage: result = family_5(9); list(result.keys())
        [(a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^8 + a^5 + a^4)*x^72 + (a^8 + a^7 + a^4 + a)*x^65 + a^5*x^18 + a^2*x^9 + x^3,
        (a^8 + a^5 + a^3 + a^2)*x^144 + (a^8 + a^6 + a^4 + a^3 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^2 + a)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + 1)*x^65 + (a^5 + a)*x^18 + a^4*x^9 + x^3,
        ...
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3]
        
        sage: len(result)
        511
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
    
    k = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    exponents = tuple((3 * 2**(3*i),  (9  * 2**(3*i)) % (2**n - 1), 6 * 2**(3*i),  (18 * 2**(3*i)) % (2**n - 1)) for i in range(k))

    def _poly(a_val):
        trace = sum(a_val**e_a1 * x**e_x1 + a_val**e_a2 * x**e_x2 for e_a1, e_x1, e_a2, e_x2 in exponents)
        return x**3 + (F(1) / a_val) * trace

    if a is None:
        a_vals = (a_val for a_val in F if a_val != 0)
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    else:
        a_vals = (a,)
    
    pairs = (
        (_poly(a_val), {'a': a_val})
        for a_val in a_vals
    )
    return aggregate_results(pairs)


def family_6(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^6 * x^18 + a^12 * x^36)`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_6
        sage: F.<a> = GF(2^9)
        sage: family_6(9, a^8 + a^5 + a^3 + 1)
        (a^7 + a^6 + 1)*x^288 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + a)*x^260 + (a^5 + a^4)*x^144 + (a^7 + a^5 + a^4 + a^3 + a^2)*x^130 + (a^8 + a^6 + a^4 + a^2 + a)*x^36 + (a^7 + a^6 + a^2 + 1)*x^18 + x^3

        sage: family_6(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3

        sage: result = family_6(9); list(result.keys())
        [(a^6 + a^3 + 1)*x^288 + (a^7 + a^5 + a^2 + 1)*x^260 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^6 + a^2)*x^36 + a^5*x^18 + x^3,
        (a^7 + a^6 + a^3 + 1)*x^288 + a*x^260 + (a^8 + a^5 + a^3 + a^2)*x^144 + (a^8 + a^6 + a^4 + a^3 + a)*x^130 + (a^7 + a^4 + a^3)*x^36 + (a^5 + a)*x^18 + x^3,
        ...
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3]

        sage: len(result)
        511
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
    
    k = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    exponents = tuple((6  * 2**(3*i), (18 * 2**(3*i)) % (2**n - 1), 12 * 2**(3*i), (36 * 2**(3*i)) % (2**n - 1)) for i in range(k))

    def _poly(a_val):
        trace = sum(a_val**e_a1 * x**e_x1 + a_val**e_a2 * x**e_x2 for e_a1, e_x1, e_a2, e_x2 in exponents)
        return x**3 + (F(1) / a_val) * trace

    if a is None:
        a_vals = (a_val for a_val in F if a_val != 0)
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    else:
        a_vals = (a,)
    
    pairs = (
        (_poly(a_val), {'a': a_val})
        for a_val in a_vals
    )
    return aggregate_results(pairs)


def family_7_9(n, s=None, u=None, v=None, w=None):
    r"""
    Return the Bracken-Byrne-Markin-McGuire construction from 2011.
    Defined by `f(x) = ux^(2^s + 1) + u^(2^k) * x^(2^-k + 2^(k + s)) + vx^(2^-k + 1) + wu^(2^k + 1) * x^(2^s + 2^(k + s))`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k, gcd(k, 3) = 1
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and 3|(k + s); if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)
    - ``v`` -- (optional) element of GF(2^k) with v*w != 1; if None, returns a list for all valid v in GF(2^k)
    - ``w`` -- (optional) element of GF(2^k) with v*w != 1; if None, returns a list for all valid w in GF(2^k)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_7_9
        sage: F.<a> = GF(2^12)
        sage: K = F.subfield(4)
        sage: family_7_9(12, 5, a, K(1), K(0))
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + x^257 + a*x^33
        
        sage: v, w = K.random_element(), K.random_element(); v, w
        (a4^3 + a4 + 1, a4)
        sage: family_7_9(12, 11, None, v, w)
        {(a^11 + a^10 + a^9 + a^6 + a^4 + a^2 + a)*x^2056 + a*x^2049 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^264 + (a^10 + a^9 + a^6 + a^5 + a^3)*x^257: 
            [{'s': 11, 'u': a, 'v': a4^3 + a4 + 1, 'w': a4}],
        ...
        (a^10 + a^7 + a^3 + a^2)*x^2056 + (a^11 + a^6 + a^5 + a^4 + a^2 + 1)*x^2049 + (a^10 + a^9 + a^7 + a^3 + a^2 + a)*x^264 + (a^10 + a^9 + a^6 + a^5 + a^3)*x^257: 
            [{'s': 11, 'u': a^11 + a^6 + a^5 + a^4 + a^2 + 1, 'v': a4^3 + a4 + 1, 'w': a4}]}

        sage: family_7_9(12, 5, a)
        {((a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^10 + a^9 + a^3)*x^544 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257 + a*x^33,
            [{'s': 5, 'u': a, 'v': a4^2, 'w': a4^3 + a4^2}])
        ...
        ((a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^10 + a^9 + a^8 + a^7 + a^6 + a^4)*x^544 + (a^8 + a^6 + a^5 + a^4 + a^2 + 1)*x^257 + a*x^33,
            [{'s': 5, 'u': a, 'v': a4^3, 'w': a4^2 + a4 + 1}])}

        sage: result = family_7_9(12, 5); list(result.keys())
        [(a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + a*x^33,
        (a^11 + a^9 + a^7 + a^6 + a^2 + a)*x^768 + a^2*x^33
        ...
        (a^10 + a^9 + a^7 + a^3 + a^2 + a)*x^768 + (a^7 + a + 1)*x^544 + x^257 + (a^11 + a^6 + a^5 + a^4 + a^2 + 1)*x^33]
        
        sage: len(result)
        416448
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
 
    k = n // 3
    if gcd(k, 3) != 1:
        raise TypeError("k must be coprime to 3")

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    K = F.subfield(k)

    if u is None: 
        g = F.multiplicative_generator()
        u_vals = (g**i for i in range(2**n - 1) if gcd(i, 2**n - 1) == 1)
    elif not is_primitive_element(F, u) or u == 0:
        raise TypeError("u must be a primitive element of GF(2^n)")
    else:
        u_vals = (u,)

    if v is not None and v not in K:
        raise TypeError("v must be an element of GF(2^k)")
    if w is not None and w not in K:
        raise TypeError("w must be an element of GF(2^k)")

    if v is not None and w is not None:
        if v * w == K(1):
            raise TypeError("v and w must satisfy v*w != 1")
        pair_vw = ((v, w),)
    else:
        v_vals = K if v is None else (v,)
        w_vals = K if w is None else (w,)
        pair_vw = tuple((v_val, w_val) for v_val in v_vals for w_val in w_vals if v_val * w_val != K(1))

    if s is None:
        s_vals = tuple(s_val for s_val in range(1, n) if (gcd(s_val, 3*k) == 1 and (k + s_val) % 3 == 0))
    elif gcd(s, 3*k) != 1 or (k + s) % 3 != 0:
        raise TypeError("s must satisfy gcd(s, 3k) = 1 and 3|(k+s)")
    else:
        s_vals = (s,)
        
    def _poly(s_val, u_val, w_val, v_val):
        e1 = (2**s_val + 1)
        e2 = (2**(n - k) + 2**(k + s_val)) % (2**n - 1)
        e3 = (2**(n - k) + 1) % (2**n - 1)
        e4 = (2**s_val + 2**(k + s_val)) % (2**n - 1)
        return (u_val * x**e1 + u_val**(2**k) * x**e2 + F(v_val) * x**e3 + F(w_val) * u_val**(2**k + 1) * x**e4)

    pairs = (
        (_poly(s_val, u_val, w_val, v_val), {'s': s_val, 'u': u_val, 'v': v_val, 'w': w_val})
        for u_val in u_vals
        for s_val in s_vals
        for v_val, w_val in pair_vw
    )
    return aggregate_results(pairs)


def family_11(n, k=None, i=None, a=None):
    r"""
    Return the Budaghyan-Helleseth-Kaleyski construction from 2020.
    Defined by `f(x) = x^3 + a(x^2^i + 1) + bx^(3 * 2^m) + c(x^(2^(i + m) + 2^m))^2^k`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``k`` -- (optional) positive integer; if None, returns a list for all valid k in {0, ... , n - 1}
    - ``i`` -- (optional) positive integer; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``a`` -- (optional) primitive element of GF(2^2); if None, returns a list for all primitive elements of GF(2^2)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_11
        sage: F.<a> = GF(2^10)
        sage: family_11(10, 2, 3, a^5 + a^3 + a)
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3

        sage: family_11(10, 2, None, a^5 + a^3 + a)
        {x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3: 
            [{'k': 2, 'i': 9, 'a': a^5 + a^3 + a}],
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3: 
            [{'k': 2, 'i': 3, 'a': a^5 + a^3 + a}],
        (a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3: 
            [{'k': 2, 'i': 5, 'a': a^5 + a^3 + a}],
        (a^5 + a^3 + a)*x^516 + x^144 + (a^5 + a^3 + a + 1)*x^96 + x^3: 
            [{'k': 2, 'i': 7, 'a': a^5 + a^3 + a}]}

        sage: family_11(10, 2, 5)
        {(a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3: 
            [{'k': 2, 'i': 5, 'a': a2}],
        (a^5 + a^3 + a)*x^132 + (a^5 + a^3 + a)*x^96 + x^3: 
            [{'k': 2, 'i': 5, 'a': a2 + 1}]}

        sage: result = family_11(10); list(result.keys())
        [(a^5 + a^3 + a)*x^513 + (a^5 + a^3 + a + 1)*x^96 + x^48 + x^3
        (a^5 + a^3 + a + 1)*x^513 + (a^5 + a^3 + a)*x^96 + x^48 + x^3
        ...
        (a^5 + a^3 + a + 1)*x^384 + (a^5 + a^3 + a)*x^96 + x^12 + x^3]

        sage: len(result)
        39
    """
    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    if m % 2 == 0 or m % 3 == 0:
        raise TypeError("m must be odd and not divisible by 3")
    
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    K = F.subfield(2)

    if k is None:
        k_vals = range(n)
    elif not (0 <= k < n):
        raise TypeError("k must be in {0, ... , n - 1}")
    else:
        k_vals = (k,)
    
    def _valid_i(i_val, k_val):
        """
        If is k even: i in {m-2, m, n-1, (m-2)^(-1) mod n}, where (m-2)^(-1) exists if gcd(m-2, n) = 1
        If is k odd: i in {m+2, m, (m+2)^(-1) mod n}, where (m+2)^(-1) exists if gcd(m+2, n) = 1
        """
        p = (m - 2) if k_val % 2 == 0 else (m + 2)
        base = {m - 2, m, n - 1} if k_val % 2 == 0 else {m + 2, m}

        # If i_val is given, check that it is in the base set or is the inverse of p mod n
        if i_val is not None:
            if i_val in base or (gcd(p, n) == 1 and (i_val * p) % n == 1):
                return (i_val,)
            return ()
        
        # Compute the inverse of p mod n if it exists and add it to the base set
        inv_set = {pow(p, -1, n)} if gcd(p, n) == 1 else set()
        return tuple(base | inv_set) # Return the union of the base set and the inverse if it exists

    pair_ki = [(k_val, i_val) for k_val in k_vals for i_val in _valid_i(i, k_val)]
    if not pair_ki:
        raise TypeError("No valid pairs of k and i")
    
    if a is None:
        g = K.multiplicative_generator()
        a_vals = (g, g**2)
    elif a == 0 or not is_primitive_element(K, a):
        raise TypeError("a must be a primitive element of GF(2^2)")
    else:
        a_vals = (a,)
    
    def _poly(k_val, i_val, a_val):
        e1 = ((2**i_val + 1) * 2**k_val) % (2**n - 1)
        e2 = 3 * 2**m
        e3 = ((2**(i_val + m) + 2**m) * 2**k_val) % (2**n - 1)
        return x**3 + a_val * x**e1 + a_val**2 * x**e2 + F(1) * x**e3

    pairs = (
        (_poly(k_val, i_val, a_val), {'k': k_val, 'i': i_val, 'a': a_val})
        for k_val, i_val in pair_ki
        for a_val in a_vals
    )
    return aggregate_results(pairs)


def family_12(n, i=None, s=None, a=None, b=None, c=None):
    r"""
    Return the Zheng-Kan-Li-Peng-Tang construction from 2022.
    Defined by `f(x) = a * Tr^n_m(bx^(2^i + 1)) + a^q * Tr^n_m(cx^(2^s + 1))`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``i`` -- (optional) integer with gcd(i, n) = 1; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``s`` -- (optional) integer; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``a`` -- (optional) element of GF(2^n) not in GF(q) with a + a^q != 0; if None, returns a list for all valid a in GF(2^n)
    - ``b`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero b in GF(2^n)
    - ``c`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero c in GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_12
        sage: F.<a> = GF(2^10)
        sage: family_12(10, 1, 5, a^8 + a^5 + a^4 + a^3 + a^2, a^9 + a^8 + a^7 + a^3 + a^2 + a, a^8 + a^7 + a^6 + a^2 + a)
        (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^7 + a^2 + a)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^3
        
        sage: family_12(10, 1, 5, a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1)
        {(a^7 + a^4 + a^3 + a^2 + a)*x^96 + (a^8 + a^6 + a^4 + a^3 + a + 1)*x^33 + (a^9 + a^8 + a^5 + a^4 + a^2 + 1)*x^3: 
            [{'i': 1, 's': 5, 'a': a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1, 'b': a, 'c': a},
            ...
            {'i': 1, 's': 5, 'a': a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1, 'b': a, 'c': a^9 + a^6 + a^4 + a}],
        ...
        (a^8 + a^6 + 1)*x^96 + (a^7 + a^6 + a^4 + a + 1)*x^33 + (a^9 + a^8 + a^7 + a^6 + a^5 + 1)*x^3: 
            [{'i': 1, 's': 5, 'a': a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1, 'b': a^9 + a^5 + a^4 + a^2 + a + 1, 'c': a^7 + a^3 + a^2 + a + 1},
            ...
            {'i': 1, 's': 5, 'a': a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1, 'b': a^9 + a^5 + a^4 + a^2 + a + 1, 'c': a^9 + a^6 + a^5 + a^2}]}

        sage: family_12(10, 1, None, a^8 + a^5 + a^4 + a^3 + a^2,  a^9 + a^8 + a^7 + a^3 + a^2 + a)
        {(a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^9 + a^7 + a^6 + a^4)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^3: 
            [{'i': 1, 's': 5, 'a': a^8 + a^5 + a^4 + a^3 + a^2, 'b': a^9 + a^8 + a^7 + a^3 + a^2 + a, 'c': a},
            ...
            {'i': 1, 's': 5, 'a': a^8 + a^5 + a^4 + a^3 + a^2, 'b': a^9 + a^8 + a^7 + a^3 + a^2 + a, 'c': a^9 + a^6 + a^4 + a}],
        ...
        (a^9 + a^5 + a^4 + a + 1)*x^288 + (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^9 + a^5 + a^4 + a + 1)*x^9 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^3: 
            [{'i': 1, 's': 3, 'a': a^8 + a^5 + a^4 + a^3 + a^2, 'b': a^9 + a^8 + a^7 + a^3 + a^2 + a, 'c': 1}]}

        sage: result = family_12(10, 1, None, a^6 + a^4 + 1); list(result.keys())
        [(a^8 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^9 + a^8 + a^7 + a^4 + a + 1)*x^33 + (a^7 + a^5 + a)*x^3,
        (a^8 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^9 + a^3 + 1)*x^33 + (a^7 + a^5 + a)*x^3,
        ...
        (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^513 + (a^6 + a^4 + 1)*x^96 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^48 + (a^6 + a^4 + 1)*x^3]

        sage: len(result)
        137423

        sage: result = family_12(6); list(result.keys())
        [(a^5 + a^4 + a^2 + 1)*x^24 + (a^5 + a^2 + a)*x^9 + a^2*x^3,
        (a^5 + a^4 + 1)*x^24 + (a^3 + a^2 + 1)*x^9 + a^3*x^3,
        ...
        (a^5 + a^4 + a^2 + a)*x^33 + x^12 + (a^5 + a^3 + a^2)*x^9]

        sage: len(result)
        12432
    """
    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    q = 2**m
    if m % 2 == 0:
        raise TypeError("m must be odd")

    F = GF(2**n, 'a')
    Q = GF(q)
    R = PolynomialRing(F, 'x')
    x = R.gen()

    if a is None:
        a_vals = tuple(a_val for a_val in F if a_val not in Q and a_val + a_val**q != F(0))
    elif a not in F or a in Q or a + a**q == F(0):
        raise TypeError("a must be an element of GF(2^n) not in GF(q) with a + a^q != 0")
    else:
        a_vals = (a,)
    
    if i is None:
        i_vals = tuple(i_val for i_val in range(1, n) if gcd(i_val, n) == 1)
    elif gcd(i, n) != 1:
        raise TypeError("i must satisfy gcd(i, n) = 1")
    else:
        i_vals = (i,)
    
    if b is None:
        b_vals = tuple(b_val for b_val in F if b_val != F(0))
    elif b == 0 or b not in F:
        raise TypeError("b must be a nonzero element of GF(2^n)")
    else:
        b_vals = (b,)
    
    if c is None:
        c_vals = tuple(c_val for c_val in F if c_val != F(0))
    elif c == 0 or c not in F:
        raise TypeError("c must be a nonzero element of GF(2^n)")
    else:
        c_vals = (c,)
    
    pair_ibcs = [(i_val, b_val, c_val, s_val) for i_val in i_vals for b_val in b_vals for c_val in c_vals for s_val in (family12_s_candidates(F, i_val, b_val, c_val) if s is None else ([s] if family12_validates_s(F, i_val, b_val, c_val, s) else []))]
    if not pair_ibcs:
        raise ValueError("No valid combinations of i, b, c, s found")

    def _poly(i_val, s_val, a_val, b_val, c_val):
        e1 = 2**i_val + 1
        e2 = (q * e1) % (2**n - 1)
        e3 = 2**s_val + 1
        e4 = (q * e3) % (2**n - 1)
        return a_val * b_val * x**e1 + a_val * b_val**q * x**e2 + a_val**q * c_val * x**e3 + a_val**q * c_val**q * x**e4

    pairs = (
        (_poly(i_val, s_val, a_val, b_val, c_val), {'i': i_val, 's': s_val, 'a': a_val, 'b': b_val, 'c': c_val})
        for i_val, b_val, c_val, s_val in pair_ibcs
        for a_val in a_vals
    )
    return aggregate_results(pairs)


def family_13(n, s=None, v=None, mu=None):
    r"""
    Return the Li-Zhou-Li-Qu construction from 2022.
    Defined by `f(x) = L(z)^(2^m + 1) + cz^(2^m + 1)`.

    NOTE: When all optional parameters are given, it returns a single polynomial. When one or more are None, returns a dict mapping each valid polynomial found to a list of parameter sets that yield that polynomial.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``s`` -- (optional) integer with gcd(s, n/3) = 1; if None, returns a list for all valid s in {1, ... , m}
    - ``v`` -- (optional) nonzero element of GF(2^(n/3)); if None, returns a list for all nonzero v in GF(2^m)
    - ``mu`` -- (optional) nonzero element of GF(2^n) satisfying mu^(2^(2*(n/3)) + 2^(n/3) + 1) != 1; if None, returns a list for all valid mu in GF(2^n)
    
    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import family_13
        sage: F.<a> = GF(2^9)
        sage: M = F.subfield(3)
        sage: family_13(9, 1, M(1), a^7 + a^5 + a^3 + a + 1)
        x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10
        
        sage: family_13(9, 1, M(1))
        {x^144 + a^5*x^130 + x^129 + (a^5 + a^4 + a^2 + a)*x^32 + x^24 + (a^7 + a^6 + a^5 + a^4 + a + 1)*x^18 + (a^5 + a^4 + a^2 + a)*x^17 + a^5*x^10: 
            [{'s': 1, 'v': 1, 'mu': a^5}],
        ...
        x^144 + (a^6 + a)*x^130 + x^129 + (a^7 + a^5 + a^3 + a + 1)*x^32 + x^24 + (a^7 + a^6 + a^4 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a + 1)*x^17 + (a^6 + a)*x^10: 
            [{'s': 1, 'v': 1, 'mu': a^6 + a}]}
        
        sage: family_13(9, 1, None, a^7 + a^5 + a^3 + a + 1)
        {x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10 + (a^8 + a^6 + a^4 + 1)*x^9: 
            [{'s': 1, 'v': a^3, 'mu': a^7 + a^5 + a^3 + a + 1}],
        ...
        x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10: 
            [{'s': 1, 'v': 1, 'mu': a^7 + a^5 + a^3 + a + 1}]}

        sage: family_13(9, 1)
        {x^144 + a^5*x^130 + x^129 + (a^5 + a^4 + a^2 + a)*x^32 + x^24 + (a^7 + a^6 + a^5 + a^4 + a + 1)*x^18 + (a^5 + a^4 + a^2 + a)*x^17 + a^5*x^10 + (a^8 + a^6 + a^4 + 1)*x^9: 
            [{'s': 1, 'v': 1, 'mu': a^5}],
        ...
        x^144 + (a^6 + a)*x^130 + x^129 + (a^7 + a^5 + a^3 + a + 1)*x^32 + x^24 + (a^7 + a^6 + a^4 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a + 1)*x^17 + (a^6 + a)*x^10: 
            [{'s': 1, 'v': 1, 'mu': a^6 + a}]}

        sage: result = family_13(12); list(result.keys())
        [x^544 + a*x^514 + x^513 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^64 + x^48 + (a^11 + a^10 + a^8 + a^7 + a^3 + a + 1)*x^34 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^33 + a*x^18 + (a^10 + a^9 + a^8 + a^4 + a^3 + a^2 + 1)*x^17,
        x^544 + a*x^514 + x^513 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^64 + x^48 + (a^11 + a^10 + a^8 + a^7 + a^3 + a + 1)*x^34 + (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^33 + a*x^18 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a)*x^17,
        ...
        x^2176 + (a^10 + a^9 + a^8 + a^7 + a^4 + a^2 + a + 1)*x^2056 + x^2049 + (a^11 + a^9 + a^3 + a)*x^256 + x^144 + (a^9 + a^7 + a^5 + a^3 + a^2 + a)*x^136 + (a^11 + a^9 + a^3 + a)*x^129 + (a^10 + a^9 + a^8 + a^7 + a^4 + a^2 + a + 1)*x^24]

        sage: len(result)
        32760
    """        
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")

    m = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    M = F.subfield(m)

    if v is None:
        v_vals = tuple(v_val for v_val in M if v_val != 0)
    elif v == 0 or v not in M:
        raise TypeError("v must be a nonzero element of GF(2^m)")
    else:
        v_vals = (v,)

    if s is None:
        s_vals = tuple(s_val for s_val in range(1, m + 1) if gcd(s_val, m) == 1)    
    elif gcd(s, m) != 1:
        raise TypeError("gcd(s, m) must be 1")
    else:
        s_vals = (s,)
    
    if mu is None:
        mu_vals = tuple(mu_val for mu_val in F if mu_val != 0 and mu_val**(2**(2*m) + 2**m + 1) != 1)
    elif mu not in F or mu == 0 or mu**(2**(2*m) + 2**m + 1) == 1:
        raise TypeError("mu must be a nonzero element of GF(2^n) satisfying mu^(2^(2*m)+2^m+1) != 1")
    else:
        mu_vals = (mu,)
    
    def _permutes(s_val, mu_val):
        """
        Check that L(z) = z^(2^(m+s)) + mu*z^(2^s) + z is a permutation of GF(2^n).
        A linearized polynomial permutes GF(2^n) if it has no nonzero roots, equivalently if gcd(L(x), x^(2^n) - x) = x.
        """
        L = x**(2**(m + s_val)) + mu_val * x**(2**s_val) + x
        return L % (x**(2**n) - x) if L.gcd(x**(2**n) - x) == x else None

    pair_smu = tuple((s_val, mu_val, L) for mu_val in mu_vals for s_val in s_vals if (L := _permutes(s_val, mu_val)) is not None)

    def _poly(L, v_val):
        return L**(2**m + 1) + v_val * x**(2**m + 1)
    
    pairs = (
        (_poly(L, v_val), {'s': s_val, 'v': v_val, 'mu': mu_val})
        for s_val, mu_val, L in pair_smu
        for v_val in v_vals
    )
    return aggregate_results(pairs)
