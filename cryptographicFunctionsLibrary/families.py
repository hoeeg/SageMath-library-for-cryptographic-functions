from sage.all import *
from helpers import is_primitive_element, family12_conditions, family12_check_s


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
        u = [u_val for u_val in F if  u_val != 0 and is_primitive_element(F, u_val)]
    elif not is_primitive_element(F, u):
        raise TypeError("u must be a primitive element of GF(2^n)")
    else:
        u = [u]

    if s is None:
        s = [s_val for s_val in range(1, n) if gcd(s_val, 3 * k) == 1]
    elif gcd(s, 3*k) != 1:
        raise TypeError("gcd(s, 3*k) must be 1")
    else:
        s = [s]

    def _poly(s_val, u_val):
        i = (s_val * k) % p
        m = p - i
        e_x = (2**(i * k) + 2**(m * k + s_val)) % (2**n - 1)
        return (x**(2**s_val + 1) + u_val**(2**k - 1) * x**e_x).mod(x**(2**n) - x)
    
    res = set()
    for u_val in u:
        for s_val in s:    
            res.add(_poly(s_val, u_val))
    
    if not res:
        raise TypeError("No valid polynomials found")
    
    return list(res) if len(res) > 1 else list(res)[0]


def family_1(n, s=None, u=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2008 for `p = 3`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and gcd(k, 3) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_1
        sage: F.<a> = GF(2^12)
        sage: family_1(12, 5, a^2)
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33

        sage: family_1(12, None, a^6 + a^5 + 1)
        [x^2049 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^264,
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33,
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^528 + x^3,
        x^129 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^24]

        sage: family_1(12, 1)
        [(a^9 + a^8 + a^7 + a^5 + 1)*x^528 + x^3,
        (a^11 + a^9 + a^7 + a^4 + a^3 + a^2 + a)*x^528 + x^3,
        ...
        (a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^528 + x^3]

        sage: result = family_1(12); result
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

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 4k
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and gcd(k, 3) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_2
        sage: F.<a> = GF(2^16)
        sage: family_2(16, 5, a^15 + a^14 + a^13 + a^11 + a^3 + a^2 + a)
        (a^15 + a^14 + a^13 + a^12 + a^10 + a^9 + a^8 + a^4 + a^2 + a + 1)*x^33

        sage: family_2(16, None, a^15 + a^14 + a^7 + a^6 + a^3 + a)
        [(a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^129,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^8193,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^2049,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^33]
 
        sage: family_2(16, 1)
        [(a^15 + a^12 + a^11 + a^8 + a^7 + a^6 + a)*x^3,
        (a^15 + a^14 + a^13 + a^12 + a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^3,
        ...
        (a^15 + a^13 + a^10 + a^8 + a^6 + a^5 + a^3 + a^2)*x^3]
        
        sage: result = family_2(16); result
        [(a^15 + a^12 + a^11 + a^8 + a^7 + a^6 + a)*x^3,
        (a^14 + a^12 + a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^2049,
        ...
        (a^7 + a^6 + a^5 + a^3 + a^2)*x^129]
        
        sage: len(result)
        20480
    """
    return _family_1_2(n, 4, s, u)


def family_3(n, i=None, s=None, c=None):
    r"""
    Return the Budaghyan-Carlet construction from 2008.
    Defined by `f(x) = sx^(q + 1) + x^(2^i + 1) + x^(q * (2^i + 1)) + cx^(2^i * q + 1) + c^q * x^(2^i + q))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be even 
    - ``i`` -- (optional) positive integer with gcd(i, n/2) = 1; if None, returns a list over all valid i in {1, ... , m - 1}
    - ``s`` -- (optional) element of GF(2^n) not in GF(q); if None, returns a list over all valid s in GF(2^n) \ GF(q)
    - ``c`` -- (optional) element of GF(2^n); if None, returns a list over all valid c in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_3
        sage: F.<a> = GF(2^6)
        sage: family_3(6, 1, a, a)
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3

        sage: family_3(6, None, a^2 + a + 1, a^5 + a^4 + a^2 + a + 1)
        [x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5,
        x^24 + (a^5 + a^4 + a^2 + a + 1)*x^17 + a*x^10 + (a^2 + a + 1)*x^9 + x^3]
        
        sage: family_3(6, 1)
        [x^24 + (a^4 + a^3 + a^2)*x^17 + (a^3 + a)*x^10 + (a^3 + 1)*x^9 + x^3,
        x^24 + (a^5 + a^4 + a^3 + 1)*x^17 + (a^3 + a^2 + 1)*x^10 + (a^4 + a^3 + a^2 + a)*x^9 + x^3,
        ...
        x^24 + a^4*x^17 + (a^5 + a^4 + a)*x^10 + (a^5 + a^4 + a)*x^9 + x^3]

        sage: result = family_3(6); result
        [x^24 + (a^4 + a^3 + a^2)*x^17 + (a^3 + a)*x^10 + (a^5 + a^4 + a^2)*x^9 + x^3,
        x^24 + (a^3 + a)*x^17 + (a^4 + a^3 + a^2)*x^10 + (a^5 + a^4 + a^2)*x^9 + x^3,
        ...
        x^40 + (a^4 + a^2 + 1)*x^33 + a^5*x^12 + (a^5 + a^4 + a^3)*x^9 + x^5]

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

    if c is None:
        c = [c_val for c_val in F]
    elif c not in F:
        raise TypeError("c must be an element of GF(2^n)")
    else:
        c = [c]
    
    if s is None:
        s = [s_val for s_val in F if s_val not in K]
    elif s in K or s not in F:
        raise TypeError("s must be in GF(2^n) but not in GF(q)")
    else:
        s = [s]
    
    if i is None:
        i = [i_val for i_val in range(1, m) if gcd(i_val, m) == 1]
    elif gcd(i, m) != 1:
        raise TypeError("gcd(i, m) must be 1")
    else:
        i = [i]
    
    def _poly(i_val, s_val, c_val):
        P = x**(2**i_val + 1) + c_val * x**(2**i_val) + c_val**q * x + 1
        K_gen = F.gen()**(q - 1)
        v = F(1)
        for _ in range(q + 1):
            if P(v) == 0:
                raise TypeError("The polynomial x^{2^i+1} + cx^{2^i} + c^q x + 1 has a root satisfying x^{q+1} = 1")
            v *= K_gen

        e_xq = (q * (2**i_val + 1)) % (2**n - 1)
        e_c = (2**i_val * q + 1) % (2**n - 1)
        return (s_val * x**(q + 1) + x**(2**i_val + 1) + x**e_xq + c_val * x**e_c + c_val**q * x**(2**i_val + q)).mod(x**(2**n) - x)
    
    res = set()
    for i_val in i:
        for s_val in s:
            for c_val in c:
                try:
                    res.add(_poly(i_val, s_val, c_val))
                except TypeError:
                    continue

    if not res:
        raise TypeError("No valid polynomials found")
    
    return list(res) if len(res) > 1 else list(res)[0]


def family_4(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr_n(a^3 * x^9)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_4
        sage: F.<a> = GF(2^9)
        sage: family_4(9, a^8 + a^7 + a^5 + a^3 + 1)
        (a^8 + a^7 + a^6 + a^4 + a^3)*x^288 + (a^8 + a^7 + a^5 + a^3 + a^2 + a + 1)*x^260 + (a^8 + a^7 + a^5 + 1)*x^144 + (a^7 + a^6 + 1)*x^130 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^72 + (a^6 + a^5 + a^4 + a^2 + a)*x^65 + (a^7 + a^6 + a^4 + a^3 + a)*x^36 + (a^8 + a^6 + a^4 + a^3 + 1)*x^18 + (a^5 + a^4 + a^2 + 1)*x^9 + x^3

        sage: family_4(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3

        sage: result = family4(9); result
        [(a^6 + a^3 + 1)*x^288 + (a^7 + a^5 + a^2 + 1)*x^260 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^8 + a^5 + a^4)*x^72 + (a^8 + a^7 + a^4 + a)*x^65 + (a^6 + a^2)*x^36 + a^5*x^18 + a^2*x^9 + x^3,
        (a^7 + a^6 + a^3 + 1)*x^288 + a*x^260 + (a^8 + a^5 + a^3 + a^2)*x^144 + (a^8 + a^6 + a^4 + a^3 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^2 + a)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + 1)*x^65 + (a^7 + a^4 + a^3)*x^36 + (a^5 + a)*x^18 + a^4*x^9 + x^3,
        ...
        x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3]

        sage: len(result)
        511
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
  
    def _poly(a_val):
        trace = sum(a_val**(3 * 2**i) * x**((9 * 2**i) % (2**(n) - 1)) for i in range(n))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)
    
    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    
    return _poly(a)


def family_5(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^3 * x^9 + a^6 * x^18)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family5
        sage: F.<a> = GF(2^9)
        sage: family_5(9, a^6 + a^5)
        (a^3 + a^2 + 1)*x^144 + (a^8 + a^5 + a^3 + a^2)*x^130 + (a^4 + a^3 + a + 1)*x^72 + (a^3 + a + 1)*x^65 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a)*x^9 + x^3
        
        sage: family_5(9, F(1))
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3

        sage: result = family_5(9); result
        [(a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^8 + a^5 + a^4)*x^72 + (a^8 + a^7 + a^4 + a)*x^65 + a^5*x^18 + a^2*x^9 + x^3,
        (a^8 + a^5 + a^3 + a^2)*x^144 + (a^8 + a^6 + a^4 + a^3 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^2 + a)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + 1)*x^65 + (a^5 + a)*x^18 + a^4*x^9 + x^3,        ...
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

    def _poly(a_val):
        trace = sum(a_val**(3 * 2**(3 * i)) * x**((9 * 2**(3 * i)) % (2**(n) - 1))+ a_val**(6 * 2**(3 * i)) * x**((18 * 2**(3 * i) % (2**(n) - 1))) for i in range(k))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)

    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    
    return _poly(a)


def family_6(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^6 * x^18 + a^12 * x^36)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_6
        sage: F.<a> = GF(2^9)
        sage: family_6(9, a^8 + a^5 + a^3 + 1)
        (a^8 + a^7 + a^6 + a^4 + a^3 + a)*x^288 + (a^8 + a^7 + a^6 + a^2 + 1)*x^260 + (a^7 + a^6 + a^5 + a^4)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^3 + 1)*x^36 + (a^8 + a^7 + a^5 + a^4 + a^2 + a + 1)*x^18 + x^3

        sage: family_6(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3

        sage: result = family_6(9); result
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

    def _poly(a_val):
        trace = sum(a_val**(6 * 2**(3 * i)) * x**((18 * 2**(3 * i)) % (2**(n) - 1))+ a_val**(12 * 2**(3 * i)) * x**((36 * 2**(3 * i)) % (2**(n) - 1)) for i in range(k))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)

    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    elif a == 0 or a not in F:
        raise TypeError("a must be a nonzero element of GF(2^n)")
    
    return _poly(a)


def family_7_9(n, s=None, u=None, v=None, w=None):
    r"""
    Return the Bracken-Byrne-Markin-McGuire construction from 2011.
    Defined by `f(x) = ux^(2^s + 1) + u^(2^k) * x^(2^-k + 2^(k + s)) + vx^(2^-k + 1) + wu^(2^k + 1) * x^(2^s + 2^(k + s))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k, gcd(k, 3) = 1
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and 3|(k + s); if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)
    - ``v`` -- (optional) element of GF(2^k) with v*w != 1; if None, returns a list for all valid v in GF(2^k)
    - ``w`` -- (optional) element of GF(2^k) with v*w != 1; if None, returns a list for all valid w in GF(2^k)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_7_9
        sage: F.<a> = GF(2^12)
        sage: u, v, w = a, a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1, a^10 + a^9 + a^8 + a^4 + a^3 + a^2
        sage: family_7_9(12, 5, u, v, w)
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^11 + a^10 + a^9 + a^6 + a^4 + a^2 + a)*x^544 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257 + a*x^33
        
        sage: v, w = a^11 + a^9 + a^8 + a^6 + a^3 + a + 1, a^10 + a^9 + a^6 + a^5 + a^3 + 1
        sage: family_7_9(12, 11, None, v, w)
        [(a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a)*x^2056 + (a^10 + a^9 + a^8 + a^6 + a^2 + 1)*x^2049 + (a^10 + a^8 + a^5 + a^4 + a^3 + a + 1)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        (a^10 + a^9 + a^7 + a^6 + a^4 + a^3)*x^2056 + (a^11 + a^10 + a^8 + a^6 + a^4 + a^3 + a^2 + a)*x^2049 + (a^11 + a^7 + a^6 + a^3)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        ...
        (a^7 + a^5 + a^4 + a^3 + a)*x^2056 + (a^11 + a^9 + a^6 + a^5 + a^4 + a + 1)*x^2049 + (a^9 + a^8 + a^7 + a^6 + a^2)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257]
        
        sage: family_7_9(12, 11, None, v)
        [(a^11 + a^6 + a^4 + a^3 + a^2 + a)*x^2056 + (a^11 + a^9 + a^7 + a^5 + a^3 + a^2 + 1)*x^2049 + (a^11 + a^10 + a^8 + a^6 + a^4 + a^3)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        (a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^4 + 1)*x^2056 + (a^11 + a^10 + a^9 + a^8 + a^6 + a^4 + a^3 + a^2 + a)*x^2049 + (a^11 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        ...
        (a^7 + a^5 + a^4 + a)*x^2056 + (a^11 + a^8 + a^3 + 1)*x^2049 + (a^11 + a^10 + a^7 + a^4 + 1)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257]

        sage: family_7_9(12, 11, None, None, w)
        [(a^10 + a^9 + a^8 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^2056 + (a^8 + a^7 + a^5 + a)*x^2049 + (a^10 + a^7 + a + 1)*x^264 + (a^10 + a^9 + a^8 + a^4 + a^3 + a^2 + 1)*x^257,
        (a^10 + a^9 + a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^2056 + (a^8 + a^7 + a^3 + 1)*x^2049 + (a^11 + a^10 + a^6 + a^4 + a + 1)*x^264 + (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^257,
        ...
        (a^11 + a^10 + a^8 + a^6 + a^5 + a^3)*x^2056 + (a^11 + a^8 + a^7 + a^6 + a^2)*x^2049 + (a^11 + a^10 + a^8 + a^7 + a^5 + a^2 + a + 1)*x^264 + (a^10 + a^9 + a^8 + a^4 + a^3 + a^2)*x^257]

        sage: family_7_9(12, 5, a)
        [(a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^11 + a^10 + a^8 + a^7 + a^3 + a + 1)*x^544 + (a^8 + a^6 + a^5 + a^4 + a^2 + 1)*x^257 + a*x^33,
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^9 + a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^544 + (a^11 + a^10 + a^8 + a^5 + a + 1)*x^257 + a*x^33
        ...
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^11 + a^10 + a^6 + a^4 + a + 1)*x^544 + a*x^33]

        sage: result = family_7_9(12, 5); result
        [(a^11 + a^9 + a^8 + a^7 + a^5 + a^2 + a + 1)*x^768 + (a^7 + a^3 + a)*x^544 + (a^8 + a^6 + a^5 + a^4 + a^2)*x^257 + (a^10 + a^9 + a^8 + a^7 + a^4 + a^3 + a + 1)*x^33,
        (a^11 + a^10 + a^3 + a^2)*x^768 + (a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^3 + a + 1)*x^544 + x^257 + (a^7 + a^4 + a^3 + a)*x^33,
        ...
        (a^11 + a^10 + a^9 + a^8 + a^4 + a^2 + a)*x^768 + (a^11 + a^10 + a^8 + a^7 + a^6 + a^4 + a^3 + a)*x^544 + (a^11 + a^10 + a^8 + a^5 + a)*x^257 + (a^9 + a^7 + a^6 + a)*x^33]
        
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
        u = [u_val for u_val in F if  u_val != 0 and is_primitive_element(F, u_val)]
    elif not is_primitive_element(F, u) or u == 0:
        raise TypeError("u must be a primitive element of GF(2^n)")
    else:
        u = [u]

    if v is not None and v not in K:
        raise TypeError("v must be an element of GF(2^k)")
    if w is not None and w not in K:
        raise TypeError("w must be an element of GF(2^k)")

    if v is not None and w is not None:
        if v * w == K(1):
            raise TypeError("v and w must be elements of GF(2^k)")
        pair = {(v, w)}
    else:
        v_vals = K if v is None else [v]
        w_vals = K if w is None else [w]

        pair = set((v_val, w_val) for v_val in v_vals for w_val in w_vals if v_val * w_val != K(1))

    if s is None:
        s = [s_val for s_val in range(1, n) if (gcd(s_val, 3*k) == 1 and (k + s_val) % 3 == 0)]
    elif gcd(s, 3*k) != 1 or (k + s) % 3 != 0:
        raise TypeError("s must satisfy gcd(s, 3k) = 1 and 3|(k+s)")
    else:
        s = [s]
        
    def _poly(s_val, w_val, v_val, u_val):
        e_ux = (2**(n - k) + 2**(k + s_val)) % (2**n - 1)
        e_v = (2**(n - k) + 1) % (2**n - 1)
        e_wu = (2**s_val + 2**(k + s_val)) % (2**n - 1)
        return (u_val * x**(2**s_val + 1) + u_val**(2**k) * x**e_ux + v_val * x**e_v + w_val * u_val**(2**k + 1) * x**e_wu).mod(x**(2**n) - x)

    res = set()
    for u_val in u:
        for v_val, w_val in pair:
            for s_val in s:
                res.add(_poly(s_val, w_val, v_val, u_val))

    if not res:
        raise TypeError("No valid polynomials found")

    return list(res) if len(res) > 1 else list(res)[0]


def family_11(n, k=None, i=None, a=None):
    r"""
    Return the Budaghyan-Helleseth-Kaleyski construction from 2020.
    Defined by `f(x) = x^3 + a(x^2^i + 1) + bx^(3 * 2^m) + c(x^(2^(i + m) + 2^m))^2^k`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``k`` -- (optional) positive integer; if None, returns a list for all valid k in {0, ... , n - 1}
    - ``i`` -- (optional) positive integer; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``a`` -- (optional) primitive element of GF(2^2); if None, returns a list for all primitive elements of GF(2^2)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_11
        sage: F.<a> = GF(2^10)
        sage: family_11(10, 2, 3, a^5 + a^3 + a)
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3

        sage: family_11(10, 2, None, a^5 + a^3 + a)
        [(a^5 + a^3 + a)*x^516 + x^144 + (a^5 + a^3 + a + 1)*x^96 + x^3,
        x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3,
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3,
        (a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3]

        sage: family_11(10, 2, 5)
        [(a^5 + a^3 + a)*x^132 + (a^5 + a^3 + a)*x^96 + x^3,
        (a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3]

        sage: result = family_11(10); result
        [(a^5 + a^3 + a + 1)*x^576 + (a^5 + a^3 + a)*x^96 + x^18 + x^3,
        (a^5 + a^3 + a)*x^516 + x^144 + (a^5 + a^3 + a + 1)*x^96 + x^3,
        ...
        x^576 + (a^5 + a^3 + a)*x^96 + (a^5 + a^3 + a + 1)*x^18 + x^3]

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
        k = list(range(0, n))
    else:
        k = [k]
    
    def _valid_i(i_arg, k_val):
        p = (m - 2) if k_val % 2 == 0 else (m + 2)
        base = {m - 2, m, n - 1} if k_val % 2 == 0 else {m + 2, m}

        if i_arg is not None:
            if i_arg not in base and (i_arg * p) % n != 1:
                raise TypeError("i is not valid")
            return [i_arg]
        
        inv_set = {pow(p, -1, n)} if gcd(p, n) == 1 else set()
        return base | inv_set

    i = {k_val: _valid_i(i, k_val) for k_val in k}
    
    if a is None:
        a = [a_val for a_val in K if a_val != 0 and is_primitive_element(K, a_val)]
    elif a == 0 or not is_primitive_element(K, a) or a not in K:
        raise TypeError("a must be a primitive element of GF(2)")
    else:
        a = [a]
    
    def _poly(k_val, i_val, a_val, b_val, c_val):
        e_a = (2**i_val + 1) * 2**k_val % (2**n - 1)
        e_b = 3 * 2**m
        e_c = ((2**(i_val + m) + 2**m) * 2**k_val) % (2**n - 1)
        return (x**3 + a_val*x**e_a + b_val*x**e_b + c_val*x**e_c).mod(x**(2**n) - x)

    res = set()
    for a_val in a:
        b_val, c_val = a_val**2, F(1)
        for k_val in k:
            for i_val in i[k_val]:
                res.add(_poly(k_val, i_val, a_val, b_val, c_val))

    if not res:
        raise TypeError("No valid polynomials found")

    return list(res) if len(res) > 1 else list(res)[0]


def family_12(n, i=None, s=None, a=None, b=None, c=None):
    r"""
    Return the Zheng-Kan-Li-Peng-Tang construction from 2022.
    Defined by `f(x) = a * Tr^n_m(bx^(2^i + 1)) + a^q * Tr^n_m(cx^(2^s + 1))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``i`` -- (optional) integer with gcd(i, n) = 1; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``s`` -- (optional) integer; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``a`` -- (optional) element of GF(2^n) not in GF(q) with a + a^q != 0; if None, returns a list for all valid a in GF(2^n)
    - ``b`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero b in GF(2^n)
    - ``c`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero c in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family12
        sage: F.<a> = GF(2^10)
        sage: family_12(10, 1, 5, a^8 + a^5 + a^4 + a^3 + a^2,  a^9 + a^8 + a^7 + a^3 + a^2 + a, a^8 + a^7 + a^6 + a^2 + a)
        (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^7 + a^2 + a)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^3

        sage: family_12(10, 1, 5, a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1)
        [(a^6 + a^3 + a^2 + a + 1)*x^96 + (a^9 + a^8 + a^7 + a^3 + a + 1)*x^33 + (a^9 + a^7 + a)*x^3,
        (a^7 + a^6 + a^5 + 1)*x^96 + (a^9 + a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^33 + (a^9 + a^4 + a^3 + a + 1)*x^3,
        ...
        (a^9 + a^8 + a^6 + a^3 + a^2 + 1)*x^96 + (a^9 + a^6 + a^5 + a^2 + a)*x^33 + (a^9 + a^5 + a^3 + a + 1)*x^3]

        sage: family_12(10, None, None, a^8 + a^5 + a^4 + a^3 + a^2,  a^9 + a^8 + a^7 + a^3 + a^2 + a, a^8 + a^7 + a^6 + a^2 + a)
        [(a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^129 + (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^36 + (a^7 + a^2 + a)*x^33,
        (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^96 + (a^7 + a^2 + a)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^3,
        (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^288 + (a^7 + a^2 + a)*x^33 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^9,
        (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^513 + (a^7 + a^6 + a^5 + a^4 + a^3 + a)*x^48 + (a^7 + a^2 + a)*x^33]

        sage: result = family_12(10, 1, None, a^9 + a^8 + a^7 + a^5 + a^3 + a^2 + 1); result
        [(a^9 + a^8 + a^7 + a^3 + a^2)*x^96 + (a^9 + a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^33 + (a^8 + a^5 + a^3 + a^2)*x^3,
        (a^7 + a^3)*x^129 + (a^6 + a^2 + a + 1)*x^96 + (a^8 + a^7 + a^5 + a^4 + a^2 + a)*x^36 + (a^9 + a^7 + a^5 + a^4 + a^3 + a + 1)*x^3,
        ...
        (a^9 + a^8 + a^5 + a^4 + a^3 + a + 1)*x^513 + (a^7 + a^3 + a^2 + a)*x^96 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^3)*x^48 + (a^8 + a^7 + a^2 + 1)*x^3]

        sage: len(result)
        137423
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
        a = [a_val for a_val in F if a_val not in Q and a_val + a_val**q != F(0)]
    elif a not in F or a in Q:
        raise TypeError("a must be an element of GF(2^n) not in GF(q) with a + a^q != 0")
    else:
        a = [a]
    
    if i is None:
        i = [i_val for i_val in range(1, n) if gcd(i_val, n) == 1]
    elif gcd(i, n) != 1:
        raise TypeError("i must satisfy gcd(i, n) = 1")
    else:
        i = [i]
    
    b = [b] if b is not None else [b_val for b_val in F if b_val != F(0)]
    c = [c] if c is not None else [c_val for c_val in F if c_val != F(0)]

    def _poly(i_val, a_val, b_val, c_val, s_val):
        e_x1 = (2**i_val + 1)
        e_b = (q * e_x1) % (2**n - 1)
        e_x2 = (2**s_val + 1)
        e_c = (q * e_x2) % (2**n - 1)
        return (a_val * b_val * x**e_x1 + a_val * b_val**q * x**e_b + a_val**q * c_val * x**e_x2 + a_val**q * c_val**q * x**e_c).mod(x**(2**n) - x)
    
    res = set()
    for a_val in a:
        for i_val in i:
            for b_val in b:
                for c_val in c:
                    if s is not None:
                        if family12_check_s(F, i_val, b_val, c_val, s):
                            res.add(_poly(i_val, a_val, b_val, c_val, s))
                    else:
                        for s_val in family12_conditions(F, i_val, b_val, c_val):
                            res.add(_poly(i_val, a_val, b_val, c_val, s_val))
                    
    if not res:
        raise TypeError("No valid polynomials found")
    
    return list(res) if len(res) > 1 else list(res)[0]


def family_13(n, s=None, v=None, mu=None):
    r"""
    Return the Li-Zhou-Li-Qu construction from 2022.
    Defined by `f(x) = L(z)^(2^m + 1) + cz^(2^m + 1)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``s`` -- (optional) integer with gcd(s, n/3) = 1; if None, returns a list for all valid s in {1, ... , m}
    - ``v`` -- (optional) nonzero element of GF(2^(n/3)); if None, returns a list for all nonzero v in GF(2^m)
    - ``mu`` -- (optional) nonzero element of GF(2^n) satisfying mu^(2^(2*(n/3)) + 2^(n/3) + 1) != 1; if None, returns a list for all valid mu in GF(2^n)
    
    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family_13
        sage: F.<a> = GF(2^9)
        sage: Fm = F.subfield(3)
        sage: family_13(9, 1, Fm(1), a^7 + a^5 + a^3 + a + 1)
        x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10

        sage: family_13(9, 1, Fm(1))
        [x^144 + (a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^130 + x^129 + (a^4 + a)*x^32 + x^24 + (a^7 + a^6 + a^5 + a^3 + 1)*x^18 + (a^4 + a)*x^17 + (a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^10,
        x^144 + (a^8 + a^6 + a^5 + a^4 + a^2)*x^130 + x^129 + (a^8 + a^7 + a^5 + a)*x^32 + x^24 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^18 + (a^8 + a^7 + a^5 + a)*x^17 + (a^8 + a^6 + a^5 + a^4 + a^2)*x^10,
        ...
        x^144 + (a^7 + a^6 + a^5 + a^3 + a)*x^130 + x^129 + (a^3 + a^2 + a + 1)*x^32 + x^24 + (a^8 + a^7 + a^6 + a^5 + a^2)*x^18 + (a^3 + a^2 + a + 1)*x^17 + (a^7 + a^6 + a^5 + a^3 + a)*x^10]

        sage: family_13(9, 1, None, a^7 + a^5 + a^3 + a + 1)
        [x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10 + (a^8 + a^6 + a^3 + a^2)*x^9,
        x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10 + (a^4 + a^3 + a^2 + 1)*x^9,
        ...        
        x^144 + (a^7 + a^5 + a^3 + a + 1)*x^130 + x^129 + (a^8 + a^7 + a^5 + a^2 + 1)*x^32 + x^24 + (a^7 + a^6 + a^5)*x^18 + (a^8 + a^7 + a^5 + a^2 + 1)*x^17 + (a^7 + a^5 + a^3 + a + 1)*x^10]

        sage: family_13(9, 1)
        [x^144 + (a^3 + a^2)*x^130 + x^129 + (a^7 + a^5 + a^4 + a^2 + 1)*x^32 + x^24 + (a^8 + a^6 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^4 + a^2 + 1)*x^17 + (a^3 + a^2)*x^10 + (a^8 + a^6 + a^3 + a^2 + 1)*x^9,
        x^144 + (a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^130 + x^129 + (a^6 + a^4 + a^3 + a^2 + 1)*x^32 + x^24 + (a^6 + a^5 + a^4 + a^2 + 1)*x^18 + (a^6 + a^4 + a^3 + a^2 + 1)*x^17 + (a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^10 + (a^4 + a^3 + a^2)*x^9,
        ...
        x^144 + (a^7 + a^5 + 1)*x^130 + x^129 + (a^7 + a^6 + a^4 + a^2)*x^32 + x^24 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^18 + (a^7 + a^6 + a^4 + a^2)*x^17 + (a^7 + a^5 + 1)*x^10 + (a^8 + a^6 + a^3 + a^2)*x^9]

        sage: family_13(9)
        [x^144 + (a^3 + a^2)*x^130 + x^129 + (a^7 + a^5 + a^4 + a^2 + 1)*x^32 + x^24 + (a^8 + a^6 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^4 + a^2 + 1)*x^17 + (a^3 + a^2)*x^10 + (a^8 + a^6 + a^3 + a^2 + 1)*x^9,
        x^288 + (a^8 + a^3 + a^2 + a + 1)*x^260 + x^257 + (a^8 + a^7 + a^6 + a^5 + a^2 + a)*x^64 + x^40 + (a^8 + a^6 + a^5 + a + 1)*x^36 + (a^8 + a^7 + a^6 + a^5 + a^2 + a)*x^33 + (a^8 + a^3 + a^2 + a + 1)*x^12 + (a^8 + a^6 + a^4 + 1)*x^9,
        ...
        x^288 + (a^7 + a^5 + a^4 + a^3 + a^2 + 1)*x^260 + x^257 + (a^7 + a^6 + a^3)*x^64 + x^40 + (a^8 + a^7 + a^5 + a + 1)*x^36 + (a^7 + a^6 + a^3)*x^33 + (a^7 + a^5 + a^4 + a^3 + a^2 + 1)*x^12 + (a^8 + a^6 + a^4 + 1)*x^9]

        sage: result = family_13(12); result
        [x^544 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^514 + x^513 + (a^7 + a^6 + a^4 + a^3 + a)*x^64 + x^48 + (a^11 + a^9 + a^5)*x^34 + (a^7 + a^6 + a^4 + a^3 + a)*x^33 + (a^10 + a^8 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^18 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^17,
        x^2176 + (a^10 + a^9 + a^5 + a^4 + a^3 + a^2)*x^2056 + x^2049 + (a^10 + a^8 + a^7 + a^6 + a^5 + a + 1)*x^256 + x^144 + (a^10 + a^8 + a^7 + a^5 + a^4 + a^3 + 1)*x^136 + (a^10 + a^8 + a^7 + a^6 + a^5 + a + 1)*x^129 + (a^10 + a^9 + a^5 + a^4 + a^3 + a^2)*x^24 + (a^8 + a^6 + a^5 + a^4 + a^2)*x^17,
        ...
        x^544 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^514 + x^513 + (a^9 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^64 + x^48 + (a^11 + a^8 + a^5 + a^3 + a^2)*x^34 + (a^9 + a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^33 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^18 + (a^10 + a^9 + a^8 + a^4 + a^3 + a^2)*x^17]

        sage: len(result)
        32760
    """
    def _permutes(s_val, mu_val):
        L = x**(2**(m + s_val)) + mu_val*x**(2**s_val) + x
        return L.gcd(x**(2**n) - x) == x
        
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")

    m = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    Fm = F.subfield(m)

    if v is None:
        v = [v_val for v_val in Fm if v_val != 0]
    elif v == 0 or v not in Fm:
        raise TypeError("v must be a nonzero element of GF(2^m)")
    else:
        v = [v]

    if s is not None:
        if gcd(s, m) != 1:
            raise TypeError("gcd(s, m) must be 1")
    
    s = [s_val for s_val in range(1, m + 1) if gcd(s_val, m) == 1] if s is None else [s]
    
    if mu is None:
        mu = [mu_val for mu_val in F if mu_val != 0 and mu_val**(2**(2*m) + 2**m + 1) != 1]
    elif mu not in F or mu == 0 or mu**(2**(2*m) + 2**m + 1) == 1:
        raise TypeError("mu must be a nonzero element of GF(2^n) satisfying mu^(2^(2*m)+2^m+1) != 1")
    else:
        mu = [mu]

    def _poly(s_val, mu_val, v_val):
        L = x**(2**(m + s_val)) + mu_val*x**(2**s_val) + x
        return (L**(2**m + 1) + v_val*x**(2**m + 1)).mod(x**(2**n) - x)

    res = set()

    for mu_val in mu:
        for s_val in s:
            if not _permutes(s_val, mu_val):
                continue
            for v_val in v:
                res.add(_poly(s_val, mu_val, v_val))
    
    if not res:
        raise TypeError("No valid polynomials found")
    
    return list(res) if len(res) > 1 else list(res)[0]
