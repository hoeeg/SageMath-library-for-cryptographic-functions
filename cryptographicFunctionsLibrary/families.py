from sage.all import *
from helpers import is_primitive_element


def _family1_2(n, p, s, u):
    """
    Shared implementation for family1 (p=3) and family2 (p=4).
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


def family1(n, s=None, u=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2008 for `p = 3`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k, gcd(k, 3) = 1
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family1
        sage: F.<a> = GF(2^12)
        sage: family1(12, 5, a^2)
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33

        sage: family1(12, None, a^6 + a^5 + 1)
        [x^2049 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^264,
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^768 + x^33,
        (a^11 + a^9 + a^7 + a^6 + a^2)*x^528 + x^3,
        x^129 + (a^11 + a^9 + a^7 + a^6 + a^2)*x^24]

        sage: family1(12, 1)
        [(a^9 + a^8 + a^7 + a^5 + 1)*x^528 + x^3,
        (a^11 + a^9 + a^7 + a^4 + a^3 + a^2 + a)*x^528 + x^3,
        ...
        (a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^528 + x^3]

        sage: result = family1(12); result
        [(a^9 + a^6 + a^5 + a^4 + a^2 + 1)*x^768 + x^33,
        (a^11 + a^10 + a^8 + a^7 + a^6 + a^3 + a^2 + a)*x^768 + x^33,
        ...
        x^129 + (a^11 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^24]

        sage: len(result)
        576
    """
    return _family1_2(n, 3, s, u)


def family2(n, s=None, u=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2008 for `p = 4`.
    Defined by `f(x) = x^(2^s + 1) + u^(2^k - 1) * x^(2^(ik) + 2^(mk + s))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k, gcd(k, 3) = 1
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1; if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n); if None, returns a list for all primitive elements of GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family2
        sage: F.<a> = GF(2^16)
        sage: family2(16, 5, a^15 + a^14 + a^13 + a^11 + a^3 + a^2 + a)
        (a^15 + a^14 + a^13 + a^12 + a^10 + a^9 + a^8 + a^4 + a^2 + a + 1)*x^33

        sage: family2(16, None, a^15 + a^14 + a^7 + a^6 + a^3 + a)
        [(a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^3,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^129,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^8193,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^2049,
        (a^9 + a^8 + a^7 + a^3 + a^2 + a)*x^33]
 
        sage: family2(16, 1)
        [(a^15 + a^12 + a^11 + a^8 + a^7 + a^6 + a)*x^3,
        (a^15 + a^14 + a^13 + a^12 + a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^3,
        ...
        (a^15 + a^13 + a^10 + a^8 + a^6 + a^5 + a^3 + a^2)*x^3]
        
        sage: result = family2(16); result
        [(a^15 + a^12 + a^11 + a^8 + a^7 + a^6 + a)*x^3,
        (a^14 + a^12 + a^11 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1)*x^2049,
        ...
        (a^7 + a^6 + a^5 + a^3 + a^2)*x^129]
        
        sage: len(result)
        20480
    """
    return _family1_2(n, 4, s, u)


def family3(n, i=None, s=None, c=None):
    r"""
    Return the Budaghyan-Carlet construction from 2008.
    Defined by `f(x) = sx^(q + 1) + x^(2^i + 1) + x^(q * (2^i + 1)) + cx^(2^i * q + 1) + c^q * x^(2^i + q))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be even 
    - ``i`` -- (optional) positive integer with gcd(i, n/2) = 1; if None, returns a list over all valid i in {1, ... , m - 1}
    - ``s`` -- (optional) element of GF(2^n) not in GF(q); if None, returns a list over all valid s in {1, ... , n - 1}
    - ``c`` -- (optional) element of GF(2^n); if None, returns a list over all valid c in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family3
        sage: F.<a> = GF(2^6)
        sage: family3(6, 1, a, a)
        x^24 + a*x^17 + (a^5 + a^4 + a^2 + a + 1)*x^10 + a*x^9 + x^3

        sage: family3(6, None, a^2 + a + 1, a^5 + a^4 + a^2 + a + 1)
        [x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5,
        x^24 + (a^5 + a^4 + a^2 + a + 1)*x^17 + a*x^10 + (a^2 + a + 1)*x^9 + x^3]
        
        sage: family3(6, 1)
        [x^24 + (a^4 + a^3 + a^2)*x^17 + (a^3 + a)*x^10 + (a^3 + 1)*x^9 + x^3,
        x^24 + (a^5 + a^4 + a^3 + 1)*x^17 + (a^3 + a^2 + 1)*x^10 + (a^4 + a^3 + a^2 + a)*x^9 + x^3,
        ...
        x^24 + a^4*x^17 + (a^5 + a^4 + a)*x^10 + (a^5 + a^4 + a)*x^9 + x^3]

        sage: result = family3(6); result
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


def family4(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr_n(a^3 * x^9)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family4
        sage: F.<a> = GF(2^9)
        sage: a = a^8 + a^7 + a^5 + a^3 + 1
        sage: family4(9, a)
        (a^7 + a^5 + a^4 + 1)*x^288 + (a^7 + a^6 + a^4 + a^2)*x^260 + (a^8 + a^7 + a^5 + a^4 + a + 1)*x^144 + (a^8 + a^7 + a^5 + 1)*x^130 + (a^8 + a^6 + a + 1)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^65 + (a^7 + a^5 + a^2 + 1)*x^36 + (a^8 + a^6 + a^4 + a^2 + a)*x^18 + (a^7 + a^4 + a^2 + a)*x^9 + x^3

        sage: F.<a> = GF(2^9)
        sage: family4(9, F(1))
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


def family5(n, a=None):
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
        sage: family5(9, a^6 + a^5)
        (a^3 + a^2 + 1)*x^144 + (a^8 + a^5 + a^3 + a^2)*x^130 + (a^4 + a^3 + a + 1)*x^72 + (a^3 + a + 1)*x^65 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a)*x^9 + x^3
        
        sage: family5(9, F(1))
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3

        sage: result = family5(9); result
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


def family6(n, a=None):
    r"""
    Return the Budaghyan-Carlet-Leander construction from 2009.
    Defined by `f(x) = x^3 + a^-1 * Tr^n_3(a^6 * x^18 + a^12 * x^36)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n); if None, returns a list over all valid a in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family6
        sage: F.<a> = GF(2^9)
        sage: family6(9, a^8 + a^5 + a^3 + 1)
        (a^8 + a^7 + a^6 + a^4 + a^3 + a)*x^288 + (a^8 + a^7 + a^6 + a^2 + 1)*x^260 + (a^7 + a^6 + a^5 + a^4)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^3 + 1)*x^36 + (a^8 + a^7 + a^5 + a^4 + a^2 + a + 1)*x^18 + x^3

        sage: family6(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3

        sage: result = family6(9); result
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


def family7_9(n, s=None, u=None, v=None, w=None):
    r"""
    Return the Bracken-Byrne-Markin-McGuire construction from 2011.
    Defined by `f(x) = ux^(2^s + 1) + u^(2^k) * x^(2^-k + 2^(k + s)) + vx^(2^-k + 1) + wu^(2^k + 1) * x^(2^s + 2^(k + s))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must satisfy n = 3k, gcd(k, 3) = 1
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1 and 3|(k + s); if None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)
    - ``v`` -- (optional) element of GF(2^n) with v*w != 1; if None, returns a list for all valid v in GF(2^n)
    - ``w`` -- (optional) element of GF(2^n) with v*w != 1; if None, returns a list for all valid w in GF(2^n)

    NOTE: v and w must be supplied together or not at all

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family7_9
        sage: F.<a> = GF(2^12)
        sage: u, v, w = a, a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1, a^10 + a^9 + a^8 + a^4 + a^3 + a^2
        sage: family7_9(12, 5, u, v, w)
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^11 + a^10 + a^9 + a^6 + a^4 + a^2 + a)*x^544 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257 + a*x^33
        
        sage: family7_9(12, 5, a)
        [(a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^9 + a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^544 + (a^11 + a^9 + a^8 + a^6 + a^3 + a)*x^257 + a*x^33,
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^8 + a^7 + a^6 + a^4 + a^3)*x^544 + x^257 + a*x^33,
        ...
        (a^11 + a^10 + a^9 + a^7 + a^5 + a^4)*x^768 + (a^11 + a^10 + a^9 + a^6 + a^4 + a^2 + a)*x^544 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^257 + a*x^33]

        sage: v, w = a^11 + a^9 + a^8 + a^6 + a^3 + a + 1, a^10 + a^9 + a^6 + a^5 + a^3 + 1
        sage: family7_9(12, 11, None, v, w)
        [(a^10 + a^9 + a^7 + a^6 + a^4 + a^3)*x^2056 + (a^11 + a^10 + a^8 + a^6 + a^4 + a^3 + a^2 + a)*x^2049 + (a^11 + a^7 + a^6 + a^3)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        (a^11 + a^10 + a^7 + a^6 + a^2 + 1)*x^2056 + (a^10 + a^9 + a^7 + a^5 + a^4 + a^3 + 1)*x^2049 + (a^10 + a^9 + a^8 + a^7 + a^6 + a^5 + a^3)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257,
        ...
        (a^11 + a^10 + a^9 + a^8 + a^7 + a^4 + a^3 + a)*x^2056 + (a^11 + a^8 + a^4 + a + 1)*x^2049 + (a^11 + a^6 + a^4 + 1)*x^264 + (a^11 + a^9 + a^8 + a^6 + a^3 + a + 1)*x^257]

        sage: result = family7_9(12, 5); result
        [(a^6 + a^4 + a)*x^768 + (a^9 + a^7 + a^5 + a^4 + a^3 + a^2 + a)*x^544 + (a^8 + a^6 + a^5 + a^4 + a^2 + 1)*x^257 + (a^11 + a^9 + a^7 + a^3 + a^2 + 1)*x^33,
        (a^10 + a^6 + a^5 + a^3 + a)*x^768 + (a^9 + a^8 + a^6 + a^5 + a^2)*x^544 + (a^10 + a^9 + a^6 + a^5 + a^3 + 1)*x^257 + (a^7 + a^5 + a^3 + a^2 + 1)*x^33,
        ...
        (a^9 + a^7 + a^5 + a^3 + a + 1)*x^768 + (a^9 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^544 + (a^11 + a^9 + a^5 + a^4 + a^3 + a^2 + a)*x^257 + (a^9 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^33]
        
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

    # TODO: can have just one of them
    if (v is None) != (w is None):
        raise TypeError("Supply both v and w, or neither")
    
    if v is None and w is None:
        pair = set((v_val, w_val) for v_val in K for w_val in K if v_val * w_val != K(1))
    elif v not in K or w not in K:
        raise TypeError("v and w must be elements of GF(2^k)")
    elif v*w == K(1):
        raise TypeError("vw must not be 1")
    else:
        pair = {(v, w)}
    
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
     

def family10():
    pass


def family11(n, k=None, i=None, a=None):
    r"""
    Return the Budaghyan-Helleseth-Kaleyski construction from 2020.
    Defined by `f(x) = x^3 + a(x^2^i + 1) + bx^(3 * 2^m) + c(x^(2^(i + m) + 2^m))^2^k`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``k`` -- (optional) positive integer; if None, returns a list for all valid k in {0, ... , n - 1}
    - ``ì`` -- (optional) positive integer; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``a`` -- (optional) primitive element of GF(2^2); if None, returns a list for all primitive elements of GF(2^2)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family11
        sage: F.<a> = GF(2^10)
        sage: K = F.subfield(2)
        sage: a = K.gen(); a
        a2
        sage: family11(10, 2, 3, a)
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3

        sage: family11(10, 2, None, a)
        [(a^5 + a^3 + a)*x^516 + x^144 + (a^5 + a^3 + a + 1)*x^96 + x^3,
        x^192 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^6 + x^3,
        x^129 + (a^5 + a^3 + a + 1)*x^96 + (a^5 + a^3 + a)*x^36 + x^3,
        (a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3]

        sage: family11(10, 2, 5)
        [(a^5 + a^3 + a)*x^132 + (a^5 + a^3 + a)*x^96 + x^3,
        (a^5 + a^3 + a + 1)*x^132 + (a^5 + a^3 + a + 1)*x^96 + x^3]

        sage: result = family11(10); result
        [(a^5 + a^3 + a + 1)*x^576 + (a^5 + a^3 + a)*x^96 + x^18 + x^3,
        (a^5 + a^3 + a)*x^516 + x^144 + (a^5 + a^3 + a + 1)*x^96 + x^3,
        ...
        x^576 + (a^5 + a^3 + a)*x^96 + (a^5 + a^3 + a + 1)*x^18 + x^3]
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
    
    if a is None:
        a = [a_val for a_val in K if a_val != 0 and is_primitive_element(K, a_val)]
    elif a == 0 or not is_primitive_element(K, a) or a not in K:
        raise TypeError("a must be a primitive element of GF(2)")
    else:
        a = [a]
   
    def find_i(i_arg, k_val):
        if i_arg is not None:
            if k_val % 2 == 0:
                if i_arg not in {m-2,m,n-1} and (i_arg * (m - 2)) % n != 1:
                    raise TypeError("i is not valid")
            else:
                if i_arg not in {m+2, m} and (i_arg * (m + 2)) % n != 1:
                    raise TypeError("i is not valid")
            return [i_arg]
        
        if k_val % 2 == 0:
            return list({m-2, m, n-1} | {j for j in range(1, n) if (j*(m-2)) % n == 1})
        else:
            return list({m+2, m} | {j for j in range(1, n) if (j*(m+2)) % n == 1})

    def _poly(k_val, i_val, a_val, b_val, c_val):
        e_a = (2**i_val + 1) * 2**k_val % (2**n - 1)
        e_b = 3 * 2**m
        e_c = ((2**(i_val + m) + 2**m) * 2**k_val) % (2**n - 1)
        return (x**3 + a_val*x**e_a + b_val*x**e_b + c_val*x**e_c).mod(x**(2**n) - x)

    res = set()
    for a_val in a:
        b_val, c_val = a_val**2, K(1)
        for k_val in k:
            for i_val in find_i(i, k_val):
                res.add(_poly(k_val, i_val, a_val, b_val, c_val))

    if not res:
        raise TypeError("No valid polynomials found")

    return list(res) if len(res) > 1 else list(res)[0]


def _is_cube(x, F):
    """
    Check if an element is a cube in GF(2^n)
    """
    if x == F(0): 
        return True
    order = F.order() - 1
    return x**(order // 3) == F(1)


def family12(n, i=None, a=None, b=None, c=None):
    r"""
    Return the Zheng-Kan-Li-Peng-Tang construction from 2022.
    Defined by `f(x) = a * Tr^n_m(bx^(2^i + 1)) + a^q * Tr^n_m(cx^(2^s + 1))`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``i`` -- (optional) integer with gcd(i, n) = 1; if None, returns a list for all valid i in {1, ... , n - 1}
    - ``a`` -- (optional) element of GF(2^n) not in GF(q) with a + a^q != 0; if None, returns a list for all valid a in GF(2^n)
    - ``b`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero b in GF(2^n)
    - ``c`` -- (optional) nonzero element of GF(2^n); if None, returns a list for all nonzero c in GF(2^n)

    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family12
        sage: F = GF(2^10)
        sage: a = F.primitive_element()
        sage: family12(10, 1, a, a^2, a^2)
        [(z10^8 + z10^7 + z10^2 + z10 + 1)*x^96 + (z10^9 + z10^5 + z10^4 + z10^3 + z10)*x^33 + z10^3*x^3]

        sage: family12(10, 3, a, a^2, a^2)
        [(z10^8 + z10^7 + z10^2 + z10 + 1)*x^288 + (z10^9 + z10^5 + z10^4 + z10^3 + z10)*x^33 + z10^3*x^9]
        
        sage: family12(10, 1)
        [(z10^9 + z10^7 + z10^6 + z10^5 + z10^4 + z10^2)*x^1056 + (z10^8 + z10^6 + z10^2)*x^96 + (z10^8 + z10^6 + z10^2)*x^33 + z10^2*x^3]
        ...
    """
    def _conditions(i_val, b_val, c_val):
        s_list = []
 
        if not _is_cube(b_val, F):
            term = b_val**(2**(2 * i_val) - 2**i_val + 1) * c_val**(-1)
            if term**q == term:
                s_list.append(3 * i_val)
            
            term = b_val**(2**i_val - 1) * c_val**(-(2**(2 * i_val)))
            if term**q == term and (m - 2 * i_val) > 0:
                s_list.append(m - 2 * i_val)

            term = c_val * b_val**(2**i_val - 1)
            if term**q == term:
                s_list.append(m + 2 * i_val)

            if i_val == 1:
                for s in (1, n):
                    if (s * (m - 2)) % n == 1:
                        term = b_val**(2**(2 * s)) * c_val**(-(2**s - 1))
                        if term**q == term:
                            s_list.append(s)
            
            if c_val**q != c_val:
                s_list.append(m)
            
        term = b_val * c_val**(-(2**i_val))
        if term**q == term:
            s_list.append(n - i_val)

        return list(set(s_list))

    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    q = 2**m
    if m % 2 == 0:
        raise TypeError("m must be odd")

    F = GF(2**n)
    Q = GF(q)
    R = PolynomialRing(F, 'x')
    x = R.gen()

    def _poly(i_val, a_val, b_val, c_val, s_val):
        bx = b_val * x**(2**i_val + 1)
        cx = c_val * x**(2**s_val + 1)
        poly = a_val*(bx + bx**q) + a_val**q * (cx + cx**q)
        return poly.mod(x**(2**n) - x)
    
    i_list = [i] if i is not None else [i for i in range(1, n) if gcd(i, n) == 1]
    a_list = [F(a)] if a is not None else [a for a in F if a not in Q and a + a**q != F(0)]
    b_list = [F(b)] if b is not None else [b for b in F if b != F(0)]
    c_list = [F(c)] if c is not None else [c for c in F if c != F(0)]

    res = []
    for i_val in i_list:
        for a_val in a_list:
            for b_val in b_list:
                for c_val in c_list:
                   for s_val in _conditions(i_val, b_val, c_val):
                        res.append(_poly(i_val, a_val, b_val, c_val, s_val))
                    
    if not res:
        raise TypeError("No valid polynomials found")
    
    return res


def family13(n, s=None, v=None, mu=None):
    r"""
    Return the Li-Zhou-Li-Qu construction from 2022.
    Defined by `f(x) = L(z)^(2^m + 1) + cz^(2^m + 1)`.

    NOTE: When all optional parameters are given, returns a single polynomial. When one or more are None, returns a list of all valid polynomials found.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``s`` -- (optional) integer with gcd(s, n/3) = 1; if None, returns a list for all valid s in {1, ... , m}
    - ``v`` -- (optional) nonzero element of GF(2^(n/3)); randomised if None
    - ``mu`` -- (optional) nonzero element of GF(2^n) satisfying mu^(2^(2*(n/3)) + 2^(n/3) + 1) != 1; randomised if None
    
    EXAMPLES::

        sage: from cryptographicFunctionsLibrary import family13
        sage: F = GF(2^9)
        sage: a = F.gen()
        sage: Fm = GF(2^3)
        sage: family13(9, 1, Fm(1), a^7 + a^5 + a^3 + a + 1)
        x^144 + (z9^7 + z9^5 + z9^3 + z9 + 1)*x^130 + x^129 + (z9^8 + z9^7 + z9^5 + z9^2 + 1)*x^32 + x^24 + (z9^7 + z9^6 + z9^5)*x^18 + (z9^8 + z9^7 + z9^5 + z9^2 + 1)*x^17 + (z9^7 + z9^5 + z9^3 + z9 + 1)*x^10

        sage: Fm = GF(2^3)
        sage: family13(9, 1, Fm(1))
        x^144 + (z9^8 + z9^7 + z9^6 + z9^5 + z9^2 + 1)*x^130 + x^129 + (z9^8 + z9^7 + z9^6 + z9^5 + z9^3 + 1)*x^32 + x^24 + (z9^7 + z9^5 + z9^4 + z9^2 + z9)*x^18 + (z9^8 + z9^7 + z9^6 + z9^5 + z9^3 + 1)*x^17 + (z9^8 + z9^7 + z9^6 + z9^5 + z9^2 + 1)*x^10

        sage: family13(9, 1, None, a^7 + a^5 + a^3 + a + 1)
        x^144 + (z9^7 + z9^5 + z9^3 + z9 + 1)*x^130 + x^129 + (z9^8 + z9^7 + z9^5 + z9^2 + 1)*x^32 + x^24 + (z9^7 + z9^6 + z9^5)*x^18 + (z9^8 + z9^7 + z9^5 + z9^2 + 1)*x^17 + (z9^7 + z9^5 + z9^3 + z9 + 1)*x^10 + (z9^8 + z9^6 + z9^4 + 1)*x^9

        sage: family13(9, 1)
        x^144 + (z9^8 + z9^7 + z9^2 + 1)*x^130 + x^129 + (z9^6 + z9^5 + z9^4 + z9^2)*x^32 + x^24 + (z9^8 + z9^7 + z9^6 + z9^5 + z9)*x^18 + (z9^6 + z9^5 + z9^4 + z9^2)*x^17 + (z9^8 + z9^7 + z9^2 + 1)*x^10 + (z9^8 + z9^6 + z9^4)*x^9

        sage: family13(9)
        [x^144 + (z9^6 + z9^5 + 1)*x^130 + x^129 + (z9^8 + z9^7 + z9^4 + z9^3 + z9^2)*x^32 + x^24 + (z9^8 + z9^7 + z9^2 + z9 + 1)*x^18 + (z9^8 + z9^7 + z9^4 + z9^3 + z9^2)*x^17 + (z9^6 + z9^5 + 1)*x^10 + (z9^8 + z9^6 + z9^4 + 1)*x^9,
        x^288 + (z9^6 + z9^5 + 1)*x^260 + x^257 + (z9^8 + z9^7 + z9^4 + z9^3 + z9^2)*x^64 + x^40 + (z9^8 + z9^7 + z9^2 + z9 + 1)*x^36 + (z9^8 + z9^7 + z9^4 + z9^3 + z9^2)*x^33 + (z9^6 + z9^5 + 1)*x^12 + (z9^8 + z9^6 + z9^4 + 1)*x^9]

        sage: family13(12)
        [x^544 + (z12^10 + z12^9 + z12^8 + z12^7 + z12^5 + z12^3 + z12^2 + z12 + 1)*x^514 + x^513 + (z12^11 + z12^10 + z12^8 + z12^6 + z12^5 + z12^4)*x^64 + x^48 + (z12^10 + z12^9 + z12^5 + z12^4 + z12)*x^34 + (z12^11 + z12^10 + z12^8 + z12^6 + z12^5 + z12^4)*x^33 + (z12^10 + z12^9 + z12^8 + z12^7 + z12^5 + z12^3 + z12^2 + z12 + 1)*x^18 + (z12^10 + z12^9 + z12^6 + z12^5 + z12^3 + 1)*x^17,
        x^2176 + (z12^10 + z12^9 + z12^8 + z12^7 + z12^5 + z12^3 + z12^2 + z12 + 1)*x^2056 + x^2049 + (z12^11 + z12^10 + z12^8 + z12^6 + z12^5 + z12^4)*x^256 + x^144 + (z12^10 + z12^9 + z12^5 + z12^4 + z12)*x^136 + (z12^11 + z12^10 + z12^8 + z12^6 + z12^5 + z12^4)*x^129 + (z12^10 + z12^9 + z12^8 + z12^7 + z12^5 + z12^3 + z12^2 + z12 + 1)*x^24 + (z12^10 + z12^9 + z12^6 + z12^5 + z12^3 + 1)*x^17]
    """
    def _permutes(s_val, mu_val):
        L = x**(2**(m + s_val)) + mu_val*x**(2**s_val) + x
        return all(L(a) != 0 for a in F if a != 0)
        
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")

    m = n // 3
    F = GF(2**n)
    R = PolynomialRing(F, 'x')
    x = R.gen()
    Fm = GF(2**m)

    if v is None:
        while True:
            v = Fm.random_element() 
            if v != 0:
                break
    elif v == 0 or v**(2**m - 1) != 1:
        raise TypeError("v must be a nonzero element of GF(2^m)")

    if s is not None:
        if gcd(s, m) != 1:
            raise TypeError("gcd(s, m) must be 1")
    
    s = [s_val for s_val in range(1, m + 1) if gcd(s_val, m) == 1] if s is None else [s]
    
    if mu is None:
        while True:
            mu = F.random_element()
            if mu == 0 or mu**(2**(2*m) + 2**m + 1) == 1:
                continue
            if all(_permutes(s_val, mu) for s_val in s):
                break
    else:
        if mu not in F or mu == 0:
            raise TypeError("mu must be a nonzero element of GF(2^n)")
        if mu**(2**(2*m) + 2**m + 1) == 1:
            raise TypeError("mu must satisfy mu^(2^(2*m)+2^m+1) != 1")
        for s_val in s:
            if not _permutes(s_val, mu):
                raise TypeError("L does not permute GF(2^n) for s and the given mu")

    def _poly(s_val):
        L = x**(2**(m + s_val)) + mu*x**(2**s_val) + x
        return (L**(2**m + 1) + v*x**(2**m + 1)).mod(x**(2**n) - x)

    if len(s) > 1:
        return [_poly(s_val) for s_val in s]
    
    return _poly(s[0])
