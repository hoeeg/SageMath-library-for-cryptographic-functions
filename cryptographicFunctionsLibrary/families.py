from sage.all import *
from sage.crypto.sbox import SBox


def _family1_2(n, p, s, u):
    """
    Shared implementation for family1 (p=3) and family2 (p=4)

    INPUT:
    - ``n`` -- the degree of the field GF(2^n)
    - ``p`` -- an integer in {3, 4}
    - ``s`` -- an integer coprime to n and 3
    - ``u`` -- a primitive element of GF(2^n)
    """
    if n < 12:
        raise TypeError("n must be >= 12")
    if n % p != 0:
        raise TypeError(f"n must be divisible by p = {p}")
    
    k = n // p
    if gcd(k, 3) != 1:
        raise TypeError("gcd(k, 3) must be 1")

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    if u is None:
        u = F.primitive_element()
    
    def _poly(s_val):
        if gcd(s_val, 3*k) != 1:
            raise TypeError("gcd(s, 3*k) must be 1")
    
        i = (s_val * k) % p
        m = p - i
        e = (2**(i * k) + 2**(m * k + s_val)) % (2**n - 1)
        return x**(2**s_val + 1) + u**(2**k - 1) * x**e
    
    if s is None:
        return [_poly(s_val) for s_val in range(1, n) if gcd(s_val, 3 * k) == 1]
    
    return _poly(s)


def family1(n, s=None, u=None):
    """
    Return the Budaghyan-Carlet-Leander construction from 2008 for p = 3.
    If ``s`` is None every valid s in {1,...,n-1} is tried and the list of resulting polynomials is returned. Otherwise the single polynomial for the given s is returned.

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must satisfy n = 3k, gcd(k,3)=1, n >= 12
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1;
               If None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family1
        sage: family1(12, 1)
        (a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^528 + x^3

        sage: family1(12)
        [(a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^528 + x^3,
        (a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^768 + x^33,
        x^129 + (a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^24,
        x^2049 + (a^10 + a^9 + a^8 + a^6 + a^4 + a^3)*x^264]
    """
    return _family1_2(n, 3, s, u)


def family2(n, s=None, u=None):
    """
    Return the Budaghyan-Carlet-Leander construction from 2008 for p = 4
    If ``s`` is None every valid s in {1,...,n-1} is tried and the list of resulting polynomials is returned. Otherwise the single polynomial for the given s is returned

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must satisfy n = 3k, gcd(k,3)=1, n >= 12
    - ``s`` -- (optional) integer with gcd(s, 3k) = 1;
               If None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family2
        sage: family2(16, 1)
        (a^15 + 1)*x^3
        
        sage: family2(16)
        [(a^15 + 1)*x^3,
        (a^15 + 1)*x^33,
        (a^15 + 1)*x^129,
        (a^15 + 1)*x^2049,
        (a^15 + 1)*x^8193]
    """
    return _family1_2(n, 4, s, u)


def family3(n, i=None, s=None, c=None):
    """
    Return the Budaghyan-Carlet construction from 2008

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must be even 
    - ``i`` -- (optional) positive integer with gcd(i, n/2) = 1;
               If None, returns a list over all valid i in {1, ... , m - 1}
    - ``s`` -- (optional) element of GF(2^n) not in GF(q);
               Randomised if None
    - ``c`` -- (optional) element of GF(2^n);
               Randomised if None

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family3
        sage: family3(6, 1)
        x^24 + (a^4 + a^3 + a^2)*x^17 + (a^3 + a)*x^10 + (a^5 + a^4 + a)*x^9 + x^3

        sage: family3(6)
        [x^24 + (a^5 + a^4 + a)*x^17 + a^4*x^10 + (a^5 + a^4 + a^2)*x^9 + x^3,
        x^40 + (a^5 + a^4 + a)*x^33 + a^4*x^12 + (a^5 + a^4 + a^2)*x^9 + x^5]

        sage: F.<a> = GF(2^6)
        sage: c = a^5 + a^4 + a^2 + a + 1
        sage: s = a^2 + a + 1
        sage: family3(6, 1, s, c)
        x^24 + (a^5 + a^4 + a^2 + a + 1)*x^17 + a*x^10 + (a^2 + a + 1)*x^9 + x^3

        sage: family3(6, None, s, c)
        [x^24 + (a^5 + a^4 + a^2 + a + 1)*x^17 + a*x^10 + (a^2 + a + 1)*x^9 + x^3,
        x^40 + (a^5 + a^4 + a^2 + a + 1)*x^33 + a*x^12 + (a^2 + a + 1)*x^9 + x^5]
    """
    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()
    K = GF(2**m)

    if c is None:
        c = F.random_element()

    if s is None:
        while True:
            candidate = F.random_element()
            if candidate not in K:
                s = candidate
                break
    elif s in K:
        raise TypeError("s must not be in GF(2^m)")
    
    def _poly(i_val):
        if gcd(i_val, m) != 1:
            raise TypeError("gcd(i, m) must be 1")

        P = x**(2**i_val + 1) + c * x**(2**i_val) + c**q * x + 1
        K_gen = F.gen()**(q - 1)
        v = F(1)
        for _ in range(q + 1):
            if P(v) == 0:
                raise ValueError("The polynomial x^{2^i+1} + cx^{2^i} + c^q x + 1 has a root satisfying x^{q+1} = 1")
            v *= K_gen
    
        return (s * x**(q + 1) + x**(2**i_val + 1) + x**(q * (2**i_val + 1)) + c * x**(2**i_val * q + 1) + c**q * x**(2**i_val + q))
    
    if i is None:
        return [_poly(i_val) for i_val in range(1, m) if gcd(i_val, m) == 1]
    
    return _poly(i)


def family4(n, a=None):
    """
    Return the Budaghyan-Carlet-Leander construction from 2009

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``a`` -- (optional) nonzero element of GF(2^n);
               If None, returns a list over all valid a in GF(2^n)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family4
        sage: family4(9)
        [(a^6 + a^3 + 1)*x^288 + (a^7 + a^5 + a^2 + 1)*x^260 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^8 + a^5 + a^4)*x^72 + (a^8 + a^7 + a^4 + a)*x^65 + (a^6 + a^2)*x^36 + a^5*x^18 + a^2*x^9 + x^3,
        ...
        x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3]
        
        sage: F.<a> = GF(2^9)
        sage: family4(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^72 + x^65 + x^36 + x^18 + x^9 + x^3

        sage: F.<a> = GF(2^9)
        sage: a = a^8 + a^7 + a^5 + a^3 + 1
        sage: family4(9, a)
        (a^7 + a^5 + a^4 + 1)*x^288 + (a^7 + a^6 + a^4 + a^2)*x^260 + (a^8 + a^7 + a^5 + a^4 + a + 1)*x^144 + (a^8 + a^7 + a^5 + 1)*x^130 + (a^8 + a^6 + a + 1)*x^72 + (a^8 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^65 + (a^7 + a^5 + a^2 + 1)*x^36 + (a^8 + a^6 + a^4 + a^2 + a)*x^18 + (a^7 + a^4 + a^2 + a)*x^9 + x^3
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    
    def _poly(a_val):
        trace = sum((a_val**(3 * 2**i) * x**(9 * 2**i)) for i in range(n))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)
    
    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    
    if a == 0:
        raise TypeError("a must be nonzero")
    if a not in F:
        raise TypeError("a must be an element of GF(2^n)")
    
    return _poly(a)



def family5(n, a=None):
    """
    Return the Budaghyan-Carlet-Leander construction from 2009

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n);
               If None, returns a list over all valid a in GF(2^n)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family5
        sage: family5(9)
        [(a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^8 + a^5 + a^4)*x^72 + (a^8 + a^7 + a^4 + a)*x^65 + a^5*x^18 + a^2*x^9 + x^3,
        ...
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3]
        
        sage: F.<a> = GF(2^9)
        sage: family5(9, F(1))
        x^144 + x^130 + x^72 + x^65 + x^18 + x^9 + x^3

        sage: F.<a> = GF(2^9)
        sage: a = a^6 + a^5
        sage: family5(9, a)
        (a^3 + a^2 + 1)*x^144 + (a^8 + a^5 + a^3 + a^2)*x^130 + (a^4 + a^3 + a + 1)*x^72 + (a^3 + a + 1)*x^65 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + a + 1)*x^18 + (a^7 + a^5 + a^3 + a)*x^9 + x^3
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
    
    k = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    def _poly(a_val):
        trace = sum(a_val**(3 * 2**(3 * i)) * x**(9 * 2**(3 * i)) + a_val**(6 * 2**(3 * i)) * x**(18 * 2**(3 * i)) for i in range(k))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)

    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    
    if a == 0:
        raise TypeError("a must be nonzero")
    if a not in F:
        raise TypeError("a must be an element of GF(2^n)")
    
    return _poly(a)


def family6(n, a=None):
    """
    Return the Budaghyan-Carlet-Leander construction from 2009

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must be divisible by 3
    - ``a`` -- (optional) nonzero element of GF(2^n);
               If None, returns a list over all valid a in GF(2^n)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family6
        sage: family6(9)
        [(a^6 + a^3 + 1)*x^288 + (a^7 + a^5 + a^2 + 1)*x^260 + (a^8 + a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^4 + a + 1)*x^130 + (a^6 + a^2)*x^36 + a^5*x^18 + x^3,
        ...
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3]

        sage: F.<a> = GF(2^9)
        sage: family6(9, F(1))
        x^288 + x^260 + x^144 + x^130 + x^36 + x^18 + x^3

        sage: F.<a> = GF(2^9)
        sage: a = a^8 + a^5 + a^3 + 1
        sage: family6(9, a)
        (a^8 + a^7 + a^6 + a^4 + a^3 + a)*x^288 + (a^8 + a^7 + a^6 + a^2 + 1)*x^260 + (a^7 + a^6 + a^5 + a^4)*x^144 + (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^130 + (a^8 + a^7 + a^6 + a^5 + a^3 + 1)*x^36 + (a^8 + a^7 + a^5 + a^4 + a^2 + a + 1)*x^18 + x^3
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
    
    k = n // 3
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    def _poly(a_val):
        trace = sum(a_val**(6 * 2**(3 * i)) * x**(18 * 2**(3 * i)) + a_val**(12 * 2**(3 * i)) * x**(36 * 2**(3 * i)) for i in range(k))
        return (x**3 + (F(1) / a_val) * trace).mod(x**(2**n) - x)

    if a is None:
        return [_poly(a_val) for a_val in F if a_val != 0]
    
    if a == 0:
        raise TypeError("a must be nonzero")
    if a not in F:
        raise TypeError("a must be an element of GF(2^n)")
    
    return _poly(a)


def family7(n, s=None, u=None, v=None, w=None):
    """
    Return the Bracken-Byrne-Markin-McGuire construction from 2011

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n);
               Must satisfy n = 3k, gcd(k,3) = 1
    - ``s`` -- (optional) integer with gcd(s,3k)=1 and 3|(k+s);
               If None, returns a list for all valid s in {1, ... , n - 1}
    - ``u`` -- (optional) primitive element of GF(2^n)
    - ``v`` -- (optional) element of GF(2^n) with v*w != 1
               Randomised if None
    - ``w`` -- (optional) element of GF(2^n) with v*w != 1
               Randomised if None

    NOTE: v and w must be supplied together or not at all

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family7
        sage: family7(12, 5)
        (z12^11 + z12^10 + z12^9 + z12^7 + z12^5 + z12^4)*x^768 + (z12^11 + z12^9 + z12^8 + z12^7 + z12 + 1)*x^544 + (z12^11 + z12^10 + z12^6 + z12^4 + z12^2 + z12 + 1)*x^257 + z12*x^33

        sage: family7(12)
        [(z12^11 + z12^10 + z12^9 + z12^7 + z12^5 + z12^4)*x^768 + (z12^11 + z12^10 + z12^8 + z12^5 + z12 + 1)*x^257 + z12*x^33,
        z12*x^2049 + (z12^11 + z12^10 + z12^9 + z12^7 + z12^5 + z12^4)*x^264 + (z12^11 + z12^10 + z12^8 + z12^5 + z12 + 1)*x^257]
    """
    if n % 3 != 0:
        raise TypeError("n must be divisible by 3")
 
    k = n // 3
    K = GF(2**k)

    if gcd(k, 3) != 1:
        raise TypeError("k must be coprime to 3")

    F = GF(2**n)
    R = PolynomialRing(F, 'x')
    x = R.gen()

    if u is None:
        u = F.primitive_element()
    
    if (v is None) != (w is None):
        raise TypeError("Supply both v and w, or neither")
    
    if v is None and w is None:
        while True:
            v = K.random_element()
            w = K.random_element()
            if v*w != K(1):
                break
        v, w = F(v), F(w)
    else:
        if v*w == K(1):
            raise TypeError("vw must not be 1")
    
    def _poly(s_val):
        if gcd(s_val, 3*k) != 1 or (k + s_val) % 3 != 0:
            raise TypeError("s must satisfy gcd(s, 3k) = 1 and 3|(k+s)")
        
        return (u * x**(2**s_val + 1) + u**(2**k) * x**(2**(n - k) + 2**(k + s_val)) + v * x**(2**(n - k) + 1) + w * u**(2**k + 1) * x**(2**s_val + 2**(k + s_val))).mod(x**(2**n) - x)

    if s is None:
        return [_poly(s_val) for s_val in range(1, n) if gcd(s_val, 3 * k) == 1 and (k + s_val) % 3 == 0]
    
    return _poly(s)
     

def family10():
    pass


def family11(n, k, i=None, a=None):
    """
    Return the Budaghyan-Helleseth-Kaleyski construction from 2020

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``k`` -- positive integer with gcd(k, n) = 1
    - ``ì`` -- (optional) positive integer
    - ``a`` -- (optional) primitive element of GF(2^2)

    EXAMPLE USAGE::
        sage: from cryptographicFunctionsLibrary import family11
        sage: family11(10, 2, 3)
        (z10^5 + z10^3 + z10)*x^36 + x^3

        sage: family11(10, 2)
        [x^192 + (z10^5 + z10^3 + z10 + 1)*x^96 + (z10^5 + z10^3 + z10)*x^6 + x^3,
        x^129 + (z10^5 + z10^3 + z10 + 1)*x^96 + (z10^5 + z10^3 + z10)*x^36 + x^3,
        (z10^5 + z10^3 + z10 + 1)*x^132 + (z10^5 + z10^3 + z10 + 1)*x^96 + x^3,
        (z10^5 + z10^3 + z10)*x^516 + x^144 + (z10^5 + z10^3 + z10 + 1)*x^96 + x^3]
    """
    if n % 2 != 0:
        raise TypeError("n must be even")
    
    m = n // 2
    if m % 2 == 0 or m % 3 == 0:
        raise TypeError("m must be odd and not divisible by 3")
    
    F = GF(2**n)
    R = PolynomialRing(F, 'x')
    x = R.gen()
    K = GF(2**2)

    if a is None:
        a = K.primitive_element()
   
    if n == 10 and i == 3 and k == 2:
        b = K(0)
        c = K(0)
    else:
        b = a**2
        c = K(1)

    def f(i):
        poly = x**3 + F(a)*(x**(2**i+1))**(2**k) + F(b)*x**(3*2**m) + F(c)*(x**(2**(i+m)+2**m))**(2**k)
        return poly.mod(x**(2**n) - x)
    
    if i is not None:
        if k % 2 == 0:
            if i not in {m-2,m,n-1} and (i * (m - 2)) % n != 1:
                raise TypeError("i is not valid")
        else:
            if i not in {m+2, m} and (i * (m + 2)) % n != 1:
                raise TypeError("i is not valid")
    
        return f(i)
    
    if k % 2 == 0:
        candidates = {m-2, m, n-1} | {i for i in range(1, n) if (i * (m - 2)) % n == 1}
    else:
        candidates = {m+2, m} | {i for i in range(1, n) if (i * (m + 2)) % n == 1}
    
    return [f(i) for i in candidates]


def _is_cube(x, F):
    """
    Check if an element is a cube in GF(2^n)
    """
    if x == F(0): 
        return True
    order = F.order() - 1
    return x**(order // 3) == F(1)


def family12(n, i=None, a=None, b=None, c=None):
    """
    Return the Zheng-Kan-Li-Peng-Tang construction from 2022

    INPUT:
    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``i`` -- (optional) integer with gcd(i, n) = 1;
               If None, all valid i in {1,...,n-1} are tried
    - ``a`` -- (optional) element of GF(2^n) not in GF(q) with a + a^q != 0;
               If None, all valid a in GF(2^n) are tried
    - ``b`` -- (optional) nonzero element of GF(2^n);
               If None, all nonzero b are tried
    - ``c`` -- (optional) nonzero element of GF(2^n);
               If None, all nonzero c are tried
    
    NOTE: When all four parameters are given, returns a single polynomial (or raises an error if no case of Theorem 2 is satisfied).
    When one or more are None, returns a list of all valid polynomials found.

    EXAMPLE USAGE::
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
                        print(res)
                    
    if not res:
        raise TypeError("No valid polynomials found")
    
    return res

def family13():
    pass

