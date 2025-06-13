#!/usr/bin/env python3

import sys
import gmpy2
import math
from sage.all import random_prime, is_prime, GF
from sage.rings.factorint import factor_trial_division
from lib_eta import factors_metrics, factors_str

gmpy2.get_context().precision = 256

ZETA_gmpy2 = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT_gmpy2 = lambda n, p: pow(n%p, (p + 1) // 4, p)

def cornacchia_gmpy2(d, p):
    """Standard Cornacchia algorithm for x^2 + d*y^2 = p"""
    assert p % 12 == 7
    p = gmpy2.mpz(p)
    if d <= 0 or d >= p:
        raise ValueError("invalid input")
    if ZETA_gmpy2(2, -d, p) != 1:
        raise ValueError("no solution")
    x0 = MODSQRT_gmpy2(-d, p)
    # Choose the larger square root
    if x0 < p // 2:
        x0 = p - x0
    # Extended Euclidean algorithm
    a, b = p, x0
    limit = math.isqrt(p) # Python 3.8+
    while b > limit:
        a, b = b, a % b
        assert a > 0 and b > 0 # guaranteed positive integers
    remainder = p - b * b
    assert remainder % d == 0  # guaranteed congruence
    c = remainder // d
    t = math.isqrt(c) # Python 3.8+
    assert t * t == c  # guaranteed exact squares
    return b, t

def make_terms_cd(c,d):
    return [(x_0*c) + (x_1*d) for x_0, x_1 in [
        (-2,1), (-1,-1), (1,-2), (2,-1), (1,1), (-1,2)
    ]]

def make_norms_cd(c,d,p):
    return [p + 1 + _ for _ in make_terms_cd(c,d)]

ALL_SUB_PATTERNS = [
    (0, 1, 2, 3, 4, 5),
    (0, 5, 4, 3, 2, 1),
    (3, 4, 5, 0, 1, 2),
    (3, 2, 1, 0, 5, 4),
]

def calculate_curve_orders(p, g, a, b):
    assert p % 12 == 7
    g = gmpy2.mpz(g)
    p = gmpy2.mpz(p)
    a = gmpy2.mpz(a)
    b = gmpy2.mpz(b)
    assert a%2 == 0 and b%2 == 1
    c, d = a + b, 2 * b
    assert c**2 - c*d + d**2 == p
    assert p + 1 + a - (3*b) == p + 1 + c - (2*d)
    u0 = int(((c+d) % 3) == 2)
    u1 = int(gmpy2.mod(gmpy2.fma(ZETA_gmpy2(3,g,p), c, d), p) == 0)
    idx = (u0 * 2) + u1
    result = [0] * 6
    norms_cd = make_norms_cd(c,d,p)
    for i,j in enumerate(ALL_SUB_PATTERNS[idx]):
        result[i] = norms_cd[j]
    return result

def main(bitlen):
    bitlen = int(bitlen)
    i,j = 0, 0
    avg_i = set()
    max_fs = 0
    while True:
        i += 1
        p = random_prime(2**bitlen)
        if (p-1) % 6 != 0:
            continue
        pd6 = p//6
        if not is_prime(pd6):
            continue
        avg_i.add(i)
        j += 1
        i = 0
        a,b = cornacchia_gmpy2(3, p)
        Fp = GF(p)
        g = int(Fp.multiplicative_generator())
        orders = [int(_) for _ in calculate_curve_orders(p, g, a, b)]
        if is_prime(orders[1]) and is_prime(orders[5]):
            twists_factors = []
            print((i,j), p, 0 if len(avg_i) == 0 else sum(avg_i) / len(avg_i))
            for i,(ip,n) in enumerate([(is_prime(int(_)),_) for _ in orders]):
                factors = [(int(prime),int(power)) for prime,power in factor(n)]
                twists_factors.append(factors_metrics(factors)[0])
                factors_metrics(factors)
                print("\t", ip, n, factors_str)

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
