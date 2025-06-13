from sage.all import random_prime, GF, EllipticCurve, GF, is_prime
import math
import time
import sys
ZETA = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT = lambda n, p: pow(n%p, (p + 1) // 4, p)

def cornacchia(d, p):
    """Standard Cornacchia algorithm for x^2 + d*y^2 = p"""
    assert p % 12 == 7
    if d <= 0 or d >= p:
        raise ValueError("invalid input")
    #if ZETA(2, -d, p) != 1: raise ValueError("no solution")
    x0 = MODSQRT(-d, p)
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

def sample_primes(bitlen):
    while True:
        p = random_prime(2**bitlen, lbound=2**(bitlen-1))
        if p % 12 == 7:
            yield p

def make_norms_cd(c,d,p):
    return [p + 1 + (x_0*c) + (x_1*d) for x_0, x_1 in [
        (-2,1), (-1,-1), (1,-2), (2,-1), (1,1), (-1,2)
    ]]

def curve_orders_eisenstein_coords(c,d):
    return [(c+i, d+j) for i,j in [(-1,0), (-1,-1), (0,-1), (1,0), (1,1), (0,1)]]

def eisenstein_norm(c,d):
    return (c**2) + (d**2) - (c * d)

def cornacchia_norm(a,b):
    return (a**2) + (3 * (b**2))

ALL_SUB_PATTERNS = [
    (0, 1, 2, 3, 4, 5),
    (0, 5, 4, 3, 2, 1),
    (3, 4, 5, 0, 1, 2),
    (3, 2, 1, 0, 5, 4),
]

def calculate_curve_orders(p, g):
    assert p % 12 == 7
    a, b = cornacchia(3,p)
    assert a%2 == 0 and b%2 == 1
    c, d = a + b, 2 * b
    assert eisenstein_norm(c,d) == p # c**2 - c*d + d**2
    assert p + 1 + a - (3*b) == p + 1 + c - (2*d)
    u0 = int(((c+d) % 3) == 2)
    u1 = int((ZETA(3,g,p) * c + d) == 0)
    idx = (u0 * 2) + u1
    result = [0] * 6
    norms_cd = make_norms_cd(c,d,p)
    prime_neighbor_norms = [eisenstein_norm(*_) for _ in curve_orders_eisenstein_coords(c,d)]
    assert prime_neighbor_norms == norms_cd
    for i,j in enumerate(ALL_SUB_PATTERNS[idx]):
        result[i] = norms_cd[j]
    return result

def find_generator(g,p):
    g = int(g)
    x = 1
    while True:
        yy = (pow(x,3,p) + g) % p
        y = MODSQRT(yy, p)
        if (y*y) == yy:
            if y & 1:
                y = p - y
            E = EllipticCurve(GF(p), [0, g])
            if E.point((x,y)).order() == E.order():
                return x,y
        x += 1

def example_secp256k1():
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    F = GF(p)
    g = F.multiplicative_generator()
    q = 115792089237316195423570985008687907852837564279074904382605163141518161494337
    G = find_generator(g, p)
    orders = calculate_curve_orders(p,g)
    assert orders.index(q) == 5
    print('secpk256k1')
    print('         p', p, hex(p))
    print('         g', g)
    print('         G', G)
    print('    orders')
    for i, n in enumerate(orders):
        k = GF(n)(p).multiplicative_order() if is_prime(n) else None
        print(f'         {i}', n, hex(n), 'prime' if is_prime(n) else '', f'log2(embedding degree)={int(math.log2(k))}' if k is not None else '')
    print('       idx', orders.index(q))
    print()

def bitlen_tests():
    print('Running Tests', end='')
    with open('data/2-vs-sage.csv', 'w') as handle:
        handle.write("bitlen,count,sage,ours\n")
        for bitlen in range(32,192):
            i = 0
            totals = [0,0,0]
            for p in sample_primes(bitlen):
                i += 1
                if i > 50:
                    break
                F = GF(p)
                g = F.multiplicative_generator()
                curves = [EllipticCurve(F, [0, g**i]) for i in range(6)]
                start1 = time.perf_counter()
                orders = [_.order() for _ in curves]
                total1 = time.perf_counter() - start1
                start2 = time.perf_counter()
                predicted_orders = calculate_curve_orders(p,g)
                total2 = time.perf_counter() - start2
                assert predicted_orders == orders
                totals[0] += 1
                totals[1] += total1
                totals[2] += total2
            print('.', end='', flush=True)
            handle.write(f"{bitlen},{totals[0]},{totals[1]},{totals[2]}\n")
            handle.flush()
    print('TESTS OK ^_^')

def main():
    example_secp256k1()
    #bitlen_tests()


if __name__ == "__main__":
    main()
