"""
The problem with this approach is that `p` doesn't fit into our 2**256 - 2**32 - c
window, where c is less than 2**31.
"""

from sage.all import random_prime, is_prime, GF, EllipticCurve, factor
from sage.rings.factorint import factor_trial_division
from math import gcd, isqrt
import gmpy2
gmpy2.get_context().precision = 256

ZETA = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT = lambda n, p: pow(n%p, (p + 1) // 4, p)

ZETA_gmpy2 = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT_gmpy2 = lambda n, p: pow(n%p, (p + 1) // 4, p)

MAX_C = 2**31
#MAX_C = 2**16
bitlen = 255

def check_glv_endomorphism2(curve, p, n, g_base, generator=None):
    # Use provided generator or get default
    if generator is None:
        try:
            generator = curve.gens()[0]
        except:
            return {"supports_glv": False, "reason": "No generator specified and unable to find default generator"}

    """
    # Find a cube root of unity in the base field
    beta = ZETA(3, g_base, p) #g_base**((p-1)//3)
    # Verify beta is a non-trivial cube root of unity
    if beta == 1 or beta**3 != 1:
        assert False
    # Verify beta satisfies the minimal polynomial
    if beta**2 + beta + 1 != 0:
        assert False
    betas = [beta,beta**2]
    """

    #beta_val_6th = ZETA(6, g_base, p) # g_base**((n-1)//6)
    beta_val_3rd = ZETA(3, g_base, p) # g_base**((n-1)//3)
    beta_vals = [pow(beta_val_3rd,i,p) for i in range(1, 3)] # + [pow(beta_val_6th, i, p) for i in range(1, 6)]

    Fq = GF(n)
    Fq_star = Fq.unit_group()
    g_scalar = Fq(Fq_star.gen(0))
    #lambda_val_6th = g_scalar**((n-1)//6)
    lambda_val_3rd = g_scalar**((n-1)//3)
    lambda_vals = [lambda_val_3rd**i for i in range(1, 3)] # [lambda_val_6th**i for i in range(1, 6)] +
    #print('lambda 6', lambda_val_6th, [lambda_val_6th**i for i in range(1, 6)])
    #print('lambda 3', lambda_val_3rd, [lambda_val_3rd**i for i in range(1, 3)])
    # Test both values to see if either works
    #lambda_vals = [lambda_val_6th, lambda_val_3rd]

    for beta_i, beta in enumerate(beta_vals):
        try:
            for lambda_i,lambda_val in enumerate(lambda_vals):
                # Verify lambda^6 = 1
                if lambda_val**6 != 1:
                    continue

                # Check if lambda satisfies the minimal polynomial
                if lambda_val**2 + lambda_val + 1 != 0:
                    continue

                # The critical test: check if lambda*G = (beta*G_x, G_y)
                try:
                    endo_point = curve(beta * generator[0], generator[1])
                    scalar_point = lambda_val * generator
                    if endo_point == scalar_point:
                        return {
                            "supports_glv": True,
                            "beta": beta,
                            "lambda": lambda_val,
                            "verification": "lambda*G == E(beta*G[0], G[1])",
                            "lambda_i": lambda_i,
                            "beta_i": beta_i
                        }
                except Exception:
                    continue

        except Exception as e:
            return {"supports_glv": False, "reason": f"Error finding lambda value: {e}"}

    return {
        "supports_glv": False,
        "reason": "No suitable lambda value found that satisfies the endomorphism relation",
        "beta_candidates": beta_vals,
        "lambda_candidates": lambda_vals
    }

def cornacchia_gmpy2(d, p):
    """Standard Cornacchia algorithm for x^2 + d*y^2 = p"""
    assert p % 12 == 7
    p = gmpy2.mpz(p)
    if d <= 0 or d >= p:
        raise ValueError("invalid input")
    if ZETA_gmpy2(2, -d, p) != 1: raise ValueError("no solution")
    x0 = MODSQRT_gmpy2(-d, p)
    # Choose the larger square root
    if x0 < p // 2:
        x0 = p - x0
    # Extended Euclidean algorithm
    a, b = p, x0
    limit = isqrt(p) # Python 3.8+
    while b > limit:
        a, b = b, a % b
        assert a > 0 and b > 0 # guaranteed positive integers
    remainder = p - b * b
    assert remainder % d == 0  # guaranteed congruence
    c = remainder // d
    t = isqrt(c) # Python 3.8+
    assert t * t == c  # guaranteed exact squares
    return b, t

def sample_primes(bitlen):
    while True:
        p = random_prime(2**bitlen, lbound=2**(bitlen-1))
        if p % 12 == 7:
            yield p

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

def calculate_curve_orders(p, g):
    assert p % 12 == 7
    a, b = cornacchia_gmpy2(3,p)
    assert a%2 == 0 and b%2 == 1
    c, d = a + b, 2 * b
    assert c**2 - c*d + d**2 == p
    assert p + 1 + a - (3*b) == p + 1 + c - (2*d)
    u0 = ZETA_gmpy2(2,a*b*d,p)
    u1 = int((ZETA_gmpy2(3,g,p) * c + d) == 0)
    idx = (int(u0 == 1) * 2) + u1
    result = [0] * 6
    norms_cd = make_norms_cd(c,d,p)
    for i,j in enumerate(ALL_SUB_PATTERNS[idx]):
        result[i] = norms_cd[j]
    return result

def analyze_large_prime_compatibility(k, mod, target_residue, min_log2_q=250):
    """Analyze if k can produce large primes, not just any primes"""
    target_kq_mod = (target_residue - 1) % mod
    # Find all q residues that work
    valid_q_mods = []
    for q_mod in range(mod):
        if (k * q_mod) % mod == target_kq_mod:
            valid_q_mods.append(q_mod)
    # Analyze each residue class for large prime compatibility
    large_prime_compatible = []
    small_prime_only = []
    for q_mod in valid_q_mods:
        if q_mod == 0:
            continue  # q=0 is not prime
        elif q_mod == 1:
            large_prime_compatible.append(q_mod)  # Many large primes ≡ 1 (mod m)
        elif gcd(q_mod, mod) == 1:
            large_prime_compatible.append(q_mod)  # Coprime to mod, can contain large primes
        else:
            # q_mod shares a factor with mod
            # Check if this residue class contains any prime > mod
            #g = gcd(q_mod, mod) #?
            if q_mod <= mod and q_mod in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
                # This residue class only contains one small prime
                if 2**min_log2_q > mod:  # We need much larger primes
                    small_prime_only.append(q_mod)
                else:
                    large_prime_compatible.append(q_mod)
            else:
                small_prime_only.append(q_mod)
    if not large_prime_compatible:
        return False

    # Calculate effective density for large primes
    # This is approximate - we're estimating based on coprimality
    total_coprime_residues = sum(1 for r in range(mod) if gcd(r, mod) == 1)
    large_prime_density = len(large_prime_compatible) / total_coprime_residues

    if large_prime_density == 1:
        print(f"\nAnalyzing k={k} for p ≡ {target_residue} (mod {mod})")
        print(f"Need: k×q ≡ {target_kq_mod} (mod {mod})")
        print(f"All solutions for q (mod {mod}): {valid_q_mods}")
        print(f"Large-prime-compatible q residues: {large_prime_compatible}")
        print(f"Small-prime-only q residues: {small_prime_only}")
        print(f"Large prime density: {len(large_prime_compatible)}/{total_coprime_residues} = {large_prime_density:.3f}")
        return True

def get_prime_range(k, min_log2_q):
    p_min = (2**bitlen) - (2**32) - MAX_C
    p_max = (2**bitlen) - (2**32)
    q_min = max((p_min - 1) // k, int(2**min_log2_q))
    q_max = (p_max - 1) // k
    if q_min < q_max:
        return p_min, p_max, q_min, q_max

def test_k_practically(k, mod, target_residue, min_log2_q):
    """Actually test a k value to see if it finds primes quickly"""
    target_kq_mod = (target_residue - 1) % mod
    possible_q_mods = []
    for q_mod in range(mod):
        if (k * q_mod) % mod == target_kq_mod:
            possible_q_mods.append(q_mod)
    p_min, p_max,q_min,q_max = get_prime_range(k, min_log2_q)
    assert q_min < q_max

    while True:
        q = random_prime(q_max, False, lbound=q_min)

        if q % mod not in possible_q_mods:
            continue

        p = k * q + 1

        if not (p_min <= p <= p_max) or p % mod != target_residue:
            continue

        if is_prime(p):

            yield p,q

from collections import defaultdict

def main():
    viable_k_values = defaultdict(set)
    min_log2_q = bitlen-8
    for need_mod, need_residue in [(12,7),(4,3)]:
        for k in range(1,100):
            range_ok = get_prime_range(k, min_log2_q) is not None
            if not range_ok:
                break
            if analyze_large_prime_compatibility(k, need_mod, need_residue):
                viable_k_values[(need_mod,need_residue)].add(k)
    print(viable_k_values)

    ok = defaultdict(int)
    fail = defaultdict(int)
    dorp = 0
    for curve_order_k in viable_k_values[(12,7)]:
        # First, find safe primes for the order of our curve
        for curve_order, curve_order_bigfactor in test_k_practically(curve_order_k, 12, 7, min_log2_q):
            #assert (curve_order_k*curve_order_bigfactor)+1 == curve_order
            uprime,vprime = cornacchia_gmpy2(3, curve_order)
            assert uprime % 2 == 0 and vprime % 2 == 1

            # Generally speaing, the more representations of c & d we have...
            # As bitlen(q) gets higher, we need to try more closer values
            # As the tests are quite cheap, it's better to have more here
            u, v = 2*uprime, 2*vprime
            possible_values = [
                (v + 1, (u + v) // 2),
                (v + 1, (u + v + 1) // 2),
                (v + 1, (u + v + 2) // 2),
                (v + 1, (u + v + 3) // 2),
                (v - 1, (u + v) // 2),
                (v - 1, (u + v + 1) // 2),
                (v - 1, (u + v - 2) // 2),
                (v - 1, (u + v - 1) // 2),

                ((u + v) // 2, v - 1),
                ((u + v) // 2, v + 1),
                ((u + v + 3) // 2, v + 1),
                ((u + v + 2) // 2, v + 1),
                ((u + v + 1) // 2, v - 1),
                ((u + v - 2) // 2, v - 1),
            ]
            for i,(d,c) in enumerate(possible_values):
                p = int(c**2 - c*d + d**2)
                if p%12==7 and is_prime(p):
                    # We only care if p-1 has small factors below 2**16, or one big factor!
                    p_minus_one_factors = factor_trial_division(p-1, 2**16)
                    if not all([is_prime(_[0]**[1][0]) for _ in p_minus_one_factors]):
                        continue
                    blah = make_norms_cd(c,d,p)
                    assert curve_order in blah, f"Fail {i}"
                    print('i',i, 'q', hex(curve_order), 'p', hex(p), 'c', c, 'd', d, 'is_prime(p)', is_prime(p), 'curve valid?', curve_order in blah)
                    F = GF(p)
                    g_base = F.multiplicative_generator()
                    new_orders = calculate_curve_orders(p, g_base)
                    assert curve_order in new_orders
                    print(' factors', factor(p-1), ",", factor(curve_order-1))
                    E = EllipticCurve(F, [0, g_base])
                    # We get about 50/50% chance of it supporting the GLV endomorphism
                    # Either both work or neither works
                    glv_result = check_glv_endomorphism2(E, p, curve_order, g_base)
                    if glv_result['supports_glv']:
                        ok[i] += 1
                    else:
                        fail[i] += 1
                    dorp += 1
                    print(' GLV? 1', glv_result)
                    print()
                    if dorp > 0 and dorp % 1000 == 0:
                        for morp in sorted(set(list(ok.keys()) + list(fail.keys()))):
                            print(" c d pair", morp, 'success', ok[morp], 'fail', fail[morp])
                        return
                    break

if __name__ == "__main__":
    main()