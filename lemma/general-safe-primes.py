from sage.all import random_prime, is_prime
from math import gcd

def get_q_range(k, p_min, p_max):
    q_min = (p_min - 1 + k - 1) // k  # Ceiling division trick
    q_max = (p_max - 1) // k          # Floor division
    if q_min <= q_max:
        return q_min, q_max

def generate_k_values(B):
    k_values = {1}
    for p in [p for p in range(2, B + 1) if is_prime(p)]:
        new_values = set()
        for k in k_values:
            power = p
            while k * power <= B:
                new_values.add(k * power)
                power *= p
        k_values.update(new_values)
    return sorted(k_values)

def analyze_k(k, m, r):
    s = (r - 1) % m
    d = gcd(k, m)
    if s % d != 0:
        return []
    t = m // d
    u = (pow(k // d, -1, t) * (s // d)) % t
    a = [(u + i * t) % m for i in range(d)]
    return [_ for _ in a if gcd(_,m) == 1]

def emit_safeprimes(k, mod, target_residue, p_min, p_max, max_tries=1000):
    possible_q_mods = analyze_k(k, mod, target_residue)
    q_min, q_max = get_q_range(k, p_min, p_max)
    assert q_min < q_max
    i = 0
    while i < max_tries:
        q = random_prime(q_max, False, lbound=q_min)
        i += 1
        if q % mod not in possible_q_mods:
            continue
        p = k * q + 1
        assert p_min <= p <= p_max
        assert p % mod == target_residue
        if is_prime(p):
            yield p,q,k,i

def main():
    from sage.all import factor
    p_max = 2**128 - 2**32
    p_min = p_max - 2**31
    m = 12
    r = 7
    for k in generate_k_values(2**8):
        if k not in [6,18,54,162]:
            continue
        a = analyze_k(k, m, r)
        if not len(a):
            continue
        print(f"p mod {k} = p mod ({factor(k)})")
        for i,(p,q,k,j) in enumerate(emit_safeprimes(k, m, r, p_min, p_max)):
            print("\tp", p)
            print("\tq", q)
            #print("\t", i, j, "p,q,k", p, q, k, "OK")
            print("\t factor(p-1)", factor(p-1))
            #print("\t factor(k)", factor(k))
            break
        print()

if __name__ == "__main__":
    main()