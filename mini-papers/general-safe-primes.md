# k-Safe Prime Construction with Factorization Control

To construct elliptic curve where the order needs specific modular constraints (e.g. $p \equiv 7 \bmod{12}$), the classical safe prime structure $p = 2 q + 1$ may not have the congruence classes needed for primitive roots of unity in the multiplicative group. $p = k \times q + 1$

Let the parameters be:

 * $m$ and $r$, the desired congruence such that $p \equiv r \pmod{m}$
 * $B$, a bound on the small factors of $p-1$.

Then there exists a rejection sampling method to construct primes $p$, $q$ such that:

 * $p = k \times q + 1$, while $k$ may be composite, $|k| \le B$
 * $p \equiv r \pmod{m}$
 * $p$ lies within a given range, $[p_\text{min}, p_\text{max}]$


## Insight

The key efficiency gain comes from precomputing which residue classes of $q \pmod{m}$ could yield the desired congruence $p \equiv r \pmod{m}$. For a given k, we solve the congruence $kq \equiv (r-1) \pmod{m}$ to find all valid residue classes for $q$. This requires using the extended Euclidean algorithm when $gcd(k,m)$ divides $(r-1)$, yielding at most $gcd(k,m)$ distinct solution classes. We then precompute these to retain only those coprime to m (ensuring $q$ has no small factors that would compromise the prime generation). During the sampling phase, we can immediately reject any candidate q that doesn't belong to one of these precomputed residue classes, dramatically reducing the number of primality tests required compared to naive rejection sampling on the final prime $p$.

## Algorithm

 1. For each candidate $k \in [1, B]$:
    a. $s = (r - 1) \bmod{m}$
    b. $d = \gcd(k,m)$
    c. If $d \nmid s$, skip this $k$ (no solutions exist)
    d. Otherwise, solve $kq \equiv s \pmod{m}$:
       - $t = m/d$
       - $u = (k/d)^{-1} \cdot (s/d) \bmod{t}$
    e. Let $x_i \equiv u + it \pmod{m}$ for $i = 0, 1, ..., d-1$
    f. Let $x$ be the subset of solutions where $\gcd(x_i, m) = 1$
    g. Compute $q$ range: $q_{\min} = \lceil(p_{\min} - 1)/k\rceil$, $q_{\max} = \lfloor(p_{\max} - 1)/k\rfloor$
    h. If $q_{\min} \leq q_{\max}$ and valid residue classes exist, mark $k$ as viable

 2. For each viable $k$ (in order of preference):
    a. Generate random primes $q \in [q_{\min}, q_{\max}]$ with $(q \bmod{m}) \in x$
    b. Compute $p = k \times q + 1$
    c. Verify $p$ is prime and $p_{\min} \leq p \leq p_{\max}$
    d. If successful, yield $(p, q, k)$

 3. The resulting primes satisfy:
    - $\text{factor}(p-1) = \text{factor}(k) \cup \{q\}$
    - All small factors of $p-1$ are bounded by $B$
    - $p \equiv r \pmod{m}$ by construction

The factorization bound $B$ allows precise control over the discrete logarithm structure in multiplicative subgroups, while the modular constraints enable compatibility with specialized curve construction methods.

\newpage

## Example

For example. If we want to find primes with a given congruence

```
p mod 6 = p mod (2 * 3)
	p 340282366920938463463374607425593361139
	q 56713727820156410577229101237598893523
	 factor(p-1) 2 * 3 * 56713727820156410577229101237598893523

p mod 18 = p mod (2 * 3^2)
	p 340282366920938463463374607426330524943
	q 18904575940052136859076367079240584719
	 factor(p-1) 2 * 3^2 * 18904575940052136859076367079240584719

p mod 54 = p mod (2 * 3^3)
	p 340282366920938463463374607427247262063
	q 6301525313350712286358789026430504853
	 factor(p-1) 2 * 3^3 * 6301525313350712286358789026430504853

p mod 162 = p mod (2 * 3^4)
	p 340282366920938463463374607425394508779
	q 2100508437783570762119596342132064869
	 factor(p-1) 2 * 3^4 * 2100508437783570762119596342132064869
```

## Implementation

```python
from sage.all import random_prime, is_prime
from math import gcd

def get_q_range(k, p_min, p_max):
    q_min = (p_min - 1 + k - 1) // k
    q_max = (p_max - 1) // k
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
        #a = analyze_k(k, m, r)
        if not len(a):
            continue
        print(f"k={k}")
        for i,(p,q,k,j) in enumerate(emit_safeprimes(k, m, r, p_min, p_max)):
            print("\t", i, j, "p,q,k", p, q, k, "OK")
            print("\t factor(p-1)", factor(p-1))
            print("\t factor(k)", factor(k))
            break
        print()

if __name__ == "__main__":
    main()
```