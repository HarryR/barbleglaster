"""
Shows representations of GLV compatible curves in the coefficient space of the Eisenstein integers
Draws a lattice structure ontop
"""

import matplotlib.pyplot as plt
import math
from sage.all import GF, EllipticCurve, is_prime, next_prime, factor

ZETA = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT = lambda n, p: pow(n%p, (p + 1) // 4, p)

def reduce_to_sublattice(arbitrary_c, arbitrary_d):
    # Force d to be even while preserving essential properties
    reduced_d = arbitrary_d if arbitrary_d % 2 == 0 else arbitrary_d - 1
    # Adjust c accordingly to maintain some invariant
    reduced_c = arbitrary_c + (arbitrary_d - reduced_d) // 2
    return EisensteinInt(reduced_c, reduced_d)

def is_in_cornacchia_slice(self):
    """Check if Eisenstein integer is in the cornacchia output slice"""
    return (self.d >= 0 and
            self.d % 2 == 0 and
            self.c >= self.d // 2)

def check_glv_endomorphism2(curve, p, n, g_base, generator=None):
    # Use provided generator or get default
    if generator is None:
        try:
            generator = curve.gens()[0]
        except:
            return {"supports_glv": False, "reason": "No generator specified and unable to find default generator"}

    beta_val_3rd = ZETA(3, g_base, p)
    beta_vals = [pow(beta_val_3rd,i,p) for i in range(1, 3)]

    Fq = GF(n)
    #Fq_star = Fq.unit_group()
    g_scalar = Fq.multiplicative_generator() # Fq(Fq_star.gen(0))
    lambda_val_3rd = g_scalar**((n-1)//3)
    lambda_vals = [lambda_val_3rd**i for i in range(1, 3)]

    for beta_i, beta in enumerate(beta_vals):
        try:
            for lambda_i,lambda_val in enumerate(lambda_vals):
                # Verify lambda^6 = 1
                #if lambda_val**6 != 1:
                #    continue

                # Check if lambda satisfies the minimal polynomial
                #if lambda_val**2 + lambda_val + 1 != 0:
                #    continue

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


def is_perfect_square(n):
    """Check if n is a perfect square using integer square root."""
    if n < 0:
        return False
    root = math.isqrt(n)
    return root * root == n

class BinaryQuadraticForm:
    def __init__(self, a, b, c):  # ax² + bxy + cy²
        self.a, self.b, self.c = a, b, c

    def discriminant(self):
        return self.b**2 - 4*self.a*self.c

class EisensteinInt:
    def __init__(self, c, d):
        self.c = c  # coefficient of 1
        self.d = d  # coefficient of ω

    @property
    def wedge(self):
        c, d = self.c, self.d
        if c == 0 and d == 0:
            return None  # or 0, depending on convention
        # The 6 wedges are determined by comparing c, d, and c+d
        # This comes from the geometry of the hexagonal lattice
        if d > 0:
            if c > d:
                return 0  # 0° to 60°
            elif c > 0:
                return 1  # 60° to 120°
            return 2  # 120° to 180°
        if c < d:
            return 3  # 180° to 240°
        elif c < 0:
            return 4  # 240° to 300°
        return 5  # 300° to 360°

    def __add__(self, other:'EisensteinInt'):
        return EisensteinInt(self.c + other.c, self.d + other.d)

    def __sub__(self, other):
        return EisensteinInt(self.c - other.c, self.d - other.d)

    def __mul__(self, other):
        # (c₁ + d₁ω)(c₂ + d₂ω) = (c₁c₂ - d₁d₂) + (c₁d₂ + d₁c₂ - d₁d₂)ω
        # since ω² = ω - 1
        return EisensteinInt(
            self.c * other.c - self.d * other.d,
            self.c * other.d + self.d * other.c - self.d * other.d
        )

    def conjugate(self):
        # Complex conjugate of c + dω where ω = e^(2πi/3)
        return EisensteinInt(self.c - self.d, -self.d)

    def norm(self):
        return self.c*self.c - self.c*self.d + self.d*self.d

    def distance_norm(self, other:'EisensteinInt'):
        diff = EisensteinInt(self.c - other.c, self.d - other.d)
        return diff.norm()  # Eisenstein norm of difference

    def distance_euclidean_squared(self, other:'EisensteinInt'):
        # |z1 - z2|² in complex plane
        return (self.c - other.c)**2 + (self.d - other.d)**2 + (self.c - other.c)*(self.d - other.d)

    @property
    def x(self):
        return 2 * self.c - self.d

    @property
    def y(self):
        return self.d

    @classmethod
    def from_norm(cls, n):
        a, b = cornacchia(3, n)
        return cls(a + b, 2 * b)

    def rotate_60(self):
        # Multiplication by ω = (-1 + i√3)/2
        return EisensteinInt(-self.d, self.c + self.d)

    # Easy 6-fold symmetry operations
    def all_6_rotations(self):
        result = self
        rotations = [result]
        for _ in range(5):
            result = result.rotate_60()
            rotations.append(result)
        return rotations

    """
    def is_gauss_representable(self):
        # Check if an Eisenstein integer can be converted to Gaussian with same norm.
        if self.d % 2 != 0:
            return False
        b = self.d // 2
        # Need 3b² to be a perfect square for the conversion to work
        return is_perfect_square(3 * b * b)

    def to_gauss_same_norm(self):
        # Convert to Gaussian integer with same norm, if possible.
        if not self.is_gauss_representable():
            raise ValueError("Cannot convert to Gaussian integer with same norm")
        b = self.d // 2
        a = self.c - b
        # Verify the decomposition
        n = self.norm()
        if a*a + 3*b*b != n:
            raise ValueError("Transformation error")
        # Find b' such that a² + (b')² = a² + 3b²
        b_prime = math.isqrt(3 * b * b)
        return GaussInt(a, b_prime)
    """

"""
class GaussInt:
    def __init__(self, a, b):
        self.a = a  # real part
        self.b = b  # imaginary part (coefficient of i)

    def to_eisenstein(self):
        return EisensteinInt(self.a + self.b, 2 * self.b)

    @property
    def quadrant(self):
        a, b = self.a, self.b
        if a == 0 and b == 0:
            return None  # origin
        # The 4 quadrants in the complex plane
        if a > 0 and b >= 0:
            return 0  # First quadrant (0° to 90°)
        elif a <= 0 and b > 0:
            return 1  # Second quadrant (90° to 180°)
        elif a < 0 and b <= 0:
            return 2  # Third quadrant (180° to 270°)
        else:  # a >= 0 and b < 0
            return 3  # Fourth quadrant (270° to 360°)

    def __add__(self, other):
        return GaussInt(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        return GaussInt(self.a - other.a, self.b - other.b)

    def __mul__(self, other):
        # (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        return GaussInt(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a
        )

    def norm(self):
        return self.a * self.a + self.b * self.b

    def distance_norm(self, other: 'GaussInt'):
        diff = self - other
        return diff.norm()  # Gaussian norm of difference

    def distance_euclidean_squared(self, other: 'GaussInt'):
        # |z1 - z2|² in complex plane (same as norm for Gaussian integers)
        return (self.a - other.a)**2 + (self.b - other.b)**2

    @property
    def x(self):
        return self.a  # real coordinate

    @property
    def y(self):
        return self.b  # imaginary coordinate

    def conjugate(self):
        return GaussInt(self.a, -self.b)

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        elif self.a == 0:
            return f"{self.b}i" if self.b != 1 else "i"
        else:
            sign = "+" if self.b > 0 else "-"
            b_abs = abs(self.b)
            b_str = "i" if b_abs == 1 else f"{b_abs}i"
            return f"{self.a} {sign} {b_str}"

    def __repr__(self):
        return f"GaussInt({self.a}, {self.b})"
"""

omega = EisensteinInt(0, 1)  # 0 + 1*ω
one = EisensteinInt(1, 0)  # 1 + 0*ω = 1

j0curve_trace_coefficients = [(-2,1), (-1,-1), (1,-2), (2,-1), (1,1), (-1,2)]

def make_norms_cd(c,d,p):
    return [p + 1 + (x_0*c) + (x_1*d) for x_0, x_1 in [
        (-2,1), (-1,-1), (1,-2), (2,-1), (1,1), (-1,2)
    ]]

ALL_SUB_PATTERNS = [
    (0, 1, 2, 3, 4, 5),
    (0, 5, 4, 3, 2, 1),
    (3, 4, 5, 0, 1, 2),
    (3, 2, 1, 0, 5, 4),
]

def is_diagonal_intersection(c, d):
    # Must be on grid (either coordinate ≡ 2 mod 4)
    if not (c % 4 == 2 or d % 4 == 2):
        return False
    # Check if on slope 2 line: k = (2x - y + 6)/12
    if (2*c - d + 6) % 12 == 0:
        if (c + d - 6) % 6 == 0:
            return True
    return False

def calculate_curve_orders(p, g):
    assert p % 12 == 7
    a, b = cornacchia(3,p)
    assert a%2 == 0 and b%2 == 1
    c, d = a + b, 2 * b
    assert c**2 - c*d + d**2 == p
    assert p + 1 + a - (3*b) == p + 1 + c - (2*d)
    u0 = ZETA(2,a*b*d,p)
    u1 = int((ZETA(3,g,p) * c + d) == 0)
    idx = (int(u0 == 1) * 2) + u1
    result = [0] * 6
    norms_cd = make_norms_cd(c,d,p)
    for i,j in enumerate(ALL_SUB_PATTERNS[idx]):
        result[i] = norms_cd[j]

    diag_x = (2*c - d + 6) % 12 # 10 for left, 2 for right
    diag_y = (c + d - 6) % 6  # 5 for left, 1 for right
    assert diag_x in [10,2]
    assert diag_y in [5,1]
    assert (diag_x,diag_y) in [(10,5),(2,1)]
    is_left = (diag_x,diag_y) == (10,5)
    assert int(is_left) == ((c+d) % 3) - 1
    diagonal_point = (c+1,d) if is_left else (c-1,d)

    ie = EisensteinInt(*diagonal_point)
    #print("Anchor", ie.norm(), ie.c, ie.d)

    """
    assert is_diagonal_intersection(*diagonal_point)
    print("Intercept point", diagonal_point)
    print("\tSlope 2 ",
        'Y-Intercept', diagonal_point[1] - (2 * diagonal_point[0]),
        "X-Intercept", ((2 * diagonal_point[0]) - diagonal_point[1])/2)
    print("\tSlope -1 intersect",
        "Y-Intercept", sum(diagonal_point),
        "X-Intercept", sum(diagonal_point))
    print()
    """

    return result

def pp(ax, p:EisensteinInt, alpha, fmt='o', label=None):
    n = p.norm()
    if n % 6 == 1 or n%3 == 1:
        if n%12==7:
            color = 'purple'
        elif n%6 == 1:
            color = 'blue'
        else:
            assert n%3==1
            color = 'green'
    else:
        color = 'grey'
    ax.plot(p.c, p.d, fmt, color=color, alpha=alpha)
    if label is not None:
        ax.annotate(label, (p.c, p.d),
                    xytext=(0, 0), textcoords='offset points',
                    fontsize=8, alpha=alpha)

def point_id(n:EisensteinInt):
    a=n.c%12
    b=n.d%12
    pn = {
        (10,3): 0,
        (9,2): 1,
        (7,1): 2,
        (6,1): 3,
        (5,2): 4,
        (5,3): 5,
        (6,5): 6,
        (7,6): 7,
        (9,7): 8,
        (10,7): 9,
        (11,6): 10,
        (11,5): 11,

        (10,1): 6,
        (11,2): 7,

        (6,7): 0,
        (5,6): 1,
        (3,5): 2,
        (2,5): 3,
        (1,6): 4,
        (1,7): 5,
        (2,9): 6,
        (3,10): 7,
        (5,11): 8,
        (6,11): 9,
        (7,10): 10,
        (7,9): 11,

        (2,11): 0,
        (1,10): 1,

        (1,3): 8,
        (2,3): 9,
        (3,2): 10,
        (3,1): 11,

        (11,9): 2,
        (10,9): 3,
        (9,10): 4,
        (9,11): 5,
    }
    return pn[(a,b)]

def main():
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.set_xlabel('Real part: a-b/2', fontsize=12)
    ax.set_ylabel('Imaginary part: b*sqrt(3)/2', fontsize=12)
    ax.set_title(f'Points: a + b·ω where ω = e^(2πi/3)', fontsize=14)

    max_blah = 25
    max_p = 0

    found_curves = 0
    p = 2
    while found_curves < max_blah:
        p = next_prime(p)
        if p % 12 != 7:
            continue
        if not is_prime(p):
            continue
        F = GF(p)
        g = F.multiplicative_generator()
        orders = calculate_curve_orders(p, g)
        if orders[1] == p or orders[5] == p: # skip supersingular curves
            continue
        orders_prime = [is_prime(_) for _ in orders]
        g1_prime = orders_prime[1]
        g5_prime = orders_prime[5]
        if any([g1_prime,g5_prime]):
            glv_ok = 0
            glv_orders = []
            glv_idx = []
            for i,n in enumerate(orders):
                if i in [1,5]:
                    if not orders_prime[i]:
                        continue
                    E = EllipticCurve(F, [0, g**i])
                    assert E.order() == n
                    result = check_glv_endomorphism2(E, p, n, g)
                    if result['supports_glv']:
                        glv_idx.append(i)
                        glv_ok += 1
                        glv_orders.append(n)
            assert glv_ok != 0
            glv_ok = glv_ok != 0
            prime = EisensteinInt.from_norm(p)
            max_p = max(prime.c, prime.d, max_p)
            pp(ax, prime, 1.0, '+')
            prime_neighbors = [
                EisensteinInt(prime.c+i, prime.d+j)
                for i,j in [
                    (0,1),(1,1),(1,0),
                    (-1,0),(-1,-1),(0,-1)
                ]
            ]
            for neighbor in prime_neighbors:
                n = neighbor.norm()
                assert n in orders
                fmt='.' if orders.index(n) == 1 else 'o'
                oi = orders.index(n)
                if oi in [1,5] and oi in glv_idx:
                    alpha = 0.9
                    ax.plot([neighbor.c, prime.c], [neighbor.d, prime.d], linewidth=1, color='grey', alpha=0.75)
                else:
                    if oi in [1,5]:
                        alpha = 0.2
                        ax.plot([neighbor.c, prime.c], [neighbor.d, prime.d], linewidth=1, color='grey', alpha=0.5)
                    else:
                        alpha = 0.05
                pp(ax, neighbor, alpha, fmt=fmt)
            found_curves += 1

    """
    total_primes = 0
    adjacent_primes = 0
    for c in range(0,max_p):
        for d in range(0,max_p):
            e = EisensteinInt(c,d)
            n = e.norm()
            if n % 12 == 7 and is_prime(n):
                total_primes += 1
                pp(ax, e, 1, '*')
                if is_diagonal_intersection(e.c+1, e.d) or is_diagonal_intersection(e.c-1, e.d):
                    adjacent_primes += 1
    print(max_p, total_primes, adjacent_primes)
    """



    for c in range(0,max_p):
        for d in range(0,max_p):
            if is_diagonal_intersection(c,d):
                e = EisensteinInt(c,d)
                pp(ax, e, 0.2, '.')

    """
            e = EisensteinInt(c,d)
            n = e.norm()
            if n%12==7:
                label = None
                #label = point_id(e)
                pp(ax, e, 0.2, '.', label=label)
            else:
                if n%6==1:
                    ignore_points = [
                        (3,8),(4,8),(5,9),(5,8),(4,7),(3,7),(4,9),
                        (11,0),(0,1),(1,1),(1,0),(11,11),(0,11),
                        (7,4),(7,3),(8,3),(9,4),(9,5),(8,5),
                    ]
                    if (c%12,d%12) not in ignore_points:
                        pp(ax, e, 0.2 if n%6==1 or n%3==1 else 0.05, '.') #, label=f"{hex(e.c%12)[2:]}{hex(e.d%12)[2:]}")
                    else:
                        pp(ax, e, 0.07, '|')
                elif n%3==1:
                    pp(ax, e, 0.07, '_')
                elif is_diagonal_intersection(c,d):
                    pp(ax, e, 1, '.')

            # Mark the centers of each circle
            #if (c%12 == 0 and d%12 == 0) or (c%12==4 and d%12==8) or (c%12==8 and d%12==4):
            #    ax.annotate(f'{(c,d)}', (c, d),
            #        xytext=(c, d),
            #        fontsize=8, alpha=0.7)
        """

    for i in range(max_p):
        """
        if i % 12 == 0:
            ax.axhline(y=i, color='k', linewidth=0.5, alpha=0.1)
            ax.axvline(x=i, color='k', linewidth=0.5, alpha=0.1)
        elif i % 4 == 0:
            ax.axhline(y=i, color='k', linewidth=0.5, alpha=0.1, linestyle=':')
            ax.axvline(x=i, color='k', linewidth=0.5, alpha=0.1, linestyle=':')
        """
        if i % 12 == 6:
            ax.axhline(y=i, color='k', linewidth=0.5, alpha=0.1)
        if i % 18 == 6:
            ax.axvline(x=i, color='k', linewidth=0.5, alpha=0.1)

        if i % 2 == 0:
            ax.axvline(x=i, color='k', linewidth=0.1, alpha=0.1)
            ax.axhline(y=i, color='k', linewidth=0.1, alpha=0.1)
        if i%6 == 0:
            ax.axline((i,6), color='k', slope=2, alpha=0.2, linewidth=0.5)
            ax.axline((i,6), color='k', slope=-1, alpha=0.2, linewidth=0.5)

    plt.tight_layout()
    plt.savefig("graphs/eisenstein_lattice.png", dpi=300)
    plt.savefig("graphs/eisenstein_lattice.svg", dpi=300)

# Main execution
if __name__ == "__main__":
    main()
