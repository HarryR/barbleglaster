import sys
import time
import random
import math
from collections import namedtuple
from sage.all import FiniteField, EllipticCurve, sqrt as sage_sqrt, GF

import gmpy2
gmpy2.get_context().precision = 256

def ZETA_gmpy2(n,x,p):
    return gmpy2.powmod(x, ((p-1)//n), p)

def MODSQRT_gmpy2(n,p):
    return gmpy2.powmod(n, (p + 1)//4, p)

def _glv_shift_count(n):
    return int(math.log2(n)*1.5)

def rnddiv2(v):
    return v+1 if v&1 else v>>1

def fp_conj(x,p):
    return (-int(x) - 1) % p

def find_generator(g,p,E:EllipticCurve):
    p, g, x = (gmpy2.mpz(int(_)) for _ in (p, g, 1))
    while True:
        yy = (gmpy2.powmod(x,3,p) + g) % p
        y = MODSQRT_gmpy2(yy, p)
        if (y*y) % p == yy:
            if y & 1:
                y = p - y
            if E.point((x,y)).order() == E.order():
                return int(x),int(y)
        x += 1

def _glv_find_split_constants_explicit_tof(p:int, E:EllipticCurve):
    """Find constants for secp256k1_scalar_split_lamdba using the trace of Frobenius.

    See Benjamin Smith: "Easy scalar decompositions for efficient scalar multiplication on
    elliptic curves and genus 2 Jacobians" (https://eprint.iacr.org/2013/672), Example 2
    """
    assert p % 3 == 1
    assert E.j_invariant() == 0
    # TODO: calculate without Sage
    t = int(E.trace_of_frobenius())
    c = int(sage_sqrt((4*p - t**2)//3))
    b1 = c
    b2 = (1 - (t - c)//2) % E.order()
    return b1, b2

def _glv_calc_g1_g2(p, n, E:EllipticCurve, beta_val, lambda_val, shift_count):
    b1, b2 = _glv_find_split_constants_explicit_tof(p, E)

    beta_val = fp_conj(beta_val, p)
    lambda_val = fp_conj(lambda_val, n)

    # Python's round() is off by 1
    quotient_g1 = (2**shift_count)*(-b2)//n
    remainder_g1 = (2**shift_count)*(-b2)%n
    g1 = (quotient_g1 + (1 if remainder_g1 >= n//2 else 0)) % (2**(p.bit_length()))

    # Python's round() is off by one
    quotient_g2 = (2**shift_count)*(b1)//n
    remainder_g2 = (2**shift_count)*(b1)%n
    g2 = (quotient_g2 + (1 if remainder_g2 >= n//2 else 0)) % (2**(p.bit_length()))

    return b1,b2,g1,g2

def _glv_decompose(n, k, shift_count, g1, g2, b1, b2, lambda_val):
    """Decompose scalar k into k1 and k2 using GLV method."""
    # Ensure k is properly reduced modulo n
    k = k % n

    # Compute c1 and c2 using precomputed constants and bit operations
    # This replaces division with multiplication and right shift
    c1 = (k * int(g1)) >> shift_count
    c2 = (k * int(g2)) >> shift_count

    # Handle final bit for rounding as per the paper
    if (k * int(g1)) & (1 << (shift_count-1)):
        c1 += 1
    if (k * int(g2)) & (1 << (shift_count-1)):
        c2 += 1

    # Compute k2 = -c1·(-b1) - c2·(-b2)
    k2 = (c1 * int(b1) + c2 * int(b2)) % n

    # Compute k1 = k - k2·λ mod n
    k1 = (k - (k2 * int(lambda_val)) % n) % n

    # Ensure k1 and k2 are properly minimized
    if k1 > (n >> 1):
        k1 = k1 - n
    if k2 > (n >> 1):
        k2 = k2 - n

    #assert k1 + (lambda_val*k2) == k

    return k1, k2

def _glv_score(k, n):
   k1,k2 = k
   if k1 == 0 and k2 == 0:
       return 0.0
   log2_n = math.log2(n)
   target = log2_n / 2
   k1_bits = abs(k1).bit_length() if k1 != 0 else 0
   k2_bits = abs(k2).bit_length() if k2 != 0 else 0
   max_bits = max(k1_bits, k2_bits)
   return max(0.0, 1.0 - (max_bits - target) / target)

def _glv_decompose_args(E, p, n, beta_val, lambda_val, flip_b1, flip_b2):
    shift_count = _glv_shift_count(n)
    b1,b2,g1,g2 = _glv_calc_g1_g2(p, n, E, beta_val, lambda_val, shift_count)
    if flip_b1:
        b1 = -b1
    if flip_b2:
        b2 = -b2
    count, total = 0, 0
    for i in range(2,50):
        k = int((n-1)//i)
        k1,k2 = _glv_decompose(n, k, shift_count, g1, g2, b1, b2, lambda_val)
        total += (_glv_score((k1,k2),n) + _glv_score((b1,b2),n)) / 2
        count += 1
    return (total/count, (b1,b2,g1,g2))

def _glv_decompose_efficiency(E:EllipticCurve, p, n, beta_val, lambda_val):
    # Sort by efficiency, while flipping the parameters
    # return only the most efficient
    return sorted([
        _glv_decompose_args(E, p, n, beta_val, lambda_val, flip_b1, flip_b2)
        for flip_b2 in [True,False]
        for flip_b1 in [True,False]
    ], key=lambda _: _[0])[-1]

def _glv_check(curve:EllipticCurve, p, n, generator):
    generator = curve.point(generator)
    seen_betas = set()
    results = []
    for beta_i in range(2,1000):
        beta_val = ZETA_gmpy2(3, beta_i, p)
        if beta_val == 1 or beta_val in seen_betas:
            continue
        seen_betas.add(beta_val)
        seen_lambdas = set()
        for lambda_i in range(2,1000):
            lambda_val = ZETA_gmpy2(3, lambda_i, n)
            if lambda_val == 1 or lambda_val in seen_lambdas:
                continue
            seen_lambdas.add(lambda_val)
            endo_point = curve(beta_val * generator[0], generator[1])
            scalar_point = lambda_val * generator
            if endo_point == scalar_point:
                score, decompose_params = _glv_decompose_efficiency(curve, p, n, beta_val, lambda_val)
                results.append((score, (beta_i, beta_val, lambda_i, lambda_val, decompose_params)))
    if len(results) == 0:
        return None
    return sorted(results, key=lambda _:_[0])[-1][1]

class Scalar(namedtuple('_Scalar', ['value', 'n'])):
    """Class representing a scalar in the field of curve order."""
    def __add__(self, other):
        return Scalar((self.value + other.value) % self.n, self.n)
    def __sub__(self, other):
        return Scalar((self.value - other.value) % self.n, self.n)
    def __mul__(self, other):
        if isinstance(other, Scalar):
            return Scalar((self.value * other.value) % self.n, self.n)
        else:  # Assume it's an integer
            assert isinstance(other, int)
            return Scalar((self.value * (other % self.n)) % self.n, self.n)
    def __neg__(self):
        return Scalar((-self.value) % self.n, self.n)
    def __eq__(self, other):
        return self.value == other.value and self.n == other.n
    def __str__(self):
        return hex(self.value)
    def __int__(self):
        return int(self.value)

class FieldElement(namedtuple('_FieldElement', ['value', 'p'])):
    """Class representing an element in the prime field."""
    def __add__(self, other):
        return FieldElement((self.value + other.value) % self.p, self.p)
    def __sub__(self, other):
        return FieldElement((self.value - other.value) % self.p, self.p)
    def __mul__(self, other):
        return FieldElement((self.value * other.value) % self.p, self.p)
    def __pow__(self, exp):
        return FieldElement(pow(self.value, exp, self.p), self.p)
    def __eq__(self, other):
        return self.value == other.value and self.p == other.p
    def __str__(self):
        return hex(self.value)
    def __int__(self):
        return int(self.value)

class Point(namedtuple('_Point', ['x','y','a','b','p','is_infinity'])):
    """Class representing a point on an elliptic curve."""
    @classmethod
    def infinity(cls, a, b, p):
        """Return the point at infinity."""
        return cls(0, 0, a, b, p, True)

    def is_on_curve(self):
        """Check if the point lies on the curve."""
        if self.is_infinity:
            return True
        left = (self.y * self.y) % self.p
        right = (self.x * self.x * self.x + self.a * self.x + self.b) % self.p
        return left == right

    def __eq__(self, other):
        if self.a != other.a or self.b != other.b or self.p != other.p:
            return False
        if self.is_infinity:
            return other.is_infinity
        if other.is_infinity:
            return False
        return self.x == other.x and self.y == other.y

    def __str__(self):
        if self.is_infinity:
            return "Point(infinity)"
        return f"Point({hex(self.x)}, {hex(self.y)})"

    def __add__(self, other):
        """Add two points using the elliptic curve group law."""
        if self.is_infinity:
            return other
        if other.is_infinity:
            return self

        # Point doubling
        if self.x == other.x:
            if (self.y + other.y) % self.p == 0:
                return Point.infinity(self.a, self.b, self.p)
            else:
                # Compute the slope of the tangent line
                lambda_val = ((3 * self.x * self.x + self.a) * pow(2 * self.y, -1, self.p)) % self.p
        else:
            # Point addition
            lambda_val = ((other.y - self.y) * pow(other.x - self.x, -1, self.p)) % self.p

        x3 = (lambda_val * lambda_val - self.x - other.x) % self.p
        y3 = (lambda_val * (self.x - x3) - self.y) % self.p
        return Point(x3, y3, self.a, self.b, self.p, False)

    def scalar_mul(self, k):
        """Multiply point by scalar k using double-and-add algorithm."""
        result = Point.infinity(self.a, self.b, self.p)
        addend = self
        while k > 0:
            if k & 1:
                result = result + addend
            addend = addend + addend
            k >>= 1
        return result

class EndomorphismConstants:
    beta_i: int
    beta: FieldElement
    lambda_i: int
    neg_lambda: Scalar
    lambda_val: Scalar

    # precomputed decomposition constants
    b1: Scalar
    b2: Scalar
    g1: Scalar
    g2: Scalar

    def __init__(self, p, n, beta_i, beta, lambda_i, lambda_val, b1, b2, g1, g2):
        self.beta_i = beta_i
        self.beta = FieldElement(beta, p)
        self.lambda_i = lambda_i
        self.lambda_val = Scalar(lambda_val, n)
        self.b1 = Scalar(b1, n)
        self.b2 = Scalar(b2, n)
        self.g1 = Scalar(g1, n)
        self.g2 = Scalar(g2, n)

    @classmethod
    def from_params(cls, E:EllipticCurve, p:int, n:int, G:Point):
        (beta_i, beta_val, lambda_i, lambda_val, decompose_params) = _glv_check(E, p, n, G)
        b1,b2,g1,g2 = decompose_params
        return EndomorphismConstants(p, n, beta_i, beta_val, lambda_i, lambda_val, b1, b2, g1, g2)

class Curve256GLV:
    a: int
    b: int
    p: int
    n: int
    G: Point
    glv: EndomorphismConstants

    def __init__(self, a, b, p, n, G, glv:EndomorphismConstants):
        self.a = a
        self.b = b
        self.p = p
        self.n = n
        self.G = G
        self.glv = glv

    def print(self):
        print("a", self.a)
        print("b", self.b)
        print("p", self.p)
        print("n", self.n)
        print("G", self.G)
        print("beta", self.glv.beta)
        print("lambda_val", self.glv.lambda_val)
        print("neg_b1", self.glv.b1)
        print("neg_b2", self.glv.b2)
        print("g1", self.glv.g1)
        print("g2", self.glv.g2)

    def __str__(self):
        return f"Elliptic Curve defined by y^2 = x^3 + {self.b} over Finite Field of size {hex(self.p)}"

    @classmethod
    def from_params(cls, p:int, b:int) -> 'Curve256GLV':
        b = int(b)
        F = FiniteField(p)
        E = EllipticCurve([F(0), F(b)])
        n = E.order()
        G = find_generator(b, p, E)

        glv = EndomorphismConstants.from_params(E, p, n, G)
        G = Point(G[0], G[1], 0, b, p, False)
        return cls(0, b, p, n, G, glv)

    def apply_endomorphism(self, point):
        """Apply the endomorphism φ(P) = (β·x, y)."""
        if point.is_infinity:
            return Point.infinity(self.a, self.b, self.p)
        beta_x = (int(self.glv.beta) * point.x) % self.p
        return Point(beta_x, point.y, self.a, self.b, self.p, False)

    def decompose_scalar(self, k):
        return _glv_decompose(self.n, k, _glv_shift_count(self.n), self.glv.g1, self.glv.g2, self.glv.b1, self.glv.b2, self.glv.lambda_val)

    def decomposition_efficiency(self):
        return _glv_decompose_efficiency(self, self.p, self.n, self.glv.beta, self.glv.lambda_val)

    def scalar_mul_glv(self, point, k):
        """Perform scalar multiplication using GLV decomposition."""
        if k == 0 or point.is_infinity:
            return Point.infinity(self.a, self.b, self.p)

        # Decompose scalar k into k1 and k2
        k1, k2 = self.decompose_scalar(k)

        # Apply the endomorphism to get phi(P)
        phi_p = self.apply_endomorphism(point)

        # Handle negative k1 and k2
        if k1 < 0:
            k1 = -k1
            point = Point(point.x, (-point.y) % self.p, self.a, self.b, self.p, False)
        if k2 < 0:
            k2 = -k2
            phi_p = Point(phi_p.x, (-phi_p.y) % self.p, self.a, self.b, self.p, False)

        # Perform multi-scalar multiplication using interleaving method
        return self.simultaneous_scalar_mul(point, k1, phi_p, k2)

    def simultaneous_scalar_mul(self, p1, k1, p2, k2):
        """Perform simultaneous scalar multiplication k1·P1 + k2·P2."""
        max_bits = max(k1.bit_length(), k2.bit_length())
        k1_bin = bin(k1)[2:].zfill(max_bits)
        k2_bin = bin(k2)[2:].zfill(max_bits)
        p1_plus_p2 = p1 + p2
        result = Point.infinity(self.a, self.b, self.p)
        # Process bits from left to right (most to least significant)
        for i in range(len(k1_bin)):
            result = result + result  # Double
            if k1_bin[i] == '1' and k2_bin[i] == '1':
                result = result + p1_plus_p2
            elif k1_bin[i] == '1':
                result = result + p1
            elif k2_bin[i] == '1':
                result = result + p2
        return result

# Helper function to demonstrate usage
def demonstrate_glv(curve:Curve256GLV):
    """Demonstrate GLV decomposition by comparing with standard scalar multiplication."""
    k = random.randint(1, curve.n - 1)
    k1, k2 = curve.decompose_scalar(k)

    # Verify that k1 + k2 * lambda == k (mod n)
    recomposed = (k1 + (k2 * int(curve.glv.lambda_val)) % curve.n) % curve.n
    assert recomposed == k

    # Compute k*G using standard scalar multiplication
    result_standard = curve.G.scalar_mul(k)
    result_glv = curve.scalar_mul_glv(curve.G, k)
    assert result_standard == result_glv

def test_group_law(curve:Curve256GLV):
    """Test elliptic curve group law properties."""
    G = curve.G

    # Test point addition with identity
    inf = Point.infinity(curve.a, curve.b, curve.p)
    result = G + inf
    assert result == G

    # Test point addition with inverse
    G_neg = Point(G.x, (-G.y) % curve.p, G.a, G.b, G.p, False)
    result = G + G_neg
    assert result.is_infinity is True

    # Test associativity: (G + G) + G = G + (G + G)
    G2 = G + G
    left = G2 + G
    right = G + G2
    assert left == right

    # Test commutativity: G + G2 = G2 + G
    left = G + G2
    right = G2 + G
    assert left == right

def test_scalar_mul_properties(curve:Curve256GLV):
    """Test properties of scalar multiplication."""
    G = curve.G

    # Test 1: 0*G = O
    result = curve.scalar_mul_glv(G, 0)
    assert result.is_infinity is True

    # Test 2: 1*G = G
    result = curve.scalar_mul_glv(G, 1)
    assert result == G

    # Test 3: n*G = O (where n is the curve order)
    result = curve.scalar_mul_glv(G, curve.n)
    assert result.is_infinity is True

    # Test 4: (n-1)*G + G = O
    G_n_minus_1 = curve.scalar_mul_glv(G, curve.n - 1)
    result = G_n_minus_1 + G
    assert result.is_infinity is True

    # Test 5: a*G + b*G = (a+b)*G
    a = random.randint(1, curve.n)
    b = random.randint(1, curve.n)
    left = curve.scalar_mul_glv(G, a) + curve.scalar_mul_glv(G, b)
    right = curve.scalar_mul_glv(G, (a + b) % curve.n)
    assert left == right

    # Test 6: a*(b*G) = (a*b)*G
    a = random.randint(1, curve.n)
    b = random.randint(1, curve.n)
    left = curve.scalar_mul_glv(curve.scalar_mul_glv(G, b), a)
    right = curve.scalar_mul_glv(G, (a * b) % curve.n)
    assert left == right

def test_glv_edge_cases(curve:Curve256GLV):
    """Test edge cases specific to GLV decomposition."""
    G = curve.G

    #print("\n--- Testing GLV Edge Cases ---")

    # Test 1: Very small scalar
    k = 3
    standard = curve.G.scalar_mul(k)
    glv = curve.scalar_mul_glv(G, k)
    assert standard == glv

    # Test 2: Scalar near n/2
    k = curve.n // 2
    standard = curve.G.scalar_mul(k)
    glv = curve.scalar_mul_glv(G, k)
    assert standard == glv

    # Test 3: Scalar near n
    k = curve.n - 2
    standard = curve.G.scalar_mul(k)
    glv = curve.scalar_mul_glv(G, k)
    assert standard == glv

    # Test 4: Random point (not just G)
    rand_scalar = random.randint(1, 10000)
    P = curve.G.scalar_mul(rand_scalar)  # Create a random point
    k = random.randint(1, curve.n - 1)

    standard = P.scalar_mul(k)
    glv = curve.scalar_mul_glv(P, k)
    assert standard == glv

    # Test 5: Boundary case decomposition
    k = curve.n - 1
    k1, k2 = curve.decompose_scalar(k)
    recomposed = (k1 + (k2 * int(curve.glv.lambda_val)) % curve.n) % curve.n
    assert recomposed == k

def test_endomorphism_properties(curve:Curve256GLV):
    """Test properties of the curve endomorphism."""
    G = curve.G

    # Test 1: φ(P) is on the curve
    phi_G = curve.apply_endomorphism(G)
    assert phi_G.is_on_curve() is True

    # Test 2: φ(P+Q) = φ(P) + φ(Q)
    P = G
    Q = G.scalar_mul(2)  # 2*G

    left = curve.apply_endomorphism(P + Q)
    right = curve.apply_endomorphism(P) + curve.apply_endomorphism(Q)
    assert left == right

    # Test 3: φ(P) = λ*P in the group
    # This requires comparing scalar multiplication vs endomorphism
    phi_G = curve.apply_endomorphism(G)
    lambda_G = G.scalar_mul(int(curve.glv.lambda_val))
    assert phi_G == lambda_G

    # Test 4: φ³(P) = P (since λ³ ≡ 1 mod n)
    phi_1 = curve.apply_endomorphism(G)
    phi_2 = curve.apply_endomorphism(phi_1)
    phi_3 = curve.apply_endomorphism(phi_2)
    assert phi_3 == G

def test_performance_comparison(curve:Curve256GLV):
    """Compare performance of standard scalar multiplication vs GLV method."""
    G = curve.G

    total_standard = 0
    total_glv = 0
    total_score = 0
    max_count = 10
    count = 0

    for _ in range(max_count):
        # Use a large scalar for meaningful comparison
        k = random.randint(curve.n // 2, curve.n - 1)

        # Standard scalar multiplication
        start_time = time.time()
        result_standard = G.scalar_mul(k)
        standard_time = time.time() - start_time
        total_standard += standard_time

        # GLV scalar multiplication
        start_time = time.time()
        curve.decompose_scalar(k)
        result_glv = curve.scalar_mul_glv(G, k)
        glv_time = time.time() - start_time
        total_glv += glv_time
        assert result_standard == result_glv

        total_score += _glv_score(curve.decompose_scalar(k), curve.n)

        count += 1

    return (total_standard / total_glv), (total_score/count)

def test_curve(p, b, **extra) -> tuple[Curve256GLV,float]:
    curve = Curve256GLV.from_params(p=p, b=b, **extra)
    demonstrate_glv(curve)
    test_group_law(curve)
    test_scalar_mul_properties(curve)
    test_glv_edge_cases(curve)
    test_endomorphism_properties(curve)
    scores = test_performance_comparison(curve)
    return curve, scores

def main(*args):
    if len(args) == 2:
        p, b = args
        curve, scores = test_curve(int(p), int(b))
    elif len(args) == 4:
        p2,m2p,mx,b = args
        p = 2**int(p2) - 2**int(m2p) - int(mx)
        curve, scores = test_curve(p, int(b))
    else:
        print("Error! unknown args", args)
        return 1
    print(curve)
    return 0

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
