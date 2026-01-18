/-
  Eisenstein.lean - Eisenstein integers Z[ω] where ω² + ω + 1 = 0
-/

import EllipticPrimeOrder.Basic

/-! # Eisenstein Integers

The Eisenstein integers Z[ω] are complex numbers of the form a + bω where
ω = (-1 + √(-3))/2 is a primitive cube root of unity satisfying ω² + ω + 1 = 0.

For our purposes, we use Eisenstein coordinates (c, d) related to Cornacchia
coordinates (a, b) via:
- c = a + b
- d = 2b

where p = a² + 3b² (Cornacchia) corresponds to the Eisenstein norm c² - cd + d².
-/

/-- Eisenstein integer represented as re + im·ω where ω² + ω + 1 = 0 -/
structure Eisenstein where
  re : ℤ  -- coefficient of 1
  im : ℤ  -- coefficient of ω
  deriving Repr, DecidableEq

namespace Eisenstein

/-- The Eisenstein norm: N(a + bω) = a² - ab + b² -/
def norm (x : Eisenstein) : ℤ := x.re^2 - x.re * x.im + x.im^2

/-- Addition of Eisenstein integers -/
def add (x y : Eisenstein) : Eisenstein :=
  ⟨x.re + y.re, x.im + y.im⟩

/-- Multiplication of Eisenstein integers -/
def mul (x y : Eisenstein) : Eisenstein :=
  -- (a + bω)(c + dω) = ac + (ad + bc)ω + bd·ω²
  -- = ac + (ad + bc)ω + bd(-1 - ω)  [since ω² = -1 - ω]
  -- = (ac - bd) + (ad + bc - bd)ω
  ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re - x.im * y.im⟩

/-- Zero Eisenstein integer -/
def zero : Eisenstein := ⟨0, 0⟩

/-- Unit Eisenstein integer -/
def one : Eisenstein := ⟨1, 0⟩

/-- Negation -/
def neg (x : Eisenstein) : Eisenstein := ⟨-x.re, -x.im⟩

/-- Conjugate: (a + bω)* = a + b·ω² = a + b(-1-ω) = (a - b) - bω
    This is the Galois conjugate in Z[ω], replacing ω with ω² = -1 - ω -/
def conj (x : Eisenstein) : Eisenstein := ⟨x.re - x.im, -x.im⟩

instance : Add Eisenstein := ⟨add⟩
instance : Mul Eisenstein := ⟨mul⟩
instance : Neg Eisenstein := ⟨neg⟩
instance : Zero Eisenstein := ⟨zero⟩
instance : One Eisenstein := ⟨one⟩

/-- The norm is multiplicative -/
theorem norm_mul (x y : Eisenstein) : norm (x * y) = norm x * norm y := by
  show norm (mul x y) = norm x * norm y
  simp only [norm, mul]
  ring

/-- Norm of conjugate equals norm -/
theorem norm_conj (x : Eisenstein) : norm (conj x) = norm x := by
  simp only [norm, conj]
  ring

/-- x * conj x = norm x (as the re component, im is 0) -/
theorem mul_conj_eq_norm (x : Eisenstein) : (x * conj x).re = norm x ∧ (x * conj x).im = 0 := by
  show (mul x (conj x)).re = norm x ∧ (mul x (conj x)).im = 0
  simp only [mul, conj, norm]
  constructor <;> ring

end Eisenstein

/-! ## Coordinate Transformation

The relationship between Cornacchia coordinates (a, b) and Eisenstein coordinates (c, d):
- c = a + b
- d = 2b
- a = c - d/2 = (2c - d)/2
- b = d/2

The Cornacchia form a² + 3b² equals the Eisenstein norm c² - cd + d² when c = a + b, d = 2b.
-/

/-- Convert Cornacchia coordinates to Eisenstein coordinates -/
def cornacchiaToEisenstein (a b : ℤ) : ℤ × ℤ := (a + b, 2 * b)

/-- The Eisenstein norm of (c, d) where c = a + b, d = 2b equals a² + 3b² -/
theorem eisenstein_norm_eq_cornacchia (a b : ℤ) :
    let (c, d) := cornacchiaToEisenstein a b
    c^2 - c*d + d^2 = a^2 + 3*b^2 := by
  simp only [cornacchiaToEisenstein]
  ring
