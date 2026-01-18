/-
  TraceOfFrobenius.lean - Trace data structure and coefficient formulas
-/

import EllipticPrimeOrder.Eisenstein

/-! # Trace of Frobenius

For a prime p ≡ 7 (mod 12), the six sextic twists of j-invariant 0 curves have
traces of Frobenius that can be computed from Eisenstein coordinates (c, d).

The trace for curve index i is: t_i = coef_c[i] * c + coef_d[i] * d

where the coefficients are determined by 6th roots of unity structure.
-/

/-- Trace data encapsulating the Cornacchia decomposition p = a² + 3b²
    and derived Eisenstein coordinates c = a + b, d = 2b -/
structure TraceData (p : ℕ) where
  /-- Cornacchia coordinate a (even) -/
  a : ℤ
  /-- Cornacchia coordinate b (odd) -/
  b : ℤ
  /-- a is even -/
  hEvenA : Even a
  /-- b is odd -/
  hOddB : Odd b
  /-- Cornacchia representation: p = a² + 3b² -/
  hCornacchia : a^2 + 3*b^2 = p
  /-- a is nonzero (follows from p being prime > 3) -/
  hNonzeroA : a ≠ 0

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-- Eisenstein coordinate c = a + b -/
def c : ℤ := t.a + t.b

/-- Eisenstein coordinate d = 2b -/
def d : ℤ := 2 * t.b

/-- The Eisenstein norm c² - cd + d² equals p -/
theorem eisenstein_norm_eq_p : t.c^2 - t.c * t.d + t.d^2 = p := by
  simp only [c, d]
  have h := t.hCornacchia
  linarith [sq_nonneg t.a, sq_nonneg t.b, sq_nonneg (t.a + t.b), sq_nonneg (t.a - t.b)]

/-- d is even (since d = 2b) -/
theorem d_even : Even t.d := by
  simp only [d]
  exact even_two_mul t.b

end TraceData

/-! ## Trace Coefficients

The trace coefficients for each of the 6 sextic twists.
trace_i = coef_c[i] * c + coef_d[i] * d

These coefficients arise from the 6th roots of unity expansion:
t_i = α·ζ₆^i + β·ζ₆^(-i)
-/

/-- Trace coefficients (coef_c, coef_d) for each curve index i ∈ {0,1,2,3,4,5}
    The trace is computed as: t_i = coef_c * c + coef_d * d -/
def traceCoefficients : Fin 6 → ℤ × ℤ
  | 0 => (-2, 1)   -- t_0 = -2c + d = -2a
  | 1 => (-1, -1)  -- t_1 = -c - d = -a - 3b
  | 2 => (1, -2)   -- t_2 = c - 2d = a - 3b
  | 3 => (2, -1)   -- t_3 = 2c - d = 2a
  | 4 => (1, 1)    -- t_4 = c + d = a + 3b
  | 5 => (-1, 2)   -- t_5 = -c + 2d = -a + 3b

/-- Compute the trace of Frobenius for curve index i -/
def traceOfFrobenius {p : ℕ} (t : TraceData p) (i : Fin 6) : ℤ :=
  let (coef_c, coef_d) := traceCoefficients i
  coef_c * t.c + coef_d * t.d

/-- Compute the curve order for curve index i: order = p + 1 - trace -/
def curveOrder {p : ℕ} (t : TraceData p) (i : Fin 6) : ℤ :=
  (p : ℤ) + 1 - traceOfFrobenius t i

/-! ## Explicit Trace Formulas

Simplified expressions for traces in terms of Cornacchia coordinates (a, b).
-/

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-- t_0 = -2c + d = -2(a+b) + 2b = -2a -/
theorem trace_0_eq : traceOfFrobenius t 0 = -2 * t.a := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- t_1 = -c - d = -(a+b) - 2b = -a - 3b -/
theorem trace_1_eq : traceOfFrobenius t 1 = -t.a - 3*t.b := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- t_2 = c - 2d = (a+b) - 4b = a - 3b -/
theorem trace_2_eq : traceOfFrobenius t 2 = t.a - 3*t.b := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- t_3 = 2c - d = 2(a+b) - 2b = 2a -/
theorem trace_3_eq : traceOfFrobenius t 3 = 2 * t.a := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- t_4 = c + d = (a+b) + 2b = a + 3b -/
theorem trace_4_eq : traceOfFrobenius t 4 = t.a + 3*t.b := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- t_5 = -c + 2d = -(a+b) + 4b = -a + 3b -/
theorem trace_5_eq : traceOfFrobenius t 5 = -t.a + 3*t.b := by
  simp only [traceOfFrobenius, traceCoefficients, c, d]
  ring

/-- Order at index 0: p + 1 + 2a -/
theorem order_0_eq : curveOrder t 0 = (p : ℤ) + 1 + 2*t.a := by
  simp only [curveOrder, trace_0_eq]
  ring

/-- Order at index 1: p + 1 + a + 3b -/
theorem order_1_eq : curveOrder t 1 = (p : ℤ) + 1 + t.a + 3*t.b := by
  simp only [curveOrder, trace_1_eq]
  ring

/-- Order at index 2: p + 1 - a + 3b -/
theorem order_2_eq : curveOrder t 2 = (p : ℤ) + 1 - t.a + 3*t.b := by
  simp only [curveOrder, trace_2_eq]
  ring

/-- Order at index 3: p + 1 - 2a -/
theorem order_3_eq : curveOrder t 3 = (p : ℤ) + 1 - 2*t.a := by
  simp only [curveOrder, trace_3_eq]

/-- Order at index 4: p + 1 - a - 3b -/
theorem order_4_eq : curveOrder t 4 = (p : ℤ) + 1 - t.a - 3*t.b := by
  simp only [curveOrder, trace_4_eq]
  ring

/-- Order at index 5: p + 1 + a - 3b -/
theorem order_5_eq : curveOrder t 5 = (p : ℤ) + 1 + t.a - 3*t.b := by
  simp only [curveOrder, trace_5_eq]
  ring

end TraceData
