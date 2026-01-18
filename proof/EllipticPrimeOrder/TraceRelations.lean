/-
  TraceRelations.lean - Deeper structural relationships among traces and orders

  This module formalizes additional relationships implicit in the trace structure:
  1. Traces sum to zero (conservation from 6th roots of unity)
  2. Quadratic twist pairing: t_i + t_{i+3} = 0
  3. Order sum equals 6(p+1)
  4. Product formulas for complementary orders
  5. Characterization of when orders equal specific values
-/

import EllipticPrimeOrder.TraceOfFrobenius

/-! # Trace Sum Relations

The six traces of Frobenius for j-invariant 0 sextic twists satisfy a fundamental
conservation law: their sum is zero. This follows from the 6th roots of unity
structure, where ∑_{i=0}^{5} ζ₆^i = 0.
-/

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-! ## Quadratic Twist Pairing

Indices i and (i+3) mod 6 form quadratic twist pairs. Their traces are negatives:
  t_i + t_{i+3} = 0

This means their orders sum to 2(p+1):
  order_i + order_{i+3} = 2(p+1)
-/

/-- Traces at indices 0 and 3 are negatives: t₀ + t₃ = 0 -/
theorem trace_0_3_neg : traceOfFrobenius t 0 + traceOfFrobenius t 3 = 0 := by
  rw [trace_0_eq, trace_3_eq]
  ring

/-- Traces at indices 1 and 4 are negatives: t₁ + t₄ = 0 -/
theorem trace_1_4_neg : traceOfFrobenius t 1 + traceOfFrobenius t 4 = 0 := by
  rw [trace_1_eq, trace_4_eq]
  ring

/-- Traces at indices 2 and 5 are negatives: t₂ + t₅ = 0 -/
theorem trace_2_5_neg : traceOfFrobenius t 2 + traceOfFrobenius t 5 = 0 := by
  rw [trace_2_eq, trace_5_eq]
  ring

/-- Sum of all six traces is zero -/
theorem trace_sum_zero :
    traceOfFrobenius t 0 + traceOfFrobenius t 1 + traceOfFrobenius t 2 +
    traceOfFrobenius t 3 + traceOfFrobenius t 4 + traceOfFrobenius t 5 = 0 := by
  rw [trace_0_eq, trace_1_eq, trace_2_eq, trace_3_eq, trace_4_eq, trace_5_eq]
  ring

/-! ## Order Sum Relations -/

/-- Quadratic twist pair: orders at 0 and 3 sum to 2(p+1) -/
theorem order_0_3_sum : curveOrder t 0 + curveOrder t 3 = 2 * ((p : ℤ) + 1) := by
  simp only [curveOrder]
  have h := trace_0_3_neg t
  linarith

/-- Quadratic twist pair: orders at 1 and 4 sum to 2(p+1) -/
theorem order_1_4_sum : curveOrder t 1 + curveOrder t 4 = 2 * ((p : ℤ) + 1) := by
  simp only [curveOrder]
  have h := trace_1_4_neg t
  linarith

/-- Quadratic twist pair: orders at 2 and 5 sum to 2(p+1) -/
theorem order_2_5_sum : curveOrder t 2 + curveOrder t 5 = 2 * ((p : ℤ) + 1) := by
  simp only [curveOrder]
  have h := trace_2_5_neg t
  linarith

/-- Sum of all six orders equals 6(p+1) -/
theorem order_sum :
    curveOrder t 0 + curveOrder t 1 + curveOrder t 2 +
    curveOrder t 3 + curveOrder t 4 + curveOrder t 5 = 6 * ((p : ℤ) + 1) := by
  have h03 := order_0_3_sum t
  have h14 := order_1_4_sum t
  have h25 := order_2_5_sum t
  linarith

/-! ## Product Formulas

The product of orders in a quadratic twist pair has a beautiful form:
  order_i × order_{i+3} = (p+1)² - t_i²

Since t_i and t_{i+3} are negatives, t_i² = t_{i+3}².
-/

/-- Product of orders at 0 and 3: (p+1+2a)(p+1-2a) = (p+1)² - 4a² -/
theorem order_0_3_product :
    curveOrder t 0 * curveOrder t 3 = ((p : ℤ) + 1)^2 - 4 * t.a^2 := by
  rw [order_0_eq, order_3_eq]
  ring

/-- Product of orders at 1 and 4: (p+1)² - (a+3b)² -/
theorem order_1_4_product :
    curveOrder t 1 * curveOrder t 4 = ((p : ℤ) + 1)^2 - (t.a + 3*t.b)^2 := by
  rw [order_1_eq, order_4_eq]
  ring

/-- Product of orders at 2 and 5: (p+1)² - (a-3b)² -/
theorem order_2_5_product :
    curveOrder t 2 * curveOrder t 5 = ((p : ℤ) + 1)^2 - (t.a - 3*t.b)^2 := by
  rw [order_2_eq, order_5_eq]
  ring

/-! ## Trace Squared Relations

The squares of traces satisfy useful identities.
-/

/-- t₀² = 4a² -/
theorem trace_0_sq : (traceOfFrobenius t 0)^2 = 4 * t.a^2 := by
  rw [trace_0_eq]; ring

/-- t₁² = (a + 3b)² -/
theorem trace_1_sq : (traceOfFrobenius t 1)^2 = (t.a + 3*t.b)^2 := by
  rw [trace_1_eq]; ring

/-- t₂² = (a - 3b)² -/
theorem trace_2_sq : (traceOfFrobenius t 2)^2 = (t.a - 3*t.b)^2 := by
  rw [trace_2_eq]

/-- Sum of squared traces: t₀² + t₃² = 8a², and (t₁² + t₄²) + (t₂² + t₅²) = 4(a+3b)² + 4(a-3b)²
    Total: 8a² + 4(a² + 9b²) + 4(a² + 9b²) - ... Let's compute directly. -/
theorem trace_sq_sum :
    (traceOfFrobenius t 0)^2 + (traceOfFrobenius t 1)^2 + (traceOfFrobenius t 2)^2 +
    (traceOfFrobenius t 3)^2 + (traceOfFrobenius t 4)^2 + (traceOfFrobenius t 5)^2 =
    12 * t.a^2 + 36 * t.b^2 := by
  rw [trace_0_eq, trace_1_eq, trace_2_eq, trace_3_eq, trace_4_eq, trace_5_eq]
  ring

/-- Simplified: sum of squared traces = 12p -/
theorem trace_sq_sum_simplified :
    (traceOfFrobenius t 0)^2 + (traceOfFrobenius t 1)^2 + (traceOfFrobenius t 2)^2 +
    (traceOfFrobenius t 3)^2 + (traceOfFrobenius t 4)^2 + (traceOfFrobenius t 5)^2 =
    12 * (p : ℤ) := by
  have h := trace_sq_sum t
  have hcorn := t.hCornacchia
  -- a² + 3b² = p, so 12a² + 36b² = 12(a² + 3b²) = 12p
  calc (traceOfFrobenius t 0)^2 + (traceOfFrobenius t 1)^2 + (traceOfFrobenius t 2)^2 +
       (traceOfFrobenius t 3)^2 + (traceOfFrobenius t 4)^2 + (traceOfFrobenius t 5)^2
      = 12 * t.a^2 + 36 * t.b^2 := h
    _ = 12 * (t.a^2 + 3 * t.b^2) := by ring
    _ = 12 * (p : ℤ) := by rw [hcorn]

/-! ## Anomalous Curve Detection

A curve is anomalous if its order equals p (trace = 1). We can characterize
when this occurs in terms of Cornacchia coordinates.
-/

/-- Order at index i equals p iff trace equals 1 -/
theorem order_eq_p_iff_trace_one (i : Fin 6) :
    curveOrder t i = (p : ℤ) ↔ traceOfFrobenius t i = 1 := by
  simp only [curveOrder]
  constructor <;> intro h <;> linarith

/-- Order 0 = p iff 2a = -1 (impossible for integer a) -/
theorem order_0_eq_p_iff : curveOrder t 0 = (p : ℤ) ↔ 2 * t.a = -1 := by
  rw [order_eq_p_iff_trace_one, trace_0_eq]
  constructor <;> intro h <;> linarith

/-- Order 0 can never equal p (2a ≠ -1 for integers) -/
theorem order_0_ne_p : curveOrder t 0 ≠ (p : ℤ) := by
  intro h
  have h2a := (order_0_eq_p_iff t).mp h
  -- 2a = -1 has no integer solution (2 doesn't divide -1)
  have hdiv : (2 : ℤ) ∣ 2 * t.a := ⟨t.a, rfl⟩
  rw [h2a] at hdiv
  norm_num at hdiv

/-- Order 3 can never equal p (2a ≠ 1 for even a) -/
theorem order_3_ne_p : curveOrder t 3 ≠ (p : ℤ) := by
  intro h
  have ht := (order_eq_p_iff_trace_one t 3).mp h
  rw [trace_3_eq] at ht
  -- 2a = 1 means a = 1/2, but a is even
  obtain ⟨k, hk⟩ := t.hEvenA
  rw [hk] at ht
  -- 2 * (2k) = 1 → 4k = 1, impossible for integer k
  omega

/-! ## Order Relations to p

The orders relate to p through the Cornacchia form p = a² + 3b².
-/

/-- Order 0 × Order 3 in terms of p -/
theorem order_0_3_product_alt :
    curveOrder t 0 * curveOrder t 3 = ((p : ℤ) + 1)^2 - 4 * ((p : ℤ) - 3 * t.b^2) := by
  have hprod := order_0_3_product t
  have hcorn := t.hCornacchia
  calc curveOrder t 0 * curveOrder t 3
      = ((p : ℤ) + 1)^2 - 4 * t.a^2 := hprod
    _ = ((p : ℤ) + 1)^2 - 4 * ((p : ℤ) - 3 * t.b^2) := by rw [← hcorn]; ring

/-- Simplified: order_0 × order_3 = (p+1)² - 4p + 12b² = (p-1)² + 12b² -/
theorem order_0_3_product_simplified :
    curveOrder t 0 * curveOrder t 3 = ((p : ℤ) - 1)^2 + 12 * t.b^2 := by
  have h := order_0_3_product_alt t
  linarith [sq_nonneg ((p : ℤ) + 1), sq_nonneg ((p : ℤ) - 1)]

end TraceData

/-! # Corollaries About Prime Orders

These results give additional constraints on when orders can be prime.
-/

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-- If order_0 is prime, order_3 divides (p-1)² + 12b² -/
theorem order_3_divides_if_order_0_prime (_h0 : Prime (curveOrder t 0)) :
    curveOrder t 3 ∣ ((p : ℤ) - 1)^2 + 12 * t.b^2 := by
  have hprod := order_0_3_product_simplified t
  exact ⟨curveOrder t 0, by linarith⟩

/-- The quadratic twist pairing gives: if one order is > (p+1),
    then its twist has order < (p+1) -/
theorem twist_order_bound_0 (hq : curveOrder t 0 > (p : ℤ) + 1) :
    curveOrder t 3 < (p : ℤ) + 1 := by
  have hsum := order_0_3_sum t
  linarith

theorem twist_order_bound_1 (hq : curveOrder t 1 > (p : ℤ) + 1) :
    curveOrder t 4 < (p : ℤ) + 1 := by
  have hsum := order_1_4_sum t
  linarith

theorem twist_order_bound_2 (hq : curveOrder t 2 > (p : ℤ) + 1) :
    curveOrder t 5 < (p : ℤ) + 1 := by
  have hsum := order_2_5_sum t
  linarith

end TraceData
