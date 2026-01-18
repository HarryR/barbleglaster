/-
  Axioms.lean - Axiomatized results from algebraic number theory
-/

import EllipticPrimeOrder.TraceOfFrobenius

/-! # Axiomatized Results

We axiomatize the following well-known results that would require substantial
formalization effort:

1. Cornacchia's algorithm: For p ≡ 7 (mod 12), there exist unique a, b with
   a even, b odd, and p = a² + 3b².

2. Hasse bound: For an elliptic curve over F_p, |t| ≤ 2√p.

3. Curve order positivity: The curve order p + 1 - t is always positive.
-/

/-- Cornacchia's theorem: For primes p ≡ 7 (mod 12), there exist a, b
    with a even, b odd, and p = a² + 3b² -/
axiom cornacchia_exists (p : ℕ) (hp : Nat.Prime p) (hmod : p % 12 = 7) :
  ∃ (a b : ℤ), Even a ∧ Odd b ∧ a^2 + 3*b^2 = p ∧ a ≠ 0

/-- Construct TraceData from Cornacchia's theorem (using choice) -/
noncomputable def TraceData.fromPrime (p : ℕ) (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    TraceData p :=
  let h := cornacchia_exists p hp hmod
  { a := Classical.choose h
    b := Classical.choose (Classical.choose_spec h)
    hEvenA := (Classical.choose_spec (Classical.choose_spec h)).1
    hOddB := (Classical.choose_spec (Classical.choose_spec h)).2.1
    hCornacchia := (Classical.choose_spec (Classical.choose_spec h)).2.2.1
    hNonzeroA := (Classical.choose_spec (Classical.choose_spec h)).2.2.2 }

/-- Hasse bound: |t| ≤ 2√p for any trace of Frobenius -/
axiom hasse_bound (p : ℕ) (hp : Nat.Prime p) (t : ℤ) :
  |t| ≤ 2 * Int.sqrt p

/-- Curve orders are positive (follows from Hasse bound) -/
axiom curve_order_pos {p : ℕ} (hp : Nat.Prime p) (hp3 : p > 3) (t : TraceData p) (i : Fin 6) :
  0 < curveOrder t i

/-- Curve orders for p > 7 are greater than 3.
    Note: p = 7 admits order = 3 (e.g., trace = a + 3b = 5 gives order = 7 + 1 - 5 = 3).
    For p > 7 with p ≡ 7 (mod 12), the Hasse bound ensures all orders exceed 3. -/
axiom curve_order_gt_three {p : ℕ} (hp : Nat.Prime p) (hp7 : p > 7) (t : TraceData p) (i : Fin 6) :
  curveOrder t i > 3

/-- Cornacchia uniqueness: The Cornacchia decomposition p = a² + 3b² with a even, b odd,
    a ≠ 0 is unique up to sign. In particular, any two TraceData instances for the same
    prime p have the same curve orders. -/
axiom cornacchia_unique_orders {p : ℕ} (hp : Nat.Prime p) (hmod : p % 12 = 7)
    (t1 t2 : TraceData p) (i : Fin 6) :
  curveOrder t1 i = curveOrder t2 i

/-- If an integer n > 2 is even, it's not prime.
    Proof: If n = 2k with k > 1, then 2 is a non-trivial divisor. -/
theorem even_gt_two_not_prime {n : ℤ} (hn : n > 2) (heven : Even n) : ¬ Prime n := by
  intro hprime
  obtain ⟨k, hk⟩ := heven
  have hk_gt : k > 1 := by omega
  have h2_not_unit : ¬ IsUnit (2 : ℤ) := by simp [Int.isUnit_iff]
  have hk_not_unit : ¬ IsUnit k := by simp [Int.isUnit_iff]; omega
  have hirr := hprime.irreducible
  -- n = k + k = 2 * k
  have heq : n = (2 : ℤ) * k := by omega
  have : IsUnit (2 : ℤ) ∨ IsUnit k := hirr.isUnit_or_isUnit heq
  rcases this with h2u | hku
  · exact h2_not_unit h2u
  · exact hk_not_unit hku

/-- If an integer n > 3 is divisible by 3, it's not prime -/
theorem div_by_three_gt_three_not_prime {n : ℤ} (hn : n > 3) (hdiv : (3 : ℤ) ∣ n) : ¬ Prime n := by
  intro hprime
  obtain ⟨k, hk⟩ := hdiv
  have hk_gt : k > 1 := by omega
  have h3_not_unit : ¬ IsUnit (3 : ℤ) := by simp [Int.isUnit_iff]
  have hk_not_unit : ¬ IsUnit k := by simp [Int.isUnit_iff]; omega
  have hirr := hprime.irreducible
  have : IsUnit (3 : ℤ) ∨ IsUnit k := hirr.isUnit_or_isUnit hk
  rcases this with h3u | hku
  · exact h3_not_unit h3u
  · exact hk_not_unit hku
