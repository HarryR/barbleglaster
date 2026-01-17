/-
  MainTheorem.lean - Main theorem: at most 2 curve orders can be prime

  IMPORTANT: This must be consistent with the Python implementation in
  lemma/2-eisenstein-mapping-magic.py
-/

import EllipticPrimeOrder.Divisibility
import EllipticPrimeOrder.IndexMapping

/-! # Main Theorem

For a prime p ≡ 7 (mod 12), among the six sextic twists of j-invariant 0 curves,
at most two can have prime order. Specifically, curves E_{g^1} and E_{g^5} are
the only candidates for prime order.

The proof proceeds by showing:
1. Orders at curveOrder indices 0 and 3 are always even → composite
2. Among curveOrder indices {1, 2, 4, 5}, exactly two are divisible by 3 → composite
3. The permutation ensures curves E_{g^1} and E_{g^5} always get non-div-by-3 orders

## Key Insight

The permutation `curveToOrder` maps curve indices to curveOrder indices such that:
- When a ≡ 1 (mod 3): curveOrders {1, 5} are div by 3, curves {1, 5} get {2, 4}
- When a ≡ 2 (mod 3): curveOrders {2, 4} are div by 3, curves {1, 5} get {1, 5}

In both cases, curves E_{g^1} and E_{g^5} receive non-div-by-3 curveOrders!
-/

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-- Orders at curveOrder indices 0 and 3 cannot be prime (they are even and > 2) -/
theorem order_0_not_prime (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    ¬ Prime (curveOrder t 0) := by
  have heven := order_0_even t hmod
  have hgt := curve_order_gt_three hp (prime_mod7_12_gt_three ⟨hp, hmod⟩) t 0
  exact even_gt_two_not_prime (by linarith) heven

theorem order_3_not_prime (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    ¬ Prime (curveOrder t 3) := by
  have heven := order_3_even t hmod
  have hgt := curve_order_gt_three hp (prime_mod7_12_gt_three ⟨hp, hmod⟩) t 3
  exact even_gt_two_not_prime (by linarith) heven

/-- If curveOrder index i has order divisible by 3 and > 3, then order is not prime -/
theorem order_div3_not_prime (hp : Nat.Prime p) (hmod : p % 12 = 7) (i : Fin 6)
    (hdiv : 3 ∣ curveOrder t i) : ¬ Prime (curveOrder t i) := by
  have hgt := curve_order_gt_three hp (prime_mod7_12_gt_three ⟨hp, hmod⟩) t i
  exact div_by_three_gt_three_not_prime hgt hdiv

/-- At most 4 curveOrder indices have composite order -/
theorem at_least_four_composite (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    ¬ Prime (curveOrder t 0) ∧
    ¬ Prime (curveOrder t 3) ∧
    ((¬ Prime (curveOrder t 1) ∧ ¬ Prime (curveOrder t 5)) ∨
     (¬ Prime (curveOrder t 2) ∧ ¬ Prime (curveOrder t 4))) := by
  constructor
  · exact order_0_not_prime t hp hmod
  constructor
  · exact order_3_not_prime t hp hmod
  · cases some_orders_div_3 t hmod with
    | inl h15 =>
      left
      constructor
      · exact order_div3_not_prime t hp hmod 1 h15.1
      · exact order_div3_not_prime t hp hmod 5 h15.2
    | inr h24 =>
      right
      constructor
      · exact order_div3_not_prime t hp hmod 2 h24.1
      · exact order_div3_not_prime t hp hmod 4 h24.2

/-- Main theorem: At most 2 of the 6 curveOrders can be prime -/
theorem at_most_two_prime_orders (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    ∃ (S : Finset (Fin 6)), S.card ≤ 2 ∧
      ∀ i : Fin 6, Prime (curveOrder t i) → i ∈ S := by
  have h := at_least_four_composite t hp hmod
  cases h.2.2 with
  | inl h15 =>
    -- curveOrders {1, 5} are composite, potentially prime is {2, 4}
    use {2, 4}
    constructor
    · simp
    · intro i hi
      fin_cases i <;> simp at * <;> try exact absurd hi h.1
      · exact absurd hi h15.1
      · exact absurd hi h.2.1
      · exact absurd hi h15.2
  | inr h24 =>
    -- curveOrders {2, 4} are composite, potentially prime is {1, 5}
    use {1, 5}
    constructor
    · simp
    · intro i hi
      fin_cases i <;> simp at * <;> try exact absurd hi h.1
      · exact absurd hi h24.1
      · exact absurd hi h.2.1
      · exact absurd hi h24.2

end TraceData

/-! ## The Main Result: Curves E_{g^1} and E_{g^5} Can Have Prime Order

This section proves that curves E_{g^1} and E_{g^5} always receive curveOrders
that are not divisible by 3, making them the only candidates for prime order.
-/

open IndexMapping

/-- Curves E_{g^1} and E_{g^5} always get non-div-by-3 curveOrders.

    When u0 = 0 (a ≡ 1 mod 3):
    - curveOrder indices {1, 5} ARE divisible by 3
    - permIdx < 2, so curves {1, 5} get curveOrders from {2, 4}
    - curveOrders {2, 4} are NOT divisible by 3 ✓

    When u0 = 1 (a ≡ 2 mod 3):
    - curveOrder indices {2, 4} ARE divisible by 3
    - permIdx ≥ 2, so curves {1, 5} get curveOrders from {1, 5}
    - curveOrders {1, 5} are NOT divisible by 3 ✓ -/
theorem curves_1_5_get_nondiv3_orders {p : ℕ} (_hp : Nat.Prime p) (hmod : p % 12 = 7)
    (t : TraceData p) (g : ZMod p) :
    let orderIdx1 := curveToOrder (permIdx t g) 1
    let orderIdx5 := curveToOrder (permIdx t g) 5
    ¬(3 ∣ curveOrder t orderIdx1) ∧ ¬(3 ∣ curveOrder t orderIdx5) := by
  have hdiv := TraceData.exactly_two_div_3 t hmod
  cases hdiv with
  | inl h15 =>
    -- a ≡ 1 (mod 3): curveOrders {1, 5} div by 3, {2, 4} NOT div by 3
    have hu0 : u0 t = 0 := by
      rw [u0_zero_iff t hmod]
      have h1_mod := TraceData.order1_mod3_eq t hmod
      have h1_div : curveOrder t 1 % 3 = 0 := Int.emod_eq_zero_of_dvd h15.1
      rw [h1_mod] at h1_div
      have hne := TraceData.a_not_div_3 t hmod
      have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
        have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
        have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
        omega
      cases h_range with
      | inl h1 => exact h1
      | inr h2 =>
        exfalso
        have h_calc : (2 + t.a) % 3 = 1 := by
          have h_add : (2 + t.a) % 3 = (2 % 3 + t.a % 3) % 3 := Int.add_emod _ _ _
          rw [h_add, h2]; norm_num
        rw [h_calc] at h1_div
        norm_num at h1_div
    have hperm_lt2 := u0_zero_permIdx_lt2 t g hu0
    have hget := idx_lt2_curves_get_24 (permIdx t g) hperm_lt2
    -- h15.2.2 gives: ¬(3 ∣ curveOrder t 2) ∧ ¬(3 ∣ curveOrder t 4)
    have h24_not := h15.2.2
    constructor
    · -- Curve E_{g^1} gets a curveOrder from {2, 4}, which are NOT div by 3
      have hmem := hget.1
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      cases hmem with
      | inl h => rw [h]; exact h24_not.1
      | inr h => rw [h]; exact h24_not.2
    · -- Curve E_{g^5} gets a curveOrder from {2, 4}, which are NOT div by 3
      have hmem := hget.2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      cases hmem with
      | inl h => rw [h]; exact h24_not.1
      | inr h => rw [h]; exact h24_not.2
  | inr h24 =>
    -- a ≡ 2 (mod 3): curveOrders {2, 4} div by 3, {1, 5} NOT div by 3
    have hu0 : u0 t = 1 := by
      rw [u0_one_iff t hmod]
      have h2_mod := TraceData.order2_mod3_eq t hmod
      have h2_div : curveOrder t 2 % 3 = 0 := Int.emod_eq_zero_of_dvd h24.1
      rw [h2_mod] at h2_div
      have hne := TraceData.a_not_div_3 t hmod
      have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
        have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
        have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
        omega
      cases h_range with
      | inl h1 =>
        exfalso
        have h_calc : (2 - t.a) % 3 = 1 := by
          have h_sub : (2 - t.a) % 3 = (2 % 3 - t.a % 3) % 3 := Int.sub_emod _ _ _
          rw [h_sub, h1]; norm_num
        rw [h_calc] at h2_div
        norm_num at h2_div
      | inr h2 => exact h2
    have hperm_ge2 := u0_one_permIdx_ge2 t g hu0
    have hget := idx_ge2_curves_get_15 (permIdx t g) hperm_ge2
    -- h24.2.2 gives: ¬(3 ∣ curveOrder t 1) ∧ ¬(3 ∣ curveOrder t 5)
    have h15_not := h24.2.2
    constructor
    · -- Curve E_{g^1} gets a curveOrder from {1, 5}, which are NOT div by 3
      have hmem := hget.1
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      cases hmem with
      | inl h => rw [h]; exact h15_not.1
      | inr h => rw [h]; exact h15_not.2
    · -- Curve E_{g^5} gets a curveOrder from {1, 5}, which are NOT div by 3
      have hmem := hget.2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      cases hmem with
      | inl h => rw [h]; exact h15_not.1
      | inr h => rw [h]; exact h15_not.2

/-- If a curve E_{g^j} has prime order, then j ∈ {1, 5}.

    This follows from:
    1. Curves E_{g^0} and E_{g^3} always have even orders (composite)
    2. Among curves {1, 2, 4, 5}, exactly one pair has div-by-3 orders
    3. The permutation ensures curves {1, 5} get the non-div-by-3 orders -/
theorem prime_order_implies_curve_1_or_5 {p : ℕ} (hp : Nat.Prime p) (hmod : p % 12 = 7)
    (t : TraceData p) (g : ZMod p) (curveIdx : Fin 6)
    (hPrime : Prime (curveOrder t (curveToOrder (permIdx t g) curveIdx))) :
    curveIdx = 1 ∨ curveIdx = 5 := by
  -- Helper: curveToOrder _ 0 is always 0 or 3 (both even orders)
  have hcurve0_even : ∀ idx : Fin 4, Even (curveOrder t (curveToOrder idx 0)) := by
    intro idx
    fin_cases idx <;> simp [curveToOrder] <;>
      (first | exact TraceData.order_3_even t hmod | exact TraceData.order_0_even t hmod)
  -- Helper: curveToOrder _ 3 is always 0 or 3 (both even orders)
  have hcurve3_even : ∀ idx : Fin 4, Even (curveOrder t (curveToOrder idx 3)) := by
    intro idx
    fin_cases idx <;> simp [curveToOrder] <;>
      (first | exact TraceData.order_0_even t hmod | exact TraceData.order_3_even t hmod)
  -- Helper: all curve orders are > 3
  have hgt_helper : ∀ i : Fin 6, curveOrder t i > 3 := by
    intro i
    exact curve_order_gt_three hp (prime_mod7_12_gt_three ⟨hp, hmod⟩) t i
  -- First, curveIdx cannot be 0 or 3 (even orders)
  have h03 : curveIdx ≠ 0 ∧ curveIdx ≠ 3 := by
    constructor
    · intro h0
      rw [h0] at hPrime
      have heven := hcurve0_even (permIdx t g)
      have hgt := hgt_helper (curveToOrder (permIdx t g) 0)
      exact even_gt_two_not_prime (by linarith) heven hPrime
    · intro h3
      rw [h3] at hPrime
      have heven := hcurve3_even (permIdx t g)
      have hgt := hgt_helper (curveToOrder (permIdx t g) 3)
      exact even_gt_two_not_prime (by linarith) heven hPrime
  -- Now use the div-by-3 analysis
  have hndiv := curves_1_5_get_nondiv3_orders hp hmod t g
  have hdiv := TraceData.exactly_two_div_3 t hmod
  -- Helper: When permIdx < 2, curveToOrder _ 2 ∈ {1, 5}
  have hcurve2_lt2 : ∀ idx : Fin 4, idx.val < 2 →
      curveToOrder idx 2 ∈ ({1, 5} : Finset (Fin 6)) := by
    intro ⟨i, hi⟩ hidx; interval_cases i <;> simp [curveToOrder]
  -- Helper: When permIdx ≥ 2, curveToOrder _ 2 ∈ {2, 4}
  have hcurve2_ge2 : ∀ idx : Fin 4, idx.val ≥ 2 →
      curveToOrder idx 2 ∈ ({2, 4} : Finset (Fin 6)) := by
    intro ⟨i, hi⟩ hidx; interval_cases i <;> simp [curveToOrder]
  -- Helper: When permIdx < 2, curveToOrder _ 4 ∈ {1, 5}
  have hcurve4_lt2 : ∀ idx : Fin 4, idx.val < 2 →
      curveToOrder idx 4 ∈ ({1, 5} : Finset (Fin 6)) := by
    intro ⟨i, hi⟩ hidx; interval_cases i <;> simp [curveToOrder]
  -- Helper: When permIdx ≥ 2, curveToOrder _ 4 ∈ {2, 4}
  have hcurve4_ge2 : ∀ idx : Fin 4, idx.val ≥ 2 →
      curveToOrder idx 4 ∈ ({2, 4} : Finset (Fin 6)) := by
    intro ⟨i, hi⟩ hidx; interval_cases i <;> simp [curveToOrder]
  -- If curveIdx ∈ {2, 4}, show the order would be composite
  by_contra hne
  push_neg at hne
  have hcurve_not_15 : curveIdx ≠ 1 ∧ curveIdx ≠ 5 := ⟨hne.1, hne.2⟩
  -- So curveIdx ∈ {2, 4}
  have hcurve_24 : curveIdx = 2 ∨ curveIdx = 4 := by
    fin_cases curveIdx <;> simp_all
  -- The curveOrder that curveIdx receives must be div by 3
  cases hdiv with
  | inl h15 =>
    -- curveOrders {1, 5} div by 3, {2, 4} not
    -- u0 = 0, permIdx < 2
    have hu0 : u0 t = 0 := by
      rw [u0_zero_iff t hmod]
      have h1_mod := TraceData.order1_mod3_eq t hmod
      have h1_div : curveOrder t 1 % 3 = 0 := Int.emod_eq_zero_of_dvd h15.1
      rw [h1_mod] at h1_div
      have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
        have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
        have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
        omega
      cases h_range with
      | inl h1 => exact h1
      | inr h2 =>
        have h_calc : (2 + t.a) % 3 = 1 := by
          have h_add : (2 + t.a) % 3 = (2 % 3 + t.a % 3) % 3 := Int.add_emod _ _ _
          rw [h_add, h2]; norm_num
        rw [h_calc] at h1_div; norm_num at h1_div
    have hperm_lt2 := u0_zero_permIdx_lt2 t g hu0
    -- curveIdx ∈ {2, 4} with permIdx < 2 gets curveOrders from {1, 5} which ARE div by 3
    cases hcurve_24 with
    | inl h2 =>
      rw [h2] at hPrime
      have hmem := hcurve2_lt2 (permIdx t g) hperm_lt2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      have hdiv3 : 3 ∣ curveOrder t (curveToOrder (permIdx t g) 2) := by
        cases hmem with
        | inl h => rw [h]; exact h15.1
        | inr h => rw [h]; exact h15.2.1
      have hgt := hgt_helper (curveToOrder (permIdx t g) 2)
      exact div_by_three_gt_three_not_prime hgt hdiv3 hPrime
    | inr h4 =>
      rw [h4] at hPrime
      have hmem := hcurve4_lt2 (permIdx t g) hperm_lt2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      have hdiv3 : 3 ∣ curveOrder t (curveToOrder (permIdx t g) 4) := by
        cases hmem with
        | inl h => rw [h]; exact h15.1
        | inr h => rw [h]; exact h15.2.1
      have hgt := hgt_helper (curveToOrder (permIdx t g) 4)
      exact div_by_three_gt_three_not_prime hgt hdiv3 hPrime
  | inr h24 =>
    -- curveOrders {2, 4} div by 3, {1, 5} not
    -- u0 = 1, permIdx ≥ 2
    have hu0 : u0 t = 1 := by
      rw [u0_one_iff t hmod]
      have h2_mod := TraceData.order2_mod3_eq t hmod
      have h2_div : curveOrder t 2 % 3 = 0 := Int.emod_eq_zero_of_dvd h24.1
      rw [h2_mod] at h2_div
      have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
        have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
        have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
        omega
      cases h_range with
      | inl h1 =>
        have h_calc : (2 - t.a) % 3 = 1 := by
          have h_sub : (2 - t.a) % 3 = (2 % 3 - t.a % 3) % 3 := Int.sub_emod _ _ _
          rw [h_sub, h1]; norm_num
        rw [h_calc] at h2_div; norm_num at h2_div
      | inr h2 => exact h2
    have hperm_ge2 := u0_one_permIdx_ge2 t g hu0
    -- curveIdx ∈ {2, 4} with permIdx ≥ 2 gets curveOrders from {2, 4} which ARE div by 3
    cases hcurve_24 with
    | inl h2 =>
      rw [h2] at hPrime
      have hmem := hcurve2_ge2 (permIdx t g) hperm_ge2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      have hdiv3 : 3 ∣ curveOrder t (curveToOrder (permIdx t g) 2) := by
        cases hmem with
        | inl h => rw [h]; exact h24.1
        | inr h => rw [h]; exact h24.2.1
      have hgt := hgt_helper (curveToOrder (permIdx t g) 2)
      exact div_by_three_gt_three_not_prime hgt hdiv3 hPrime
    | inr h4 =>
      rw [h4] at hPrime
      have hmem := hcurve4_ge2 (permIdx t g) hperm_ge2
      simp only [Finset.mem_insert, Finset.mem_singleton] at hmem
      have hdiv3 : 3 ∣ curveOrder t (curveToOrder (permIdx t g) 4) := by
        cases hmem with
        | inl h => rw [h]; exact h24.1
        | inr h => rw [h]; exact h24.2.1
      have hgt := hgt_helper (curveToOrder (permIdx t g) 4)
      exact div_by_three_gt_three_not_prime hgt hdiv3 hPrime

/-- Corollary: For any prime p ≡ 7 (mod 12), curves E_{g^1} and E_{g^5} are
    the only candidates for prime order among the six sextic twists. -/
theorem prime_order_only_at_curves_1_5 (p : ℕ) (hp : Nat.Prime p) (hmod : p % 12 = 7) :
    ∀ (t : TraceData p) (g : ZMod p) (curveIdx : Fin 6),
      Prime (curveOrder t (curveToOrder (permIdx t g) curveIdx)) →
      curveIdx ∈ ({1, 5} : Finset (Fin 6)) := by
  intro t g curveIdx hPrime
  have h := prime_order_implies_curve_1_or_5 hp hmod t g curveIdx hPrime
  cases h with
  | inl h1 => simp [h1]
  | inr h5 => simp [h5]
