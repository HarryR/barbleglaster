/-
  Divisibility.lean - Proofs about evenness and divisibility by 3
-/

import EllipticPrimeOrder.Axioms

/-! # Divisibility Lemmas

This module proves the key divisibility properties that constrain which curves
can have prime order:

1. **Evenness**: Orders at indices 0 and 3 are always even.
   - order_0 = p + 1 + 2a (even since p+1 even and 2a even)
   - order_3 = p + 1 - 2a (even)

2. **Divisibility by 3**: Among the remaining orders {1, 2, 4, 5}, exactly two
   are divisible by 3. Which pair depends on a (mod 3):
   - If a ≡ 1 (mod 3): orders at indices 1 and 5 are divisible by 3
   - If a ≡ 2 (mod 3): orders at indices 2 and 4 are divisible by 3
-/

namespace TraceData

variable {p : ℕ} (t : TraceData p)

/-! ## Evenness of Traces -/

/-- Trace at index 0 is even (it equals -2a where a is even) -/
theorem trace_0_even : Even (traceOfFrobenius t 0) := by
  rw [trace_0_eq]
  obtain ⟨k, hk⟩ := t.hEvenA
  exact ⟨-2 * k, by linarith⟩

/-- Trace at index 3 is even (it equals 2a where a is even) -/
theorem trace_3_even : Even (traceOfFrobenius t 3) := by
  rw [trace_3_eq]
  obtain ⟨k, hk⟩ := t.hEvenA
  exact ⟨2 * k, by linarith⟩

/-! ## Evenness of Orders -/

/-- For p ≡ 7 (mod 12), p + 1 is even -/
theorem p_plus_one_even (hmod : p % 12 = 7) : Even ((p : ℤ) + 1) := by
  have hp_odd : p % 2 = 1 := by omega
  exact ⟨(p + 1) / 2, by omega⟩

/-- Order at index 0 is even.
    Proof: order_0 = p + 1 + 2a. Since p + 1 = 2m and a = 2k,
    order_0 = 2m + 4k = 2(m + 2k). -/
theorem order_0_even (hmod : p % 12 = 7) : Even (curveOrder t 0) := by
  rw [order_0_eq]
  obtain ⟨k, hk⟩ := t.hEvenA
  have ⟨m, hm⟩ := p_plus_one_even hmod
  -- a = 2k, p + 1 = 2m, so p + 1 + 2a = 2m + 4k = 2(m + 2k)
  use m + 2 * k
  calc (p : ℤ) + 1 + 2 * t.a = 2 * m + 2 * t.a := by rw [hm]; ring
    _ = 2 * m + 2 * (2 * k) := by rw [hk]; ring
    _ = (m + 2 * k) + (m + 2 * k) := by ring

/-- Order at index 3 is even.
    Proof: order_3 = p + 1 - 2a. Since p + 1 = 2m and a = 2k,
    order_3 = 2m - 4k = 2(m - 2k). -/
theorem order_3_even (hmod : p % 12 = 7) : Even (curveOrder t 3) := by
  rw [order_3_eq]
  obtain ⟨k, hk⟩ := t.hEvenA
  have ⟨m, hm⟩ := p_plus_one_even hmod
  use m - 2 * k
  calc (p : ℤ) + 1 - 2 * t.a = 2 * m - 2 * t.a := by rw [hm]; ring
    _ = 2 * m - 2 * (2 * k) := by rw [hk]; ring
    _ = (m - 2 * k) + (m - 2 * k) := by ring

/-! ## Modular Arithmetic Properties -/

/-- For p ≡ 7 (mod 12), we have p ≡ 1 (mod 3) -/
theorem p_mod_3 (hmod : p % 12 = 7) : p % 3 = 1 := by omega

/-- From p = a² + 3b², we have a ≢ 0 (mod 3).
    Proof: If 3 | a, then 3 | a² + 3b² = p, but p ≡ 1 (mod 3). -/
theorem a_not_div_3 (hmod : p % 12 = 7) : t.a % 3 ≠ 0 := by
  intro h
  have hp3 : p % 3 = 1 := p_mod_3 hmod
  have hcorn := t.hCornacchia
  have ha0 : (3 : ℤ) ∣ t.a := Int.dvd_of_emod_eq_zero h
  obtain ⟨k, hk⟩ := ha0
  have h3_dvd_sum : (3 : ℤ) ∣ t.a^2 + 3*t.b^2 := ⟨3*k^2 + t.b^2, by rw [hk]; ring⟩
  have hp_int : (p : ℤ) = t.a^2 + 3*t.b^2 := hcorn.symm
  have h3_dvd_p : (3 : ℤ) ∣ (p : ℤ) := by rw [hp_int]; exact h3_dvd_sum
  obtain ⟨q, hq⟩ := h3_dvd_p
  have h_p_mod : (p : ℤ) % 3 = 0 := by simp [hq]
  have hp3_int : (p : ℤ) % 3 = 1 := by omega
  omega

/-- a² ≡ 1 (mod 3) since a ≢ 0 (mod 3).
    This follows from Fermat's little theorem: a² ≡ 1 (mod 3) when gcd(a, 3) = 1. -/
theorem a_sq_mod_3 (hmod : p % 12 = 7) : t.a^2 % 3 = 1 := by
  have hne : t.a % 3 ≠ 0 := a_not_div_3 t hmod
  -- For integers, a % 3 ∈ {0, 1, 2} when 3 > 0
  -- Since a % 3 ≠ 0, we have a % 3 ∈ {1, 2}
  have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
    have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
    have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
    omega
  -- Key: a^2 % 3 = (a * a) % 3 = ((a % 3) * (a % 3)) % 3 = (a % 3)^2 % 3
  have hsq_mod : t.a^2 % 3 = (t.a % 3) * (t.a % 3) % 3 := by
    have : t.a^2 = t.a * t.a := sq t.a
    rw [this, Int.mul_emod]
  rcases h_range with h1 | h2
  · -- Case a ≡ 1 (mod 3): a² ≡ 1 (mod 3)
    rw [hsq_mod, h1]
    norm_num
  · -- Case a ≡ 2 (mod 3): a² ≡ 4 ≡ 1 (mod 3)
    rw [hsq_mod, h2]
    norm_num

/-! ## Divisibility by 3 for Orders -/

/-- order_1 ≡ order_5 (mod 3).
    Both orders differ by ±6b, which is divisible by 3. -/
theorem order_1_5_cong_mod3 : curveOrder t 1 % 3 = curveOrder t 5 % 3 := by
  rw [order_1_eq, order_5_eq]
  -- order_1 = (p + 1 + a) + 3b, order_5 = (p + 1 + a) - 3b
  -- Both are congruent to (p + 1 + a) mod 3
  have h1 : ((p : ℤ) + 1 + t.a + 3*t.b) % 3 = ((p : ℤ) + 1 + t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 + t.a + 3*t.b = ((p : ℤ) + 1 + t.a) + t.b * 3 by ring]
    exact Int.add_mul_emod_self_right ((p : ℤ) + 1 + t.a) t.b 3
  have h2 : ((p : ℤ) + 1 + t.a - 3*t.b) % 3 = ((p : ℤ) + 1 + t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 + t.a - 3*t.b = ((p : ℤ) + 1 + t.a) - t.b * 3 by ring]
    exact Int.sub_mul_emod_self_right ((p : ℤ) + 1 + t.a) t.b 3
  rw [h1, h2]

/-- order_2 ≡ order_4 (mod 3).
    Both orders differ by ±6b, which is divisible by 3. -/
theorem order_2_4_cong_mod3 : curveOrder t 2 % 3 = curveOrder t 4 % 3 := by
  rw [order_2_eq, order_4_eq]
  -- order_2 = (p + 1 - a) + 3b, order_4 = (p + 1 - a) - 3b
  -- Both are congruent to (p + 1 - a) mod 3
  have h1 : ((p : ℤ) + 1 - t.a + 3*t.b) % 3 = ((p : ℤ) + 1 - t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 - t.a + 3*t.b = ((p : ℤ) + 1 - t.a) + t.b * 3 by ring]
    exact Int.add_mul_emod_self_right ((p : ℤ) + 1 - t.a) t.b 3
  have h2 : ((p : ℤ) + 1 - t.a - 3*t.b) % 3 = ((p : ℤ) + 1 - t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 - t.a - 3*t.b = ((p : ℤ) + 1 - t.a) - t.b * 3 by ring]
    exact Int.sub_mul_emod_self_right ((p : ℤ) + 1 - t.a) t.b 3
  rw [h1, h2]

/-- Product of order_1 and order_2 is divisible by 3.
    Proof: (p+1+a+3b)(p+1-a+3b) = (p+1+3b)² - a² ≡ 1 - 1 = 0 (mod 3). -/
theorem order_product_div_3 (hmod : p % 12 = 7) :
    3 ∣ curveOrder t 1 * curveOrder t 2 := by
  rw [order_1_eq, order_2_eq]
  have hp3 : (p : ℤ) % 3 = 1 := by exact_mod_cast p_mod_3 hmod
  have ha2 : t.a^2 % 3 = 1 := a_sq_mod_3 t hmod
  -- Rewrite product as difference of squares: (m + a)(m - a) = m² - a² where m = p + 1 + 3b
  have heq : ((p:ℤ) + 1 + t.a + 3*t.b) * ((p:ℤ) + 1 - t.a + 3*t.b) =
             ((p:ℤ) + 1 + 3*t.b)^2 - t.a^2 := by ring
  rw [heq]
  -- m % 3 = (p + 1) % 3 = 2 since 3b ≡ 0 (mod 3)
  have hm_mod : ((p : ℤ) + 1 + 3*t.b) % 3 = 2 := by
    have h1 : ((p : ℤ) + 1 + 3*t.b) % 3 = ((p : ℤ) + 1) % 3 := by
      conv_lhs => rw [show (p : ℤ) + 1 + 3*t.b = ((p : ℤ) + 1) + t.b * 3 by ring]
      exact Int.add_mul_emod_self_right ((p : ℤ) + 1) t.b 3
    rw [h1]
    omega
  -- m² % 3 = (m % 3)² % 3 = 4 % 3 = 1
  have hm2_mod : ((p : ℤ) + 1 + 3*t.b)^2 % 3 = 1 := by
    have hsq : ((p : ℤ) + 1 + 3*t.b)^2 % 3 = (((p : ℤ) + 1 + 3*t.b) % 3) * (((p : ℤ) + 1 + 3*t.b) % 3) % 3 := by
      have : ((p : ℤ) + 1 + 3*t.b)^2 = ((p : ℤ) + 1 + 3*t.b) * ((p : ℤ) + 1 + 3*t.b) := sq _
      rw [this, Int.mul_emod]
    rw [hsq, hm_mod]
    norm_num
  -- (m² - a²) % 3 = 1 - 1 = 0, so 3 | (m² - a²)
  have hdiff : (((p : ℤ) + 1 + 3*t.b)^2 - t.a^2) % 3 = 0 := by
    have h := Int.sub_emod (((p : ℤ) + 1 + 3*t.b)^2) (t.a^2) 3
    rw [hm2_mod, ha2] at h
    omega
  exact Int.dvd_of_emod_eq_zero hdiff

/-- At least one of {order_1, order_5} or {order_2, order_4} is divisible by 3.
    Proof: Since 3 | order_1 * order_2 and 3 is prime, either 3 | order_1 or 3 | order_2.
    The congruences then give us the pairs. -/
theorem some_orders_div_3 (hmod : p % 12 = 7) :
    (3 ∣ curveOrder t 1 ∧ 3 ∣ curveOrder t 5) ∨
    (3 ∣ curveOrder t 2 ∧ 3 ∣ curveOrder t 4) := by
  have hdiv := order_product_div_3 t hmod
  have h15 := order_1_5_cong_mod3 t
  have h24 := order_2_4_cong_mod3 t
  have h3 : Prime (3 : ℤ) := Int.prime_three
  have h_or := h3.dvd_or_dvd hdiv
  cases h_or with
  | inl h1 =>
    left
    constructor
    · exact h1
    · -- 3 | order_5 follows from h1 and h15
      have h1_mod : curveOrder t 1 % 3 = 0 := Int.emod_eq_zero_of_dvd h1
      rw [h15] at h1_mod
      exact Int.dvd_of_emod_eq_zero h1_mod
  | inr h2 =>
    right
    constructor
    · exact h2
    · -- 3 | order_4 follows from h2 and h24
      have h2_mod : curveOrder t 2 % 3 = 0 := Int.emod_eq_zero_of_dvd h2
      rw [h24] at h2_mod
      exact Int.dvd_of_emod_eq_zero h2_mod

/-! ## Helper lemmas for modular arithmetic with 2 ± a -/

private lemma mod3_2_plus (a : ℤ) (ha : a % 3 = 1) : (2 + a) % 3 = 0 := by
  have h : (2 + a) % 3 = (2 % 3 + a % 3) % 3 := Int.add_emod _ _ _
  rw [h, ha]; norm_num

private lemma mod3_2_plus' (a : ℤ) (ha : a % 3 = 2) : (2 + a) % 3 = 1 := by
  have h : (2 + a) % 3 = (2 % 3 + a % 3) % 3 := Int.add_emod _ _ _
  rw [h, ha]; norm_num

private lemma mod3_2_minus (a : ℤ) (ha : a % 3 = 1) : (2 - a) % 3 = 1 := by
  have h : (2 - a) % 3 = (2 % 3 - a % 3) % 3 := Int.sub_emod _ _ _
  rw [h, ha]; norm_num

private lemma mod3_2_minus' (a : ℤ) (ha : a % 3 = 2) : (2 - a) % 3 = 0 := by
  have h : (2 - a) % 3 = (2 % 3 - a % 3) % 3 := Int.sub_emod _ _ _
  rw [h, ha]; norm_num

/-- order_1 % 3 = (2 + a) % 3 when p ≡ 1 (mod 3) -/
theorem order1_mod3_eq (hmod : p % 12 = 7) :
    curveOrder t 1 % 3 = (2 + t.a) % 3 := by
  rw [order_1_eq]
  have hp3 : (p : ℤ) % 3 = 1 := by exact_mod_cast p_mod_3 hmod
  have h1 : ((p : ℤ) + 1 + t.a + 3*t.b) % 3 = ((p : ℤ) + 1 + t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 + t.a + 3*t.b = ((p : ℤ) + 1 + t.a) + t.b * 3 by ring]
    exact Int.add_mul_emod_self_right ((p : ℤ) + 1 + t.a) t.b 3
  rw [h1]
  have hp1 : ((p : ℤ) + 1) % 3 = 2 := by omega
  have h2 : ((p : ℤ) + 1 + t.a) % 3 = (((p : ℤ) + 1) % 3 + t.a % 3) % 3 := Int.add_emod _ _ _
  have h3 : (2 + t.a) % 3 = (2 % 3 + t.a % 3) % 3 := Int.add_emod _ _ _
  rw [h2, hp1, h3]; norm_num

/-- order_2 % 3 = (2 - a) % 3 when p ≡ 1 (mod 3) -/
theorem order2_mod3_eq (hmod : p % 12 = 7) :
    curveOrder t 2 % 3 = (2 - t.a) % 3 := by
  rw [order_2_eq]
  have hp3 : (p : ℤ) % 3 = 1 := by exact_mod_cast p_mod_3 hmod
  have h1 : ((p : ℤ) + 1 - t.a + 3*t.b) % 3 = ((p : ℤ) + 1 - t.a) % 3 := by
    conv_lhs => rw [show (p : ℤ) + 1 - t.a + 3*t.b = ((p : ℤ) + 1 - t.a) + t.b * 3 by ring]
    exact Int.add_mul_emod_self_right ((p : ℤ) + 1 - t.a) t.b 3
  rw [h1]
  have hp1 : ((p : ℤ) + 1) % 3 = 2 := by omega
  have h2 : ((p : ℤ) + 1 - t.a) % 3 = (((p : ℤ) + 1) % 3 - t.a % 3) % 3 := Int.sub_emod _ _ _
  have h3 : (2 - t.a) % 3 = (2 % 3 - t.a % 3) % 3 := Int.sub_emod _ _ _
  rw [h2, hp1, h3]; norm_num

/-- Exactly two of {order_1, order_2, order_4, order_5} are divisible by 3.
    Which pair depends on a (mod 3):
    - If a ≡ 1 (mod 3): orders 1 and 5 are divisible by 3
    - If a ≡ 2 (mod 3): orders 2 and 4 are divisible by 3 -/
theorem exactly_two_div_3 (hmod : p % 12 = 7) :
    (3 ∣ curveOrder t 1 ∧ 3 ∣ curveOrder t 5 ∧ ¬(3 ∣ curveOrder t 2) ∧ ¬(3 ∣ curveOrder t 4)) ∨
    (3 ∣ curveOrder t 2 ∧ 3 ∣ curveOrder t 4 ∧ ¬(3 ∣ curveOrder t 1) ∧ ¬(3 ∣ curveOrder t 5)) := by
  have h1_mod := order1_mod3_eq t hmod
  have h2_mod := order2_mod3_eq t hmod
  have h15 := order_1_5_cong_mod3 t
  have h24 := order_2_4_cong_mod3 t
  have ha_ne := a_not_div_3 t hmod
  have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
    have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
    have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
    omega
  cases h_range with
  | inl ha1 =>
    -- a % 3 = 1: order_1, order_5 divisible by 3
    left
    have h1_div : curveOrder t 1 % 3 = 0 := by rw [h1_mod, mod3_2_plus t.a ha1]
    have h2_not : curveOrder t 2 % 3 ≠ 0 := by rw [h2_mod, mod3_2_minus t.a ha1]; norm_num
    refine ⟨?_, ?_, ?_, ?_⟩
    · exact Int.dvd_of_emod_eq_zero h1_div
    · have h5_div : curveOrder t 5 % 3 = 0 := by rw [← h15, h1_div]
      exact Int.dvd_of_emod_eq_zero h5_div
    · intro h
      have : curveOrder t 2 % 3 = 0 := Int.emod_eq_zero_of_dvd h
      exact h2_not this
    · intro h
      have h4_div : curveOrder t 4 % 3 = 0 := Int.emod_eq_zero_of_dvd h
      rw [h24, h4_div] at h2_not
      exact h2_not rfl
  | inr ha2 =>
    -- a % 3 = 2: order_2, order_4 divisible by 3
    right
    have h1_not : curveOrder t 1 % 3 ≠ 0 := by rw [h1_mod, mod3_2_plus' t.a ha2]; norm_num
    have h2_div : curveOrder t 2 % 3 = 0 := by rw [h2_mod, mod3_2_minus' t.a ha2]
    refine ⟨?_, ?_, ?_, ?_⟩
    · exact Int.dvd_of_emod_eq_zero h2_div
    · have h4_div : curveOrder t 4 % 3 = 0 := by rw [← h24, h2_div]
      exact Int.dvd_of_emod_eq_zero h4_div
    · intro h
      have : curveOrder t 1 % 3 = 0 := Int.emod_eq_zero_of_dvd h
      exact h1_not this
    · intro h
      have h5_div : curveOrder t 5 % 3 = 0 := Int.emod_eq_zero_of_dvd h
      rw [h15, h5_div] at h1_not
      exact h1_not rfl

end TraceData
