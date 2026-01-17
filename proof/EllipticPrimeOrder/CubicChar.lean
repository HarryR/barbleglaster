/-
  CubicChar.lean - Cubic character theory for the trace-to-curve-index mapping
-/

import EllipticPrimeOrder.TraceOfFrobenius

/-! # Cubic Character Theory

For p ≡ 1 (mod 3), the multiplicative group F_p* contains elements of order 3.
The cubic character χ₃(g) = g^((p-1)/3) maps elements to cube roots of unity.

## Key Mathematical Insights

From the Eisenstein norm: c² - cd + d² = p

We can show: c³ + d³ = (c + d)(c² - cd + d²) = (c + d) · p ≡ 0 (mod p)

This means: (-d/c)³ ≡ 1 (mod p)

**Therefore: -d/c is always a cube root of unity in F_p!**

The u1 classification bit checks whether ζ₃(g) = -d/c, which determines
the sub-permutation within pairs for the trace-to-curve mapping.
-/

namespace CubicChar

variable {p : ℕ}

/-- Cubic character: g^((p-1)/3) mod p -/
def cubicChar (p : ℕ) (g : ZMod p) : ZMod p := g ^ ((p - 1) / 3)

/-- Helper: For p ≡ 1 (mod 3), (p-1) is divisible by 3 -/
theorem p_sub_one_div_3 (hp : Nat.Prime p) (hmod : p % 3 = 1) : 3 ∣ (p - 1) := by
  have hp_pos : 0 < p := hp.pos
  have h : (p - 1) % 3 = 0 := by omega
  exact Nat.dvd_of_mod_eq_zero h

/-- For p ≡ 1 (mod 3), the cubic character cubed equals 1.
    Proof: g^((p-1)/3)^3 = g^(p-1) = 1 by Fermat's little theorem -/
theorem cubicChar_cubed (hp : Nat.Prime p) (hmod : p % 3 = 1) (g : ZMod p) (hg : g ≠ 0) :
    (cubicChar p g) ^ 3 = 1 := by
  haveI : Fact (Nat.Prime p) := ⟨hp⟩
  simp only [cubicChar]
  have hdiv := p_sub_one_div_3 hp hmod
  have hp1 : (p - 1) / 3 * 3 = p - 1 := Nat.div_mul_cancel hdiv
  rw [← pow_mul, hp1]
  exact ZMod.pow_card_sub_one_eq_one hg

/-- The cubic character is a cube root of unity (one of 1, ω, or ω²) -/
theorem cubicChar_is_cubeRoot (hp : Nat.Prime p) (hmod : p % 3 = 1) (g : ZMod p) (hg : g ≠ 0) :
    let χ := cubicChar p g
    χ ^ 3 = 1 := cubicChar_cubed hp hmod g hg

/-- Key theorem: From c² - cd + d² = p, we get -d/c is a cube root of unity in F_p.
    Proof: c³ + d³ = (c + d)(c² - cd + d²) = (c + d) · p ≡ 0 (mod p)
    Thus (-d/c)³ = -d³/c³ = c³/c³ = 1 (mod p) -/
theorem neg_d_div_c_cubed [Fact (Nat.Prime p)] (t : TraceData p)
    (hc : (t.c : ZMod p) ≠ 0) :
    (-(t.d : ZMod p) / (t.c : ZMod p)) ^ 3 = 1 := by
  -- c² - cd + d² = p
  have hnorm := t.eisenstein_norm_eq_p
  -- c³ + d³ = (c + d)(c² - cd + d²) = (c + d) · p
  have hfactor : t.c^3 + t.d^3 = (t.c + t.d) * (t.c^2 - t.c * t.d + t.d^2) := by ring
  have hcp : t.c^3 + t.d^3 = (t.c + t.d) * (p : ℤ) := by rw [hfactor, hnorm]
  -- In ZMod p, (c + d) * p = 0
  have hmod_p : ((t.c^3 + t.d^3 : ℤ) : ZMod p) = 0 := by
    rw [hcp]
    simp [mul_comm]
  -- So (c : ZMod p)³ + (d : ZMod p)³ = 0
  have hcube_sum : (t.c : ZMod p)^3 + (t.d : ZMod p)^3 = 0 := by
    have h1 : ((t.c^3 : ℤ) : ZMod p) = (t.c : ZMod p)^3 := by push_cast; ring
    have h2 : ((t.d^3 : ℤ) : ZMod p) = (t.d : ZMod p)^3 := by push_cast; ring
    calc (t.c : ZMod p)^3 + (t.d : ZMod p)^3
        = ((t.c^3 : ℤ) : ZMod p) + ((t.d^3 : ℤ) : ZMod p) := by rw [h1, h2]
      _ = ((t.c^3 + t.d^3 : ℤ) : ZMod p) := by push_cast; ring
      _ = 0 := hmod_p
  -- Therefore (-d/c)³ = -d³/c³ = c³/c³ = 1
  have hc_ne : (t.c : ZMod p) ≠ 0 := hc
  have hc3_ne : (t.c : ZMod p)^3 ≠ 0 := by
    intro h
    rw [pow_eq_zero_iff (by norm_num : 3 ≠ 0)] at h
    exact hc_ne h
  -- From c³ + d³ = 0, we get d³ = -c³
  have hd3 : (t.d : ZMod p)^3 = -(t.c : ZMod p)^3 := by
    have h := hcube_sum
    have : (t.d : ZMod p)^3 = 0 - (t.c : ZMod p)^3 := by
      rw [← h]; ring
    simp at this
    exact this
  -- (-d/c)³ = (-d)³/c³ = -d³/c³ = -(-c³)/c³ = c³/c³ = 1
  calc (-(t.d : ZMod p) / (t.c : ZMod p))^3
      = (-(t.d : ZMod p))^3 / (t.c : ZMod p)^3 := by rw [div_pow]
    _ = -(t.d : ZMod p)^3 / (t.c : ZMod p)^3 := by ring_nf
    _ = -(-(t.c : ZMod p)^3) / (t.c : ZMod p)^3 := by rw [hd3]
    _ = (t.c : ZMod p)^3 / (t.c : ZMod p)^3 := by ring_nf
    _ = 1 := div_self hc3_ne

/-- The u1 classification bit: checks if ζ₃(g)·c + d ≡ 0 (mod p)
    This is equivalent to checking if ζ₃(g) = -d/c -/
def u1_condition (p : ℕ) (t : TraceData p) (g : ZMod p) : Prop :=
  cubicChar p g * (t.c : ZMod p) + (t.d : ZMod p) = 0

/-- u1 condition is decidable -/
instance (p : ℕ) (t : TraceData p) (g : ZMod p) : Decidable (u1_condition p t g) :=
  inferInstanceAs (Decidable (_ = _))

/-- When c ≠ 0, u1_condition is equivalent to ζ₃(g) = -d/c -/
theorem u1_equiv_neg_d_div_c [Fact (Nat.Prime p)] (t : TraceData p) (g : ZMod p)
    (hc : (t.c : ZMod p) ≠ 0) :
    u1_condition p t g ↔ cubicChar p g = -(t.d : ZMod p) / (t.c : ZMod p) := by
  simp only [u1_condition]
  constructor
  · intro h
    -- ζ₃(g)·c + d = 0  →  ζ₃(g)·c = -d  →  ζ₃(g) = -d/c
    have heq : cubicChar p g * (t.c : ZMod p) = -(t.d : ZMod p) := by
      have : cubicChar p g * (t.c : ZMod p) + (t.d : ZMod p) = 0 := h
      calc cubicChar p g * (t.c : ZMod p)
          = cubicChar p g * (t.c : ZMod p) + (t.d : ZMod p) - (t.d : ZMod p) := by ring
        _ = 0 - (t.d : ZMod p) := by rw [this]
        _ = -(t.d : ZMod p) := by ring
    field_simp
    exact heq
  · intro h
    -- ζ₃(g) = -d/c  →  ζ₃(g)·c = -d  →  ζ₃(g)·c + d = 0
    have heq : cubicChar p g * (t.c : ZMod p) = -(t.d : ZMod p) := by
      rw [h]
      field_simp
    calc cubicChar p g * (t.c : ZMod p) + (t.d : ZMod p)
        = -(t.d : ZMod p) + (t.d : ZMod p) := by rw [heq]
      _ = 0 := by ring

end CubicChar
