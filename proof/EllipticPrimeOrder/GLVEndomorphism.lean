/-
  GLVEndomorphism.lean - GLV β parameter computation from Eisenstein structure

  This module proves that the GLV endomorphism parameter β must be one of
  the two values derivable from Eisenstein coordinates: {-d/c, (-d/c)²}.

  The specific choice between these two is determined empirically (see README).
-/

import EllipticPrimeOrder.CubicChar

/-! # GLV Endomorphism Parameter

For j-invariant 0 elliptic curves, the GLV endomorphism φ(x,y) = (βx, y)
requires β to be a primitive cube root of unity in F_p.

From the Eisenstein structure, we know -d/c is a cube root of unity.
Therefore β ∈ {-d/c, (-d/c)²} (the two primitive cube roots).

## What This Module Proves

1. -d/c is a cube root of unity (from CubicChar.lean)
2. The primitive cube roots are exactly {-d/c, (-d/c)²} (when -d/c ≠ 1)
3. Any valid GLV β must be in this set

## What Is Verified Empirically

The specific formula β = g^(2(p-1)/3) = ζ₃(g)² is verified computationally.
See README.md for details.
-/

namespace GLVEndomorphism

open CubicChar

variable {p : ℕ} [Fact (Nat.Prime p)]

/-! ## Cube Root Structure -/

/-- The three cube roots of unity are 1, ω, ω² where ω³ = 1 -/
theorem cube_roots_of_unity (x : ZMod p) (hx : x ^ 3 = 1) :
    x = 1 ∨ x ^ 2 + x + 1 = 0 := by
  -- x³ - 1 = (x - 1)(x² + x + 1)
  have hfactor : x ^ 3 - 1 = (x - 1) * (x ^ 2 + x + 1) := by ring
  have h0 : x ^ 3 - 1 = 0 := by rw [hx]; ring
  rw [hfactor] at h0
  rcases mul_eq_zero.mp h0 with hsub | hquad
  · left; exact sub_eq_zero.mp hsub
  · right; exact hquad

/-- If x is a primitive cube root (x³ = 1, x ≠ 1), then x² is also primitive -/
theorem sq_primitive_cube_root (x : ZMod p) (hx3 : x ^ 3 = 1) (hx1 : x ≠ 1) :
    x ^ 2 ≠ 1 ∧ (x ^ 2) ^ 3 = 1 := by
  constructor
  · -- x² ≠ 1
    intro hx2
    -- If x² = 1 and x³ = 1, then x = x³/x² = 1
    have : x = x ^ 3 / x ^ 2 := by field_simp
    rw [hx3, hx2] at this
    simp at this
    exact hx1 this
  · -- (x²)³ = 1
    calc (x ^ 2) ^ 3 = x ^ 6 := by ring
      _ = (x ^ 3) ^ 2 := by ring
      _ = 1 ^ 2 := by rw [hx3]
      _ = 1 := one_pow 2

/-- The two primitive cube roots are x and x² -/
theorem primitive_cube_roots_are_x_and_x_sq (x y : ZMod p)
    (hx3 : x ^ 3 = 1) (hx1 : x ≠ 1)
    (hy3 : y ^ 3 = 1) (hy1 : y ≠ 1) :
    y = x ∨ y = x ^ 2 := by
  -- Both x and y satisfy the quadratic t² + t + 1 = 0
  have hx_quad : x ^ 2 + x + 1 = 0 := by
    cases cube_roots_of_unity x hx3 with
    | inl h => exact absurd h hx1
    | inr h => exact h
  have hy_quad : y ^ 2 + y + 1 = 0 := by
    cases cube_roots_of_unity y hy3 with
    | inl h => exact absurd h hy1
    | inr h => exact h
  -- x ≠ x² (since x ≠ 1)
  have hx_ne_x2 : x ≠ x ^ 2 := by
    intro heq
    have : x * (x - 1) = 0 := by
      calc x * (x - 1) = x ^ 2 - x := by ring
        _ = x - x := by rw [← heq]
        _ = 0 := by ring
    rcases mul_eq_zero.mp this with hx0 | hxsub
    · rw [hx0] at hx3; simp at hx3
    · exact hx1 (sub_eq_zero.mp hxsub)
  -- Key: (y - x)(y - x²) = y² - y(x + x²) + x³ = y² - y(-1) + 1 = y² + y + 1 = 0
  -- Since y² + y + 1 = 0 and x + x² = -1 (from x² + x + 1 = 0)
  have hsum : x + x ^ 2 = -1 := by
    have : x ^ 2 + x + 1 = 0 := hx_quad
    calc x + x ^ 2 = x ^ 2 + x := by ring
      _ = x ^ 2 + x + 1 - 1 := by ring
      _ = 0 - 1 := by rw [this]
      _ = -1 := by ring
  have hprod : (y - x) * (y - x ^ 2) = 0 := by
    calc (y - x) * (y - x ^ 2)
        = y ^ 2 - y * x - y * x ^ 2 + x * x ^ 2 := by ring
      _ = y ^ 2 - y * (x + x ^ 2) + x ^ 3 := by ring
      _ = y ^ 2 - y * (-1) + 1 := by rw [hsum, hx3]
      _ = y ^ 2 + y + 1 := by ring
      _ = 0 := hy_quad
  -- So y = x or y = x²
  rcases mul_eq_zero.mp hprod with hy_x | hy_x2
  · left; exact sub_eq_zero.mp hy_x
  · right; exact sub_eq_zero.mp hy_x2

/-! ## GLV β Must Be Primitive -/

/-- GLV β is a primitive cube root of unity (not 1) -/
structure IsPrimitiveCubeRoot (β : ZMod p) : Prop where
  cubed_eq_one : β ^ 3 = 1
  ne_one : β ≠ 1

/-- From Eisenstein coordinates, -d/c is a cube root of unity -/
theorem neg_d_over_c_is_cube_root (t : TraceData p) (hc : (t.c : ZMod p) ≠ 0) :
    (-(t.d : ZMod p) / (t.c : ZMod p)) ^ 3 = 1 :=
  neg_d_div_c_cubed t hc

/-! ## Main Result: β ∈ {-d/c, (-d/c)²} -/

/-- Any primitive cube root of unity equals -d/c or (-d/c)² -/
theorem glv_beta_in_eisenstein_set (t : TraceData p) (β : ZMod p)
    (hc : (t.c : ZMod p) ≠ 0)
    (hβ : IsPrimitiveCubeRoot β)
    (hω : -(t.d : ZMod p) / (t.c : ZMod p) ≠ 1) :
    β = -(t.d : ZMod p) / (t.c : ZMod p) ∨
    β = (-(t.d : ZMod p) / (t.c : ZMod p)) ^ 2 := by
  let ω := -(t.d : ZMod p) / (t.c : ZMod p)
  have hω3 : ω ^ 3 = 1 := neg_d_over_c_is_cube_root t hc
  exact primitive_cube_roots_are_x_and_x_sq ω β hω3 hω hβ.cubed_eq_one hβ.ne_one

/-- Equivalent formulation: β ∈ {-d/c, (-d/c)²} -/
theorem glv_beta_computable (t : TraceData p) (β : ZMod p)
    (hc : (t.c : ZMod p) ≠ 0)
    (hβ_cubed : β ^ 3 = 1) (hβ_ne_1 : β ≠ 1)
    (hω_ne_1 : -(t.d : ZMod p) / (t.c : ZMod p) ≠ 1) :
    β = -(t.d : ZMod p) / (t.c : ZMod p) ∨
    β = (-(t.d : ZMod p) / (t.c : ZMod p)) ^ 2 :=
  glv_beta_in_eisenstein_set t β hc ⟨hβ_cubed, hβ_ne_1⟩ hω_ne_1

/-! ## Computational Formula

Empirically verified: β = g^(2(p-1)/3) for primitive root g.

This corresponds to β = (-d/c)² when ζ₃(g) = -d/c (u1 condition holds),
or β = -d/c when ζ₃(g) ≠ -d/c (u1 condition fails).

See README.md section "GLV β Direct Computation" for verification.
-/

/-- The two candidate values for GLV β -/
def glvBetaCandidates (t : TraceData p) : ZMod p × ZMod p :=
  let ω := -(t.d : ZMod p) / (t.c : ZMod p)
  (ω, ω ^ 2)

end GLVEndomorphism
