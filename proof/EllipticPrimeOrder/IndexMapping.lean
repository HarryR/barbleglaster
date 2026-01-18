/-
  IndexMapping.lean - Trace-to-curve-index mapping logic

  IMPORTANT: This module must match the Python implementation in
  lemma/2-eisenstein-mapping-magic.py exactly.
-/

import EllipticPrimeOrder.CubicChar
import EllipticPrimeOrder.Divisibility

/-! # Trace-to-Curve Index Mapping

This module formalizes the mapping from the curveOrder indices (0-5) to the
actual curve indices E_{g^j} for the six sextic twists of j-invariant 0 curves.

## Index Correspondence

The Lean curveOrder indexing differs from Python's norms_cd indexing:
  curveOrder[i] = norms_cd[(i + 3) % 6]

Specifically:
  curveOrder[0] = p + 1 + 2a = norms_cd[3]
  curveOrder[1] = p + 1 + a + 3b = norms_cd[4]
  curveOrder[2] = p + 1 - a + 3b = norms_cd[5]
  curveOrder[3] = p + 1 - 2a = norms_cd[0]
  curveOrder[4] = p + 1 - a - 3b = norms_cd[1]
  curveOrder[5] = p + 1 + a - 3b = norms_cd[2]

## Classification Bits (matching Python)

**u0 bit**: Determined by (c + d) mod 3, which equals a mod 3
- u0 = 0 when a ≡ 1 (mod 3): curveOrder indices {1, 5} ARE divisible by 3
- u0 = 1 when a ≡ 2 (mod 3): curveOrder indices {2, 4} ARE divisible by 3

**u1 bit**: Checks if ζ₃(g)·c + d ≡ 0 (mod p)

## The Four Permutations

The permutation index is: permIdx = 2 * u0 + u1 ∈ {0, 1, 2, 3}

Each permutation maps curveOrder indices to curve E_{g^j} indices such that
curves E_{g^1} and E_{g^5} always receive orders NOT divisible by 3.

Python's ALL_SUB_PATTERNS (in norms_cd indices):
  idx 0: (0, 1, 2, 3, 4, 5)
  idx 1: (0, 5, 4, 3, 2, 1)
  idx 2: (3, 4, 5, 0, 1, 2)
  idx 3: (3, 2, 1, 0, 5, 4)

Converted to Lean curveOrder indices (applying the shift):
  idx 0: curve j gets curveOrder[(j + 3) % 6]
  idx 1: curve j gets curveOrder[swapped_shifted[j]]
  idx 2: curve j gets curveOrder[j]  (identity)
  idx 3: curve j gets curveOrder[swapped[j]]
-/

open CubicChar

namespace IndexMapping

variable {p : ℕ}

/-! ## Classification Bit u0 -/

/-- u0 bit: based on (c + d) mod 3, equivalent to a mod 3.
    Returns 1 if a ≡ 2 (mod 3), returns 0 if a ≡ 1 (mod 3).
    MUST MATCH Python: u0 = int(((c+d) % 3) == 2) -/
def u0 (t : TraceData p) : Fin 2 :=
  if (t.c + t.d) % 3 = 2 then 1 else 0

/-- u0 characterization in terms of a mod 3.
    Since c = a + b and d = 2b, we have c + d = a + 3b ≡ a (mod 3). -/
theorem u0_eq_a_mod3 (t : TraceData p) :
    (t.c + t.d) % 3 = t.a % 3 := by
  -- c + d = (a + b) + 2b = a + 3b ≡ a (mod 3)
  simp only [TraceData.c, TraceData.d]
  have h : t.a + t.b + 2 * t.b = t.a + t.b * 3 := by ring
  rw [h]
  exact Int.add_mul_emod_self_right t.a t.b 3

/-- u0 = 0 iff a ≡ 1 (mod 3) -/
theorem u0_zero_iff (t : TraceData p) (hmod : p % 12 = 7) :
    u0 t = 0 ↔ t.a % 3 = 1 := by
  simp only [u0]
  have heq := u0_eq_a_mod3 t
  have hne := TraceData.a_not_div_3 t hmod
  have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
    have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
    have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
    omega
  constructor
  · intro h
    split_ifs at h with hcond
    · simp at h
    · rw [heq] at hcond
      cases h_range with
      | inl h1 => exact h1
      | inr h2 => exact absurd h2 hcond
  · intro ha1
    split_ifs with hcond
    · rw [heq, ha1] at hcond
      norm_num at hcond
    · rfl

/-- u0 = 1 iff a ≡ 2 (mod 3) -/
theorem u0_one_iff (t : TraceData p) (hmod : p % 12 = 7) :
    u0 t = 1 ↔ t.a % 3 = 2 := by
  simp only [u0]
  have heq := u0_eq_a_mod3 t
  have hne := TraceData.a_not_div_3 t hmod
  have h_range : t.a % 3 = 1 ∨ t.a % 3 = 2 := by
    have h_bound := Int.emod_lt_of_pos t.a (by norm_num : (3 : ℤ) > 0)
    have h_nonneg := Int.emod_nonneg t.a (by norm_num : (3 : ℤ) ≠ 0)
    omega
  constructor
  · intro h
    split_ifs at h with hcond
    · rw [heq] at hcond; exact hcond
    · simp at h
  · intro ha2
    split_ifs with hcond
    · rfl
    · rw [heq, ha2] at hcond; exact absurd rfl hcond

/-- Combined characterization: u0 = 0 ↔ a ≡ 1 (mod 3), u0 = 1 ↔ a ≡ 2 (mod 3) -/
theorem u0_characterization (t : TraceData p) (hmod : p % 12 = 7) :
    (u0 t = 0 ↔ t.a % 3 = 1) ∧ (u0 t = 1 ↔ t.a % 3 = 2) :=
  ⟨u0_zero_iff t hmod, u0_one_iff t hmod⟩

/-! ## Classification Bit u1 -/

/-- u1 bit: checks if ζ₃(g)·c + d ≡ 0 (mod p)
    MUST MATCH Python: u1 = int((ZETA(3,g,p) * c + d) % p == 0) -/
def u1 (t : TraceData p) (g : ZMod p) : Fin 2 :=
  if (cubicChar p g * t.c + t.d : ZMod p) = 0 then 1 else 0

/-! ## Permutation Index -/

/-- Combined permutation index: permIdx = 2 * u0 + u1 ∈ {0, 1, 2, 3}
    MUST MATCH Python: idx = (u0 * 2) + u1 -/
def permIdx (t : TraceData p) (g : ZMod p) : Fin 4 :=
  ⟨2 * (u0 t).val + (u1 t g).val, by
    have h0 : (u0 t).val < 2 := (u0 t).isLt
    have h1 : (u1 t g).val < 2 := (u1 t g).isLt
    omega⟩

/-! ## The Four Permutations

These permutations map curve indices E_{g^j} to curveOrder indices.
Given a curve E_{g^j}, the permutation tells us which curveOrder it receives.

The permutations are defined such that curves E_{g^1} and E_{g^5} always
receive curveOrders that are NOT divisible by 3 (and hence potentially prime).

Python ALL_SUB_PATTERNS converted to curveOrder indices:
| idx | u0 | u1 | a mod 3 | Div-by-3 curveOrders | Curves {1,5} get |
|-----|----|----|---------|----------------------|------------------|
| 0   | 0  | 0  | 1       | {1, 5}               | {4, 2} NOT div   |
| 1   | 0  | 1  | 1       | {1, 5}               | {2, 4} NOT div   |
| 2   | 1  | 0  | 2       | {2, 4}               | {1, 5} NOT div   |
| 3   | 1  | 1  | 2       | {2, 4}               | {5, 1} NOT div   |
-/

/-- The curve-to-curveOrder permutation for each permutation index.
    curveToOrder idx j = the curveOrder index that curve E_{g^j} receives.

    Derived from Python by converting norms_cd indices to curveOrder indices. -/
def curveToOrder : Fin 4 → Fin 6 → Fin 6
  | 0, j => ⟨(j.val + 3) % 6, Nat.mod_lt _ (by norm_num)⟩  -- shift by 3
  | 1, j => match j with | 0=>3 | 1=>2 | 2=>1 | 3=>0 | 4=>5 | 5=>4  -- shift + swap
  | 2, j => j                                               -- identity
  | 3, j => match j with | 0=>0 | 1=>5 | 2=>4 | 3=>3 | 4=>2 | 5=>1  -- swap pairs

/-! ## Key Properties -/

/-- Curves E_{g^1} and E_{g^5} get curveOrders {4, 2} when idx < 2 -/
theorem idx_lt2_curves_get_24 (idx : Fin 4) (h : idx.val < 2) :
    curveToOrder idx 1 ∈ ({2, 4} : Finset (Fin 6)) ∧
    curveToOrder idx 5 ∈ ({2, 4} : Finset (Fin 6)) := by
  fin_cases idx
  · simp [curveToOrder]
  · simp [curveToOrder]
  · simp at h
  · simp at h

/-- Curves E_{g^1} and E_{g^5} get curveOrders {1, 5} when idx ≥ 2 -/
theorem idx_ge2_curves_get_15 (idx : Fin 4) (h : idx.val ≥ 2) :
    curveToOrder idx 1 ∈ ({1, 5} : Finset (Fin 6)) ∧
    curveToOrder idx 5 ∈ ({1, 5} : Finset (Fin 6)) := by
  fin_cases idx
  · simp at h
  · simp at h
  · simp [curveToOrder]
  · simp [curveToOrder]

/-- When u0 = 0 (a ≡ 1 mod 3), permIdx < 2 -/
theorem u0_zero_permIdx_lt2 (t : TraceData p) (g : ZMod p)
    (hu0 : u0 t = 0) : (permIdx t g).val < 2 := by
  simp only [permIdx]
  have hu0_val : (u0 t).val = 0 := by rw [hu0]; rfl
  calc 2 * (u0 t).val + (u1 t g).val
      = 2 * 0 + (u1 t g).val := by rw [hu0_val]
    _ = (u1 t g).val := by ring
    _ < 2 := (u1 t g).isLt

/-- When u0 = 1 (a ≡ 2 mod 3), permIdx ≥ 2 -/
theorem u0_one_permIdx_ge2 (t : TraceData p) (g : ZMod p)
    (hu0 : u0 t = 1) : (permIdx t g).val ≥ 2 := by
  simp only [permIdx]
  have hu0_val : (u0 t).val = 1 := by rw [hu0]; rfl
  calc 2 * (u0 t).val + (u1 t g).val
      = 2 * 1 + (u1 t g).val := by rw [hu0_val]
    _ ≥ 2 := by omega

/-- Each permutation is a bijection (injective) -/
theorem perm_injective (idx : Fin 4) : Function.Injective (curveToOrder idx) := by
  intro a b hab
  fin_cases idx <;> fin_cases a <;> fin_cases b <;> simp_all [curveToOrder]

/-- Each permutation is a bijection (surjective) -/
theorem perm_surjective (idx : Fin 4) : Function.Surjective (curveToOrder idx) := by
  intro b
  fin_cases idx <;> fin_cases b <;> simp [curveToOrder] <;>
    (first | exact ⟨0, rfl⟩ | exact ⟨1, rfl⟩ | exact ⟨2, rfl⟩ |
            exact ⟨3, rfl⟩ | exact ⟨4, rfl⟩ | exact ⟨5, rfl⟩)

end IndexMapping
