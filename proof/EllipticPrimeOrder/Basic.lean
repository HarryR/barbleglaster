/-
  Basic.lean - Imports and basic definitions for prime order lemma
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Ring.Int.Parity
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Fin.Basic
import Mathlib.Tactic

/-! # Basic Definitions

This module provides foundational imports and basic definitions for the
formalization of the prime order lemma for j-invariant 0 elliptic curves.
-/

/-- A prime p is in our target class if p ≡ 7 (mod 12) -/
def IsPrimeMod7_12 (p : ℕ) : Prop := Nat.Prime p ∧ p % 12 = 7

/-- For primes p ≡ 7 (mod 12), we have p > 3 -/
lemma prime_mod7_12_gt_three {p : ℕ} (h : IsPrimeMod7_12 p) : p > 3 := by
  obtain ⟨hp, hmod⟩ := h
  omega

/-- For primes p ≡ 7 (mod 12), p is odd -/
lemma prime_mod7_12_odd {p : ℕ} (h : IsPrimeMod7_12 p) : Odd p := by
  obtain ⟨hp, hmod⟩ := h
  have : p % 2 = 1 := by omega
  exact Nat.odd_iff.mpr this

/-- For primes p ≡ 7 (mod 12), p ≡ 1 (mod 3) -/
lemma prime_mod7_12_mod3 {p : ℕ} (h : IsPrimeMod7_12 p) : p % 3 = 1 := by
  obtain ⟨_, hmod⟩ := h
  omega
