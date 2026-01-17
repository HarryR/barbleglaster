/-
  EllipticPrimeOrder - Root module for prime order lemma formalization

  This library formalizes the lemma that for primes p ≡ 7 (mod 12) and
  j-invariant 0 elliptic curves E_{g^i}: y² = x³ + g^i, at most two of
  the six sextic twists can have prime order.

  ## Mathematical Background

  For a prime p ≡ 7 (mod 12):
  - There are exactly 6 sextic twists: E_{g^i} for i ∈ {0,1,2,3,4,5}
  - Curve order: |E_{g^i}| = p + 1 - t_i where t_i is the trace of Frobenius
  - Using Cornacchia: p = a² + 3b² with a even, b odd
  - Eisenstein coordinates: c = a + b, d = 2b

  ## Key Results

  1. **Evenness**: Orders at trace indices 0 and 3 are always even
     (t_0 = -2a and t_3 = 2a, both even since a is even)

  2. **Divisibility by 3**: Among orders at indices {1,2,4,5}, exactly two
     are divisible by 3 (which pair depends on a mod 3)

  3. **Main Theorem**: At most 2 of the 6 curve orders can be prime

  ## Module Structure

  - `Basic`: Imports and basic definitions
  - `Eisenstein`: Eisenstein integers Z[ω]
  - `TraceOfFrobenius`: TraceData structure and coefficient formulas
  - `Axioms`: Axiomatized results (Cornacchia, Hasse bound)
  - `Divisibility`: Even/divisibility-by-3 lemmas
  - `CubicChar`: Cubic character theory for the u1 classification bit
  - `IndexMapping`: Trace-to-curve index mapping permutations
  - `MainTheorem`: Main theorem assembly

  ## References

  - Bröker, R. (2005). Constructing elliptic curves of prescribed order
  - Washington, L. C. (2008). Elliptic Curves: Number Theory and Cryptography
  - Silverman, J. H. (1986, 1994). The Arithmetic of Elliptic Curves
-/

import EllipticPrimeOrder.Basic
import EllipticPrimeOrder.Eisenstein
import EllipticPrimeOrder.TraceOfFrobenius
import EllipticPrimeOrder.Axioms
import EllipticPrimeOrder.Divisibility
import EllipticPrimeOrder.CubicChar
import EllipticPrimeOrder.IndexMapping
import EllipticPrimeOrder.MainTheorem
