import Lake
open Lake DSL

package «EllipticPrimeOrder» where
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.26.0"

@[default_target]
lean_lib «EllipticPrimeOrder» where
  globs := #[.submodules `EllipticPrimeOrder]
