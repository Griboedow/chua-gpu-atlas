import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

set_option autoImplicit false

namespace ChuaProof

/-- Algebraic energy balance for the stable mode plus oscillator realization.
This statement assumes the output relation; it does not construct the Chua flow. -/
theorem storage_balance (p r k q D z v phi x xd : ℝ)
    (h : x + q * xd = -p * z + k * v - D * phi) :
    p * z * (-r * z - phi) - k * v * phi + q * phi * xd =
      -x * phi - (D * phi^2 + 2 * p * z * phi + r * p * z^2) := by
  have hh := congrArg (fun t : ℝ => t * phi) h
  linear_combination hh

/-- The dissipation quadratic form is nonnegative under the stated strict margins. -/
theorem storage_dissipation (D p r z phi : ℝ)
    (hD : 0 < D) (hp : 0 < p) (hmargin : p < D * r) :
    0 ≤ D * phi^2 + 2 * p * z * phi + r * p * z^2 := by
  have hratio : p / D < r := (div_lt_iff₀ hD).2 (by nlinarith)
  have hnonneg : 0 ≤ D * (phi + p / D * z)^2 + p * (r - p / D) * z^2 :=
    add_nonneg (mul_nonneg (le_of_lt hD) (sq_nonneg _))
      (mul_nonneg (mul_nonneg (le_of_lt hp) (le_of_lt (sub_pos.mpr hratio))) (sq_nonneg _))
  have hid : D * (phi + p / D * z)^2 + p * (r - p / D) * z^2 =
      D * phi^2 + 2 * p * z * phi + r * p * z^2 := by
    field_simp
    ring
  rw [hid] at hnonneg
  exact hnonneg

/-- Sector positivity and the storage balance imply a nonpositive derivative.
Properness, the model coordinate change, and LaSalle's argument are separate. -/
theorem storage_nonincrease (D p r z phi x : ℝ)
    (hD : 0 < D) (hp : 0 < p) (hmargin : p < D * r) (hsector : 0 ≤ x * phi) :
    -x * phi - (D * phi^2 + 2 * p * z * phi + r * p * z^2) ≤ 0 := by
  have hq := storage_dissipation D p r z phi hD hp hmargin
  nlinarith

end ChuaProof

#print axioms ChuaProof.storage_balance
#print axioms ChuaProof.storage_dissipation
#print axioms ChuaProof.storage_nonincrease
