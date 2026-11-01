import Mathlib.Data.Real.Basic
import Mathlib.Tactic.NormNum

set_option autoImplicit false

namespace ChuaProof

/-- Exact rational checks of the ten graph-transform guard inequalities.
This does not assert that the Chua map satisfies the derivative bounds. -/
theorem graph_guards :
    let a : ℝ := 1.6469
    let b : ℝ := 0.00035
    let c : ℝ := 0.00024
    let d : ℝ := 0.00032
    let M : ℝ := 1500
    let Q : ℝ := 0.000041
    let L : ℝ := 0.001
    let K : ℝ := 2000
    let Γ : ℝ := 0.0001
    1 < a - b * L ∧ c * L + d < 1 ∧
    (c + d * L) / (a - b * L) < L ∧
    (b + L * d) / (a - L * c) < L ∧
    d + L * b < 1 ∧ 1 / (a - L * c) < 1 ∧
    (M * (1 + L)^3 + (d + L * b) * K) / (a - b * L)^2 < K ∧
    (M * (1 + L)^3 + K * (c * L + d)^2) / (a - L * c) < K ∧
    Q * (1 + L) + (d + L * b) * Γ < Γ ∧
    (Γ + Q * (1 + L)) / (a - L * c) < Γ := by
  norm_num

end ChuaProof
