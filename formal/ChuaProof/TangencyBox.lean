import Mathlib.Analysis.Calculus.LocalExtr.Rolle
import Mathlib.Analysis.Calculus.Deriv.Slope
import Mathlib.Topology.Order.Compact
import Mathlib.Tactic.Linarith

set_option autoImplicit false

/-! A finite-box existence argument, independent of any numerical implementation. -/

namespace ChuaProof

open Set Filter Topology

/-- At a minimum on a nontrivial closed interval, a negative left-endpoint
derivative rules out the left endpoint. -/
theorem not_min_at_left {f : ℝ → ℝ} {l r d : ℝ} (hlr : l < r)
    (hd : HasDerivAt f d l) (hdneg : d < 0) : ¬ IsMinOn f (Icc l r) l := by
  intro hmin
  have hs := (hasDerivAt_iff_tendsto_slope_left_right.mp hd).2
  have hn : 0 ≤ d := ge_of_tendsto hs (by
    filter_upwards [self_mem_nhdsWithin,
      (eventually_lt_nhds hlr).filter_mono nhdsWithin_le_nhds] with y hy hyr
    rw [slope_def_field]
    exact div_nonneg (sub_nonneg.mpr (hmin ⟨le_of_lt hy, le_of_lt hyr⟩))
      (sub_nonneg.mpr (le_of_lt hy)))
  exact (not_le_of_gt hdneg) hn

/-- The corresponding right-endpoint exclusion. -/
theorem not_min_at_right {f : ℝ → ℝ} {l r d : ℝ} (hlr : l < r)
    (hd : HasDerivAt f d r) (hdpos : 0 < d) : ¬ IsMinOn f (Icc l r) r := by
  intro hmin
  have hs := (hasDerivAt_iff_tendsto_slope_left_right.mp hd).1
  have hn : d ≤ 0 := le_of_tendsto hs (by
    filter_upwards [self_mem_nhdsWithin,
      (eventually_gt_nhds hlr).filter_mono nhdsWithin_le_nhds] with y hy hly
    rw [slope_def_field]
    exact div_nonpos_of_nonneg_of_nonpos
      (sub_nonneg.mpr (hmin ⟨le_of_lt hly, le_of_lt hy⟩))
      (sub_nonpos.mpr (le_of_lt hy)))
  exact (not_le_of_gt hdpos) hn

/-- Four uniform face signs force an interior zero which is stationary in `a`.

Continuity is required only on the closed rectangle. `ga` is connected to `g`
by an actual derivative hypothesis, so a zero of a separately named function
cannot stand in for stationarity. No convexity or parameter derivative is
needed for existence; these are additional hypotheses of `quadratic_contact`.
-/
theorem tangency_box
    (g ga : ℝ → ℝ → ℝ) {aL aR bL bR : ℝ}
    (ha : aL < aR) (hb : bL < bR)
    (hc : ContinuousOn (fun p : ℝ × ℝ => g p.1 p.2)
      (Icc aL aR ×ˢ Icc bL bR))
    (hd : ∀ a ∈ Icc aL aR, ∀ b ∈ Icc bL bR,
      HasDerivAt (fun x => g x b) (ga a b) a)
    (hbottom : ∀ a ∈ Icc aL aR, g a bL < 0)
    (htop : ∀ a ∈ Icc aL aR, 0 < g a bR)
    (hleft : ∀ b ∈ Icc bL bR, ga aL b < 0)
    (hright : ∀ b ∈ Icc bL bR, 0 < ga aR b) :
    ∃ a ∈ Ioo aL aR, ∃ b ∈ Ioo bL bR,
      g a b = 0 ∧ ga a b = 0 ∧ IsMinOn (fun x => g x b) (Icc aL aR) a := by
  let A := Icc aL aR
  let B := Icc bL bR
  letI : PreconnectedSpace B := isPreconnected_iff_preconnectedSpace.mp isPreconnected_Icc
  let f : B → A → ℝ := fun b a => g a b
  have hcont : Continuous (Function.uncurry f) :=
    hc.comp_continuous
      ((continuous_subtype_val.comp continuous_snd).prodMk
        (continuous_subtype_val.comp continuous_fst))
      (fun p => ⟨p.2.property, p.1.property⟩)
  let m : B → ℝ := fun b => sInf (f b '' (univ : Set A))
  have hm : Continuous m := isCompact_univ.continuous_sInf hcont
  have hne : (univ : Set A).Nonempty := ⟨⟨aL, le_rfl, ha.le⟩, mem_univ _⟩
  have hmin (b : B) : ∃ a : A, m b = g a b ∧
      IsMinOn (fun x => g x b) A a := by
    obtain ⟨a, _, heq, hle⟩ := isCompact_univ.exists_sInf_image_eq_and_le hne
      (hcont.comp (Continuous.prodMk_right b)).continuousOn
    exact ⟨a, heq, fun x hx => hle ⟨x, hx⟩ (mem_univ _)⟩
  let blo : B := ⟨bL, le_rfl, hb.le⟩
  let bhi : B := ⟨bR, hb.le, le_rfl⟩
  have hlo : m blo < 0 := by
    obtain ⟨a, heq, _⟩ := hmin blo
    rw [heq]
    exact hbottom a a.property
  have hhi : 0 < m bhi := by
    obtain ⟨a, heq, _⟩ := hmin bhi
    rw [heq]
    exact htop a a.property
  obtain ⟨b, hzero⟩ := intermediate_value_univ blo bhi hm ⟨hlo.le, hhi.le⟩
  have hbint : (b : ℝ) ∈ Ioo bL bR := by
    constructor
    · refine lt_of_le_of_ne b.property.1 ?_
      intro heq
      have heq' : b = blo := Subtype.ext heq.symm
      rw [heq'] at hzero
      linarith
    · refine lt_of_le_of_ne b.property.2 ?_
      intro heq
      have heq' : b = bhi := Subtype.ext heq
      rw [heq'] at hzero
      linarith
  obtain ⟨a, heq, hamin⟩ := hmin b
  have haint : (a : ℝ) ∈ Ioo aL aR := by
    constructor
    · refine lt_of_le_of_ne a.property.1 ?_
      intro heq
      have hh : IsMinOn (fun x => g x b) (Icc aL aR) aL := by
        simpa only [A, heq] using hamin
      exact not_min_at_left ha (hd aL ⟨le_rfl, ha.le⟩ b b.property)
        (hleft b b.property) hh
    · refine lt_of_le_of_ne a.property.2 ?_
      intro heq
      have hh : IsMinOn (fun x => g x b) (Icc aL aR) aR := by
        simpa only [A, heq] using hamin
      exact not_min_at_right ha (hd aR ⟨ha.le, le_rfl⟩ b b.property)
        (hright b b.property) hh
  refine ⟨a, haint, b, hbint, heq.symm.trans hzero, ?_, hamin⟩
  exact (hamin.isLocalMin (Icc_mem_nhds haint.1 haint.2)).hasDerivAt_eq_zero
    (hd a a.property b b.property)

/-- The same witness has nonzero second state derivative and parameter derivative
when the corresponding derivative bounds hold throughout the rectangle. -/
theorem quadratic_contact
    (g ga gaa gb : ℝ → ℝ → ℝ) {aL aR bL bR : ℝ}
    (ha : aL < aR) (hb : bL < bR)
    (hc : ContinuousOn (fun p : ℝ × ℝ => g p.1 p.2)
      (Icc aL aR ×ˢ Icc bL bR))
    (hd : ∀ a ∈ Icc aL aR, ∀ b ∈ Icc bL bR,
      HasDerivAt (fun x => g x b) (ga a b) a)
    (hbottom : ∀ a ∈ Icc aL aR, g a bL < 0)
    (htop : ∀ a ∈ Icc aL aR, 0 < g a bR)
    (hleft : ∀ b ∈ Icc bL bR, ga aL b < 0)
    (hright : ∀ b ∈ Icc bL bR, 0 < ga aR b)
    (hcurv : ∀ a ∈ Icc aL aR, ∀ b ∈ Icc bL bR,
      HasDerivAt (fun x => ga x b) (gaa a b) a ∧ 0 < gaa a b)
    (hparam : ∀ a ∈ Icc aL aR, ∀ b ∈ Icc bL bR,
      HasDerivAt (g a) (gb a b) b ∧ 0 < gb a b) :
    ∃ a ∈ Ioo aL aR, ∃ b ∈ Ioo bL bR,
      g a b = 0 ∧ HasDerivAt (fun x => g x b) 0 a ∧
      HasDerivAt (fun x => ga x b) (gaa a b) a ∧ 0 < gaa a b ∧
      HasDerivAt (g a) (gb a b) b ∧ 0 < gb a b := by
  obtain ⟨a, hai, b, hbi, hg, hga, _⟩ :=
    tangency_box g ga ha hb hc hd hbottom htop hleft hright
  have hac : a ∈ Icc aL aR := ⟨hai.1.le, hai.2.le⟩
  have hbc : b ∈ Icc bL bR := ⟨hbi.1.le, hbi.2.le⟩
  refine ⟨a, hai, b, hbi, hg, ?_, (hcurv a hac b hbc).1,
    (hcurv a hac b hbc).2, (hparam a hac b hbc).1, (hparam a hac b hbc).2⟩
  simpa only [hga] using hd a hac b hbc

end ChuaProof
