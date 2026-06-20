/-
Copyright (c) 2026 Nathan Humphrey. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Nathan Humphrey
-/
import Mathlib

set_option linter.style.longLine false
set_option linter.unusedVariables false

/-!
# F̂_p magnitude — Milestone 1a: the Step-5 orthogonality punchline

Part of the second-moment proof that `|G[a]|² = p^{r+1}·𝟙[a ≡ c mod p]`
(see `FHAT_EQUIDISTRIBUTION_SECOND_MOMENT_PROOF.md` in the Paper 4 folder).

This file proves the FINAL orthogonality step (Step 5). Once the autocorrelation
has been collapsed onto the order-`p` subgroup with a linear phase (Steps 2–4),
the outer transform over that subgroup is a single character-orthogonality sum on
`ZMod p`.

It mirrors the `h_orth` collapse already working in `Theorem3.lean` (which uses
`AddChar.sum_mulShift` + `ZMod.isPrimitive_stdAddChar`), here specialised to the
prime cyclic group `ZMod p`, with a subtraction in the argument so it lands
directly as `if a = c then p else 0`.

Status: this is the cleanest fully-provable piece. The novel inner-sum collapse
(Step 2, needs H1) and the LTE selection (Steps 3–4, need H2/H3) are stated in the
roadmap at the bottom and will be added in Milestone 1b.
-/

open scoped BigOperators

namespace Paper4Fhat

/-- **Step 5 (orthogonality punchline).**
On the prime cyclic group `ZMod p`, summing the standard additive character of
`(c − a)·j` over `j` gives `p` when `a = c`, and `0` otherwise.

This is the exact shape of the last line of the second-moment proof:
`∑_j e_p(j(c − a)) = p · 𝟙[a ≡ c]`. -/
theorem step5_orthogonality (p : ℕ) [NeZero p] (c a : ZMod p) :
    (∑ j : ZMod p, ZMod.stdAddChar ((c - a) * j)) = if a = c then (p : ℂ) else 0 := by
  -- put the sum into `mulShift` shape: argument `j * (c − a)`
  have hcomm : ∀ j : ZMod p,
      ZMod.stdAddChar ((c - a) * j) = ZMod.stdAddChar (j * (c - a)) := by
    intro j; rw [mul_comm]
  simp_rw [hcomm]
  -- character orthogonality — the same lemma used in Theorem3.lean's `h_orth`
  rw [AddChar.sum_mulShift (c - a) (ZMod.isPrimitive_stdAddChar p), ZMod.card]
  -- turn the `c − a = 0` test into `a = c`
  by_cases h : a = c
  · subst h; simp
  · have hne : c - a ≠ 0 := sub_ne_zero.mpr (fun hh => h hh.symm)
    simp [hne, h]

end Paper4Fhat

theorem inner_sum_collapse (p r : ℕ) [Fact p.Prime] (hodd : Odd p)
    [NeZero (p ^ r)] [NeZero (p ^ (r + 1))]
    (hpr1 : 1 < p ^ r) (y : ZMod (p ^ (r + 1))) :
    (∑ u' : ZMod (p ^ r),
        ZMod.stdAddChar (y * (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val))
      = ZMod.stdAddChar y * (if y * (p : ZMod (p ^ (r + 1))) = 0 then (p ^ r : ℂ) else 0) := by
  have hfac : ∀ s : ZMod (p ^ r),
      ZMod.stdAddChar
          (y * (1 + (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))))
        = ZMod.stdAddChar y
            * ZMod.stdAddChar
                (y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))) := by
    intro s
    rw [show y * (1 + (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1))))
          = y + y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1))) from by ring]
    rw [AddChar.map_add_eq_mul]
  have hsub :
      (∑ u' : ZMod (p ^ r),
          ZMod.stdAddChar (y * (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val))
        = ∑ s : ZMod (p ^ r),
            ZMod.stdAddChar
              (y * (1 + (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1))))) := by
    have hbij : ∃ e : ZMod (p ^ r) ≃ ZMod (p ^ r),
        ∀ u' : ZMod (p ^ r),
          (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
            = 1 + (p : ZMod (p ^ (r + 1))) * ((e u').val : ZMod (p ^ (r + 1))) := by
      classical
      -- (B) powers of (1+p) are distinct: orderOf (1+p) = p^r
      have hpow_inj : Function.Injective
          (fun u' : ZMod (p ^ r) => (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val) := by
        -- the order of (1+p) in ZMod (p^{r+1}) is exactly p^r  [LTE; next rung]
        -- 1+p is a unit (coprime to p^{r+1}); work in the units group
        have hunit : IsUnit (1 + (p : ZMod (p ^ (r + 1)))) := by
          have hcast : (1 + (p : ZMod (p ^ (r + 1)))) = ((p + 1 : ℕ) : ZMod (p ^ (r + 1))) := by
            push_cast; ring
          rw [hcast, ZMod.isUnit_iff_coprime]
          have hcop : (p + 1).Coprime p := by simp [Nat.coprime_self_add_left]
          exact hcop.pow_right _
        set u : (ZMod (p ^ (r + 1)))ˣ := hunit.unit with hu
        have huval : (u : ZMod (p ^ (r + 1))) = 1 + (p : ZMod (p ^ (r + 1))) := hunit.unit_spec
        -- the order of u in the units group is exactly p^r  [LTE; next rung]
        have hord : orderOf u = p ^ r := by
          have hpow1 : ∀ k : ℕ, u ^ k = 1 ↔ p ^ r ∣ k := by
            intro k
            rcases Nat.eq_zero_or_pos k with hk | hk
            · subst hk; simp
            have key : ((p + 1) ^ k - 1 : ℕ) % p ^ (r + 1) = 0 ↔ p ^ r ∣ k := by
              rw [← Nat.dvd_iff_mod_eq_zero]
              have hbase : 1 < p := (Fact.out : p.Prime).one_lt
              have hne : (p + 1) ^ k - 1 ≠ 0 := by
                have : 1 < (p + 1) ^ k := Nat.one_lt_pow hk.ne' (by omega)
                omega
              have hlte : padicValNat p ((p + 1) ^ k - 1) = 1 + padicValNat p k := by
                have hnd : ¬ p ∣ (p + 1) := by
                  rw [Nat.dvd_add_right (dvd_refl p), Nat.dvd_one]
                  exact (Fact.out : p.Prime).ne_one
                have h := padicValNat.pow_sub_pow (p := p) hodd
                    (x := p + 1) (y := 1) (by omega) (by simp) hnd hk.ne'
                simpa using h
              rw [padicValNat_dvd_iff_le hne, hlte, padicValNat_dvd_iff_le hk.ne']
              omega
            rw [← key]
            rw [show (((p + 1) ^ k - 1 : ℕ) % p ^ (r + 1) = 0
                  ↔ (p + 1) ^ k % p ^ (r + 1) = 1 % p ^ (r + 1)) from by
                rw [← Nat.dvd_iff_mod_eq_zero, ← Nat.modEq_iff_dvd' (Nat.one_le_pow _ _ (by omega))]
                exact ⟨fun h => h.symm, fun h => h.symm⟩]
            rw [← Nat.ModEq, Nat.ModEq.comm, ← ZMod.natCast_eq_natCast_iff, Nat.cast_pow,
                Nat.cast_one,
                show ((p + 1 : ℕ) : ZMod (p ^ (r + 1))) = 1 + (p : ZMod (p ^ (r + 1))) by
                  push_cast; ring,
                ← huval, ← Units.val_pow_eq_pow_val]
            exact ⟨fun h => by rw [h, Units.val_one], fun h => Units.val_eq_one.mp h.symm⟩
          have hpos : 0 < p ^ r := by positivity
          rw [orderOf_eq_iff hpos]
          refine ⟨(hpow1 (p ^ r)).2 dvd_rfl, ?_⟩
          intro m hm hm0 hcon
          exact absurd (Nat.le_of_dvd hm0 ((hpow1 m).1 hcon)) (Nat.not_le.2 hm)
        intro a b hab
        simp only [] at hab
        have habu : u ^ a.val = u ^ b.val := by
          apply Units.ext
          push_cast [huval]
          exact hab
        have hmod : a.val ≡ b.val [MOD p ^ r] := by
          have h := (pow_eq_pow_iff_modEq (x := u)).1 habu
          rwa [hord] at h
        have hval : a.val = b.val := by
          have ha := ZMod.val_lt a
          have hb := ZMod.val_lt b
          simp only [Nat.ModEq, Nat.mod_eq_of_lt ha, Nat.mod_eq_of_lt hb] at hmod
          exact hmod
        exact ZMod.val_injective (n := p ^ r) hval
      -- (C) every power is a principal unit:  (1+p)^k = 1 + p·ŝ  for some s
      have hpow_mem : ∀ u' : ZMod (p ^ r), ∃ s : ZMod (p ^ r),
          1 + (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))
            = (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val := by
        have hC1 : ∀ k : ℕ, ∃ w : ZMod (p ^ (r + 1)),
            (1 + (p : ZMod (p ^ (r + 1)))) ^ k = 1 + (p : ZMod (p ^ (r + 1))) * w := by
          intro k
          induction k with
          | zero => exact ⟨0, by ring⟩
          | succ k ih =>
            obtain ⟨w, hw⟩ := ih
            refine ⟨w + 1 + (p : ZMod (p ^ (r + 1))) * w, ?_⟩
            have hps : (1 + (p : ZMod (p ^ (r + 1)))) ^ (k + 1)
                = (1 + (p : ZMod (p ^ (r + 1)))) ^ k * (1 + (p : ZMod (p ^ (r + 1)))) :=
              pow_succ _ _
            linear_combination hps + (1 + (p : ZMod (p ^ (r + 1)))) * hw
        intro u'
        obtain ⟨w, hw⟩ := hC1 u'.val
        refine ⟨(w.val : ZMod (p ^ r)), ?_⟩
        conv_rhs => rw [hw]
        congr 1
        have hround : ((((w.val : ZMod (p ^ r)).val : ℕ)) : ZMod (p ^ (r + 1)))
            = (w.val : ZMod (p ^ (r + 1)))
              - (p : ZMod (p ^ (r + 1))) ^ r
                * ((w.val / p ^ r : ℕ) : ZMod (p ^ (r + 1))) := by
          have hc : ((p ^ r * (w.val / p ^ r) + (w.val : ZMod (p ^ r)).val : ℕ)
                      : ZMod (p ^ (r + 1)))
              = ((w.val : ℕ) : ZMod (p ^ (r + 1))) := by
            congr 1
            rw [ZMod.val_natCast]
            exact Nat.div_add_mod _ _
          push_cast at hc
          linear_combination hc
        rw [hround]
        rw [ZMod.natCast_val, ZMod.cast_id]
        have hpr' : (p : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r = 0 := by
          rw [← pow_succ']
          exact_mod_cast ZMod.natCast_self (p ^ (r + 1))
        linear_combination - ((w.val / p ^ r : ℕ) : ZMod (p ^ (r + 1))) * hpr'
      -- assemble e from (B)+(C): the choice function is an injective endo of a
      -- finite type, hence a bijection, and matches the powers pointwise by construction
      refine ⟨Equiv.ofBijective (fun u' => Classical.choose (hpow_mem u'))
        ((Finite.injective_iff_bijective).1 ?_), ?_⟩
      · intro a b hab
        apply hpow_inj
        simp only [] at hab ⊢
        rw [← Classical.choose_spec (hpow_mem a), ← Classical.choose_spec (hpow_mem b), hab]
      · intro u'
        rw [Equiv.ofBijective_apply]
        exact (Classical.choose_spec (hpow_mem u')).symm
    obtain ⟨e, he⟩ := hbij
    refine Fintype.sum_equiv e _ _ (fun u' => ?_)
    rw [he u']
  have horth :
      (∑ s : ZMod (p ^ r),
          ZMod.stdAddChar (y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))))
        = if y * (p : ZMod (p ^ (r + 1))) = 0 then (p ^ r : ℂ) else 0 := by
    -- the crux: (p : ZMod p^{r+1}) annihilates p^r, so s ↦ y·p·ŝ is additive
    have hpr : (p : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r = 0 := by
      rw [← pow_succ']
      exact_mod_cast ZMod.natCast_self (p ^ (r + 1))
    have hadd : ∀ s s' : ZMod (p ^ r),
        y * (p : ZMod (p ^ (r + 1))) * ((s + s').val : ZMod (p ^ (r + 1)))
          = y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))
            + y * (p : ZMod (p ^ (r + 1))) * (s'.val : ZMod (p ^ (r + 1))) := by
      intro s s'
      have hc : ((p ^ r * ((s.val + s'.val) / p ^ r) + (s + s').val : ℕ)
                  : ZMod (p ^ (r + 1)))
          = ((s.val + s'.val : ℕ) : ZMod (p ^ (r + 1))) := by
        congr 1
        rw [ZMod.val_add]
        exact Nat.div_add_mod _ _
      push_cast at hc
      linear_combination y * (p : ZMod (p ^ (r + 1))) * hc
        - y * (((s.val + s'.val) / p ^ r : ℕ) : ZMod (p ^ (r + 1))) * hpr
    let χ : AddChar (ZMod (p ^ r)) ℂ :=
      { toFun := fun s =>
          ZMod.stdAddChar (y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1))))
        map_zero_eq_one' := by simp
        map_add_eq_mul' := fun s s' => by
          rw [hadd s s', AddChar.map_add_eq_mul] }
    have hbridge :
        (∑ s : ZMod (p ^ r),
            ZMod.stdAddChar (y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))))
          = ∑ s : ZMod (p ^ r), χ s := rfl
    rw [hbridge]
    classical
    rw [AddChar.sum_eq_ite χ, ZMod.card]
    by_cases hyp : y * (p : ZMod (p ^ (r + 1))) = 0
    · -- y·p = 0: the character is trivial, so the sum is the full count p^r
      have hχ0 : χ = 0 := by
        apply AddChar.eq_zero_iff.2
        intro s
        change ZMod.stdAddChar (y * (p : ZMod (p ^ (r + 1))) * (s.val : ZMod (p ^ (r + 1)))) = 1
        rw [hyp, zero_mul, AddChar.map_zero_eq_one]
      rw [if_pos hχ0, if_pos hyp]
      push_cast
      ring
    · -- y·p ≠ 0: χ is nontrivial, so the sum vanishes
      have hχne : χ ≠ 0 := by
        rw [AddChar.ne_zero_iff]
        refine ⟨1, ?_⟩
        change ZMod.stdAddChar (y * (p : ZMod (p ^ (r + 1)))
            * ((1 : ZMod (p ^ r)).val : ZMod (p ^ (r + 1)))) ≠ 1
        haveI : Fact (1 < p ^ r) := ⟨hpr1⟩
        rw [ZMod.val_one, Nat.cast_one, mul_one]
        rw [Ne, (ZMod.isPrimitive_stdAddChar (p ^ (r + 1))).zmod_char_eq_one_iff]
        exact hyp
      rw [if_neg hχne, if_neg hyp]
  rw [hsub]
  simp_rw [hfac]
  rw [← Finset.mul_sum, horth]

theorem survivors (p r : ℕ) [Fact p.Prime] (hodd : Odd p) [NeZero (p ^ r)] (k : ℕ) :
    (1 + (p : ZMod (p ^ r))) ^ k = 1 ↔ p ^ (r - 1) ∣ k := by
  rcases Nat.eq_zero_or_pos k with hk | hk
  · subst hk; simp
  have key : ((p + 1) ^ k - 1 : ℕ) % p ^ r = 0 ↔ p ^ (r - 1) ∣ k := by
    rw [← Nat.dvd_iff_mod_eq_zero]
    have hbase : 1 < p := (Fact.out : p.Prime).one_lt
    have hne : (p + 1) ^ k - 1 ≠ 0 := by
      have : 1 < (p + 1) ^ k := Nat.one_lt_pow hk.ne' (by omega)
      omega
    have hlte : padicValNat p ((p + 1) ^ k - 1) = 1 + padicValNat p k := by
      have hnd : ¬ p ∣ (p + 1) := by
        rw [Nat.dvd_add_right (dvd_refl p), Nat.dvd_one]
        exact (Fact.out : p.Prime).ne_one
      have h := padicValNat.pow_sub_pow (p := p) hodd
          (x := p + 1) (y := 1) (by omega) (by simp) hnd hk.ne'
      simpa using h
    rw [padicValNat_dvd_iff_le hne, hlte, padicValNat_dvd_iff_le hk.ne']
    omega
  rw [← key]
  rw [show (((p + 1) ^ k - 1 : ℕ) % p ^ r = 0
        ↔ (p + 1) ^ k % p ^ r = 1 % p ^ r) from by
      rw [← Nat.dvd_iff_mod_eq_zero, ← Nat.modEq_iff_dvd' (Nat.one_le_pow _ _ (by omega))]
      exact ⟨fun h => h.symm, fun h => h.symm⟩]
  rw [← Nat.ModEq, Nat.ModEq.comm, ← ZMod.natCast_eq_natCast_iff, Nat.cast_pow,
      Nat.cast_one,
      show ((p + 1 : ℕ) : ZMod (p ^ r)) = 1 + (p : ZMod (p ^ r)) by push_cast; ring]
  exact eq_comm

theorem phase_value (p r : ℕ) [Fact p.Prime] (hodd : Odd p) [NeZero (p ^ (r + 1))]
    (hr : 1 ≤ r) (j : ℕ)
    (hbase : (1 + (p : ZMod (p ^ (r + 1)))) ^ (p ^ (r - 1))
      = 1 + (p : ZMod (p ^ (r + 1))) ^ r) :
    (1 + (p : ZMod (p ^ (r + 1)))) ^ (j * p ^ (r - 1))
      = 1 + (j : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r := by
  -- Step B: raise to the j-th power; p^r is nilpotent ((p^r)² = 0 since 2r ≥ r+1), so it linearizes
  have hnil : (p : ZMod (p ^ (r + 1))) ^ r * (p : ZMod (p ^ (r + 1))) ^ r = 0 := by
    have hz : (p : ZMod (p ^ (r + 1))) ^ (r + 1) = 0 := by
      rw [← Nat.cast_pow]; exact ZMod.natCast_self _
    calc (p : ZMod (p ^ (r + 1))) ^ r * (p : ZMod (p ^ (r + 1))) ^ r
        = (p : ZMod (p ^ (r + 1))) ^ (r + 1) * (p : ZMod (p ^ (r + 1))) ^ (r - 1) := by
          rw [← pow_add, ← pow_add]; congr 1; omega
      _ = 0 := by rw [hz, zero_mul]
  rw [mul_comm j (p ^ (r - 1)), pow_mul, hbase]
  have hpow : ∀ m : ℕ, (1 + (p : ZMod (p ^ (r + 1))) ^ r) ^ m
      = 1 + (m : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r := by
    intro m
    induction m with
    | zero => simp
    | succ m ih =>
      have hps : (1 + (p : ZMod (p ^ (r + 1))) ^ r) ^ (m + 1)
          = (1 + (p : ZMod (p ^ (r + 1))) ^ r) ^ m * (1 + (p : ZMod (p ^ (r + 1))) ^ r) :=
        pow_succ _ _
      push_cast
      linear_combination hps + (1 + (p : ZMod (p ^ (r + 1))) ^ r) * ih
        + (m : ZMod (p ^ (r + 1))) * hnil
  exact hpow j

/-- The length-`p^r` DFT of the principal-unit Gauss sum. -/
noncomputable def G (p r : ℕ) [NeZero (p ^ r)] [NeZero (p ^ (r + 1))]
    (c a : ZMod (p ^ (r + 1))) : ℂ :=
  ∑ u : ZMod (p ^ r),
    ZMod.stdAddChar (c * (1 + (p : ZMod (p ^ (r + 1)))) ^ u.val)
      * ZMod.stdAddChar (- ((p : ZMod (p ^ (r + 1))) * a * (u.val : ZMod (p ^ (r + 1)))))

theorem G_sq (p r : ℕ) [Fact p.Prime] (hodd : Odd p)
    [NeZero (p ^ r)] [NeZero (p ^ (r + 1))] (hpr1 : 1 < p ^ r) (hr : 1 ≤ r)
    (c a : ZMod (p ^ (r + 1))) (hc : IsUnit c)
    (hbase : (1 + (p : ZMod (p ^ (r + 1)))) ^ (p ^ (r - 1))
      = 1 + (p : ZMod (p ^ (r + 1))) ^ r) :
    G p r c a * (starRingEnd ℂ) (G p r c a)
      = if (ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) a
           = (ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) c
        then (p ^ (r + 1) : ℂ) else 0 := by
  classical
  -- Step 1: open the square into a double sum, change variables v = u − u'
  have step1 : G p r c a * (starRingEnd ℂ) (G p r c a)
      = ∑ v : ZMod (p ^ r),
          ZMod.stdAddChar (- ((p : ZMod (p ^ (r + 1))) * a * (v.val : ZMod (p ^ (r + 1)))))
            * (∑ u' : ZMod (p ^ r),
                ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
                  * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1)))) := by
    have hconj : ∀ x : ZMod (p ^ (r + 1)),
        (starRingEnd ℂ) (ZMod.stdAddChar x) = ZMod.stdAddChar (-x) := by
      intro x
      rw [ZMod.stdAddChar_apply, ← Circle.coe_inv_eq_conj, ← AddChar.map_neg_eq_inv,
          ← ZMod.stdAddChar_apply]
    have hopen : G p r c a * (starRingEnd ℂ) (G p r c a)
        = ∑ u : ZMod (p ^ r), ∑ u' : ZMod (p ^ r),
            ZMod.stdAddChar
              (c * (1 + (p : ZMod (p ^ (r + 1)))) ^ u.val
                  - (p : ZMod (p ^ (r + 1))) * a * (u.val : ZMod (p ^ (r + 1)))
                - (c * (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
                  - (p : ZMod (p ^ (r + 1))) * a * (u'.val : ZMod (p ^ (r + 1))))) := by
      unfold G
      rw [map_sum, Finset.sum_mul_sum]
      apply Finset.sum_congr rfl; intro u _
      apply Finset.sum_congr rfl; intro u' _
      simp only [map_mul, hconj]
      rw [← AddChar.map_add_eq_mul ZMod.stdAddChar, ← AddChar.map_add_eq_mul ZMod.stdAddChar,
          ← AddChar.map_add_eq_mul ZMod.stdAddChar]
      congr 1
      ring
    rw [hopen]
    have hg_ord : (1 + (p : ZMod (p ^ (r + 1)))) ^ (p ^ r) = 1 :=
      (survivors p (r + 1) hodd (p ^ r)).mpr (by simp)
    rw [show (∑ v : ZMod (p ^ r),
          ZMod.stdAddChar (- ((p : ZMod (p ^ (r + 1))) * a * (v.val : ZMod (p ^ (r + 1)))))
            * (∑ u' : ZMod (p ^ r),
                ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
                  * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1)))))
        = ∑ v : ZMod (p ^ r), ∑ u' : ZMod (p ^ r),
            ZMod.stdAddChar (- ((p : ZMod (p ^ (r + 1))) * a * (v.val : ZMod (p ^ (r + 1)))))
              * ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
                  * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))) from by
      apply Finset.sum_congr rfl; intro v _
      rw [Finset.mul_sum]]
    rw [Finset.sum_comm]
    conv_rhs => rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro y _
    refine (Fintype.sum_equiv (Equiv.addRight y) _ _ (fun v => ?_)).symm
    simp only [Equiv.coe_addRight]
    have hmod : ∀ n : ℕ, (1 + (p : ZMod (p ^ (r + 1)))) ^ (n % p ^ r)
        = (1 + (p : ZMod (p ^ (r + 1)))) ^ n := by
      intro n
      conv_rhs => rw [← Nat.mod_add_div n (p ^ r), pow_add, pow_mul, hg_ord, one_pow, mul_one]
    have hgfac : (1 + (p : ZMod (p ^ (r + 1)))) ^ (v + y).val
        = (1 + (p : ZMod (p ^ (r + 1)))) ^ v.val * (1 + (p : ZMod (p ^ (r + 1)))) ^ y.val := by
      rw [ZMod.val_add, hmod]
      exact pow_add _ _ _
    have hpfac : (p : ZMod (p ^ (r + 1))) * a * ((v + y).val : ZMod (p ^ (r + 1)))
        = (p : ZMod (p ^ (r + 1))) * a * (v.val : ZMod (p ^ (r + 1)))
          + (p : ZMod (p ^ (r + 1))) * a * (y.val : ZMod (p ^ (r + 1))) := by
      have hpr : (p : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r = 0 := by
        rw [← pow_succ']; exact_mod_cast ZMod.natCast_self (p ^ (r + 1))
      have hc2 : ((p ^ r * ((v.val + y.val) / p ^ r) + (v + y).val : ℕ) : ZMod (p ^ (r + 1)))
          = ((v.val + y.val : ℕ) : ZMod (p ^ (r + 1))) := by
        congr 1
        rw [ZMod.val_add]
        exact Nat.div_add_mod _ _
      push_cast at hc2
      linear_combination (p : ZMod (p ^ (r + 1))) * a * hc2
        - a * (((v.val + y.val) / p ^ r : ℕ) : ZMod (p ^ (r + 1))) * hpr
    rw [← AddChar.map_add_eq_mul ZMod.stdAddChar]
    congr 1
    rw [hgfac, hpfac]
    ring
  -- Step 2: inner sum collapses (inner_sum_collapse) to p^r · 𝟙[g^v ≡ 1 mod p^r]
  have step2 : ∀ v : ZMod (p ^ r),
      (∑ u' : ZMod (p ^ r),
          ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
            * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))))
        = ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))
            * (if p ^ (r - 1) ∣ v.val then (p ^ r : ℂ) else 0) := by
    intro v
    have hy := inner_sum_collapse p r hodd hpr1
      (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))
    rw [show (∑ u' : ZMod (p ^ r),
          ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val
            * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))))
        = ∑ u' : ZMod (p ^ r),
          ZMod.stdAddChar ((c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1))
            * (1 + (p : ZMod (p ^ (r + 1)))) ^ u'.val) from by
      apply Finset.sum_congr rfl; intro u' _; congr 1; ring]
    rw [hy]
    have hcond : (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1) * (p : ZMod (p ^ (r + 1))) = 0)
        ↔ p ^ (r - 1) ∣ v.val := by
      rw [mul_assoc, hc.mul_right_eq_zero]
      -- level reduction: y·p = 0 in ZMod(p^{r+1})  ↔  y reduces to 0 in ZMod(p^r)
      have hlevel : ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1) * (p : ZMod (p ^ (r + 1))) = 0
          ↔ (ZMod.castHom (pow_dvd_pow p (Nat.le_succ r)) (ZMod (p ^ r)))
              ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1) = 0 := by
        rw [ZMod.castHom_apply, ZMod.cast_eq_val, ZMod.natCast_eq_zero_iff]
        conv_lhs => rw [← ZMod.natCast_zmod_val ((1 + (p : ZMod (p ^ (r + 1)))) ^ v.val - 1)]
        rw [← Nat.cast_mul, ZMod.natCast_eq_zero_iff, pow_succ]
        exact Nat.mul_dvd_mul_iff_right (Fact.out : p.Prime).pos
      rw [hlevel, map_sub, map_pow, map_add, map_one, map_natCast, sub_eq_zero]
      exact survivors p r hodd v.val
    simp only [hcond]
  -- Step 3+4: only v with p^{r-1} ∣ v.val survive (survivors); on them the phase is e_p(c·j) (phase_value)
  have step34 : G p r c a * (starRingEnd ℂ) (G p r c a)
      = (p ^ r : ℂ) * ∑ j : ZMod p,
          ZMod.stdAddChar (- ((p : ZMod (p ^ (r + 1))) * a * ((j.val * p ^ (r - 1) : ℕ) : ZMod (p ^ (r + 1)))))
            * ZMod.stdAddChar (c * ((1 + (p : ZMod (p ^ (r + 1)))) ^ (j.val * p ^ (r - 1)) - 1)) := by
    rw [step1]
    simp_rw [step2]
    rw [Finset.mul_sum]
    refine Finset.sum_bij_ne_zero
      (fun x _ _ => ((x.val / p ^ (r - 1) : ℕ) : ZMod p)) ?_ ?_ ?_ ?_
    · intro x _ _; exact Finset.mem_univ _
    · intro a₁ _ h₁ a₂ _ h₂ heq
      simp only [] at heq
      have hd₁ : p ^ (r - 1) ∣ a₁.val := by
        by_contra hc'; apply h₁; simp [hc']
      have hd₂ : p ^ (r - 1) ∣ a₂.val := by
        by_contra hc'; apply h₂; simp [hc']
      have hlt₁ : a₁.val / p ^ (r - 1) < p := by
        rw [Nat.div_lt_iff_lt_mul (pow_pos (Fact.out : p.Prime).pos _)]
        have hpr : p * p ^ (r - 1) = p ^ r := by
          rw [mul_comm, ← pow_succ]; congr 1; omega
        rw [hpr]; exact ZMod.val_lt a₁
      have hlt₂ : a₂.val / p ^ (r - 1) < p := by
        rw [Nat.div_lt_iff_lt_mul (pow_pos (Fact.out : p.Prime).pos _)]
        have hpr : p * p ^ (r - 1) = p ^ r := by
          rw [mul_comm, ← pow_succ]; congr 1; omega
        rw [hpr]; exact ZMod.val_lt a₂
      have hq : a₁.val / p ^ (r - 1) = a₂.val / p ^ (r - 1) := by
        have h := (ZMod.natCast_eq_natCast_iff _ _ _).1 heq
        rwa [Nat.ModEq, Nat.mod_eq_of_lt hlt₁, Nat.mod_eq_of_lt hlt₂] at h
      apply ZMod.val_injective
      rw [← Nat.mul_div_cancel' hd₁, ← Nat.mul_div_cancel' hd₂, hq]
    · intro b _ hb
      have hblt : b.val < p := ZMod.val_lt b
      have hxlt : b.val * p ^ (r - 1) < p ^ r := by
        have hpr : p ^ r = p * p ^ (r - 1) := by rw [mul_comm, ← pow_succ]; congr 1; omega
        rw [hpr]
        exact mul_lt_mul_of_pos_right hblt (pow_pos (Fact.out : p.Prime).pos _)
      have hxval : ((b.val * p ^ (r - 1) : ℕ) : ZMod (p ^ r)).val = b.val * p ^ (r - 1) :=
        ZMod.val_natCast_of_lt hxlt
      refine ⟨((b.val * p ^ (r - 1) : ℕ) : ZMod (p ^ r)), Finset.mem_univ _, ?_, ?_⟩
      · have hdvd : p ^ (r - 1) ∣ ((b.val * p ^ (r - 1) : ℕ) : ZMod (p ^ r)).val := by
          rw [hxval]; exact dvd_mul_left _ _
        rw [if_pos hdvd, hxval]
        intro hz
        exact hb (by linear_combination hz)
      · simp only []
        rw [hxval, Nat.mul_div_cancel _ (pow_pos (Fact.out : p.Prime).pos _)]
        exact ZMod.natCast_zmod_val b
    · intro a₁ h₁ h₂
      have hd : p ^ (r - 1) ∣ a₁.val := by
        by_contra hc'; apply h₂; simp [hc']
      simp only []
      have hjlt : a₁.val / p ^ (r - 1) < p := by
        rw [Nat.div_lt_iff_lt_mul (pow_pos (Fact.out : p.Prime).pos _)]
        have hpr : p * p ^ (r - 1) = p ^ r := by rw [mul_comm, ← pow_succ]; congr 1; omega
        rw [hpr]; exact ZMod.val_lt a₁
      have hjval : ((a₁.val / p ^ (r - 1) : ℕ) : ZMod p).val = a₁.val / p ^ (r - 1) :=
        ZMod.val_natCast_of_lt hjlt
      rw [hjval, Nat.div_mul_cancel hd, if_pos hd]
      ring
  -- Step 5: final orthogonality on ZMod p (step5_orthogonality)
  rw [step34]
  have hp0 : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (Fact.out : p.Prime).pos.ne'
  have hcross : ∀ x : ZMod (p ^ (r + 1)),
      ZMod.stdAddChar (x * (p : ZMod (p ^ (r + 1))) ^ r)
        = ZMod.stdAddChar (ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p) x) := by
    intro x
    have hL : x * (p : ZMod (p ^ (r + 1))) ^ r
        = (((x.val * p ^ r : ℕ) : ℤ) : ZMod (p ^ (r + 1))) := by
      push_cast
      rw [ZMod.natCast_zmod_val]
    have hR : ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p) x
        = (((x.val : ℕ) : ℤ) : ZMod p) := by
      rw [ZMod.castHom_apply, ZMod.cast_eq_val]
      norm_cast
    rw [hL, hR, ZMod.stdAddChar_coe, ZMod.stdAddChar_coe]
    congr 1
    push_cast
    rw [pow_succ (p : ℂ) r]
    field_simp
  simp_rw [phase_value p r hodd hr _ hbase]
  have hclean : ∀ x : ZMod p,
      ZMod.stdAddChar (-((p : ZMod (p ^ (r + 1))) * a * ((x.val * p ^ (r - 1) : ℕ) : ZMod (p ^ (r + 1)))))
        * ZMod.stdAddChar (c * (1 + (x.val : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ r - 1))
      = ZMod.stdAddChar (((c - a) * (x.val : ZMod (p ^ (r + 1)))) * (p : ZMod (p ^ (r + 1))) ^ r) := by
    intro x
    rw [← AddChar.map_add_eq_mul ZMod.stdAddChar]
    congr 1
    have hpp : (p : ZMod (p ^ (r + 1))) * (p : ZMod (p ^ (r + 1))) ^ (r - 1)
        = (p : ZMod (p ^ (r + 1))) ^ r := by rw [← pow_succ']; congr 1; omega
    push_cast
    rw [← hpp]
    ring
  simp_rw [hclean]
  simp_rw [hcross]
  have hbridge : ∀ x : ZMod p,
      ZMod.stdAddChar ((ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) ((c - a) * (x.val : ZMod (p ^ (r + 1)))))
        = ZMod.stdAddChar (((ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) c
            - (ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) a) * x) := by
    intro x
    rw [map_mul, map_sub]
    congr 2
    rw [ZMod.castHom_apply, ZMod.natCast_val]
    exact ZMod.cast_cast_zmod_of_le (Nat.le_self_pow (Nat.succ_ne_zero r) p) x
  simp_rw [hbridge]
  rw [Paper4Fhat.step5_orthogonality p
        ((ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) c)
        ((ZMod.castHom (dvd_pow_self p (Nat.succ_ne_zero r)) (ZMod p)) a)]
  split_ifs with h
  · rw [pow_succ]
  · rw [mul_zero]
/-
