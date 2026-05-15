# IGUSA_G — Stationary phase / Watson on Igusa (asymptotic via residue)

## Phase 0 — verbatim (Denef Bourbaki §1.4)

**Corollary 1.4.5 (Igusa).** Suppose that C_f ∩ Supp Φ ⊂ f^{-1}(0). Then for |z| big enough, E_Φ(z) is a finite C-linear combination of functions z^λ · |z|^β · log^j |z| · ψ(ac z) with coefficients independent of z, where λ ∈ C is a pole of (s + 1) Z_Φ(s, χ_triv, f) or of Z_Φ(s, χ, f) for χ ≠ χ_triv, and β ∈ N with β ≤ (multiplicity of pole λ) - 1.

Equivalently (Watson-style): the asymptotic ∫_{K^n} Φ(x) |f(x)|^s |dx| as s → s_0 (a pole) gives the residue coefficient in the asymptotic expansion of E_Φ(z) = ∫_{K^n} Φ(x) Ψ(zf(x)) |dx| as |z| → ∞.

## Hypothesis types

- (i) Polynomial f, Schwartz-Bruhat Φ, C_f ∩ Supp Φ ⊂ f^{-1}(0).
- (ii) Identification: pole of Igusa zeta ↔ asymptotic term in oscillatory integral.

## Phase 1 — substrate check

R78 substrate g(u): C_g = {u : g(u) = 0, g'(u) = 0} ∩ Z_3.
- g(u) ≡ c mod 3, c ∈ (Z/3)*. So g(u) has **no roots in Z_3 ∩ {v_3 ≥ 1}** — g is a unit, takes values in Z_3*.
- C_g = ∅.
- The hypothesis C_f ∩ Supp Φ ⊂ f^{-1}(0) is **VACUOUSLY** satisfied (no critical points anywhere).

But then the asymptotic E_Φ(z) decays rapidly (since f has no zeros), and there are no poles of Z(s,g) to provide leading terms. **E_Φ(z) → 0 rapidly; no Watson-style asymptotic with leading term at z^{log_3(2)} arises.**

For substrate (2), cubic Postnikov phase P_a(s): pole at s=-1, residue gives the leading 1/z asymptotic of E_Φ(z) = ∫ Φ(s) Ψ(z P_a(s)) ds. This is just a finite oscillatory integral with single saddle.

## Phase 2 — conclusion shape

Watson-on-Igusa requires Igusa to FIRST produce a pole. Phase 1A/1B showed:
- Substrate 1: no poles.
- Substrate 2: pole at s = -1 only.

Watson does not produce new pole locations; it converts existing pole structure to asymptotic structure.

**Closure-relevant question:** if we had a pole at s_0 = log_3(2), Watson would give E_Φ(z) ~ z^{log_3(2)} = z^{0.631} as z → ∞ in the 3-adic norm. But z^{log_3(2)} is **growth, not decay** (positive exponent) — this would be inconsistent with E_Φ being bounded (which it is, by triangle inequality on the integral).

So **a pole at log_3(2) > 0 is structurally incompatible with the Watson-on-Igusa boundedness of E_Φ**. Confirms the categorical barrier identified in candidate B.

## Disposition: NO_FIT

Watson lemma adapted to Igusa is downstream of pole identification, which fails at Phase 1A/1B. Even granting log_3(2) as a hypothetical pole, the asymptotic shape (positive exponent in |z|^β with β > 0) is incompatible with the bounded oscillatory integral E_Φ(z).

This **reinforces** the categorical barrier: positive Re(s_0) ↔ growing asymptotic ↔ violates boundedness of the oscillatory integral. **NO_FIT.**
