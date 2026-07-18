# Probe 2c5 — GATE of the v₁ trit (Wilson's blind derivation) — P1′–P4′ ALL PASS

**Date:** 2026-07-17  **Operator:** Mᵀ (flow), tower M_tower, q=3, β*=3/5, **no v₀ trit** (baseline g₀=0).
**Target:** v≥1 tower states (γ≡0 mod 3). **Verdict: P1′ P2′ P3′ P4′ = PASS at both L=2 and L=3.**

Level class (Wilson): `g:=γ/3`; **D9** = {γ≡0 mod9}; among v₃(γ)=1, bit **b=+1 ⟺ a≡γ/3 (mod3)**;
`η := +1 ⟺ a≡2g (mod3) = −b`.

## P1′ — level sets (e mod6) × {D9, U+, U−}
- r **well-defined** on `(e mod6) × class3` at both L. Values per e-group **≤3** (L=2: max 2 — D9 is empty for L=2 since it needs mod-9 depth; L=3: max 3).
- Distinct-value set = **{2/7, 11/35, 5/14, 2/5, 3/7}** — **identical to the baseline 5-set** on the finer (a mod9, γ mod27, e mod6) key. So the coarse trit key IS the exact resolution. **PASS.**

## P2′ — bit-blindness + opposite signs (U+ − U−, per e mod6)
| e mod6 | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| U+ − U− | **0 (BLIND)** | +1/14 | **−3/35** | −1/7 | **+3/35** | +1/14 |
- **e≡0 mod6 is b-BLIND** (two levels, D9 vs U). **e≡2 and e≡4 carry the bit with OPPOSITE signs** (−3/35, +3/35 = ∓1/7·… the ±1/7 E-sector fingerprint). Odd e all live. **PASS at both L.**

## P3′ — D9 DEEP mass = R₀(s₀), s₀ = −e mod 6
DEEP = flow to v′≥1. D9 deep move-mass (raw), L-invariant:
| e | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| deep(raw) | 65/189 | 34/189 | 20/189 | 16/189 | 20/189 | 34/189 |
| R₀(−e mod6) | R₀(0)=65/189 | R₀(5)=34/189 | R₀(4)=20/189 | R₀(3)=16/189 | R₀(2)=20/189 | R₀(1)=34/189 |
Exact match, all e, both L. **PASS.**

## P4′ — unit-g DEEP mass = W_D formula
`W_D = [R₀(s₂)+R₀(s₄)]/2 − η·[R_{D/2}(s₂)−R_{D/2}(s₄)]/2`, s₂=2−e, s₄=4−e.
Matches the **raw deep move-mass** at machine precision for every (e mod6, η), both L (e.g. e=2 η=+1 → 8/27; e=1 η=+1 → 16/189; e=1 η=−1 → 34/189). **PASS.**

## The one correction to the derivation (does NOT touch any claim)
The R-autocorrelation lives on the **mod-6 BASE SHELL (D=6)**, NOT the full tower depth D=2·3^{L−1}.
At L=2 they coincide (D=6); at L=3 the full-D=18 autocorrelation gives the WRONG numbers, the D=6 one is exact. Diagnostic: **the deep masses are L-INVARIANT** (identical /189 fractions at L=2 and L=3); the "7" in the {2/7, 3/7, …} ratios traces to Z=2⁶−1=63=9·7. So R₀/R_{D/2} are the **6-element base-shell** autocorrelations, universally.
Second (expected) note: P3′/P4′ are statements about the **raw move-weight** (the R-autocorrelation mass). The h-weighted flow coincides at even e (h-ratio=1) and splits from raw at odd e by the h-parity factor {5/3,4/3} — as it must; W_D is a move-weight object, not an h-weighted one.

## Status
v₁ trit **VINDICATED** (structure P1′/P2′ exact; closed-form masses P3′/P4′ exact on the mod-6 shell).
Rung-2 now has **both trits closed-form**: v₀ (τ, `r=4/9−κ(pop)W⁻(τ)`) and v₁ (b/η, W_D on the 6-shell).
**HELD BACK (per Wilson):** the joint (g₀,g₁) corrector search — waits for Wilson's blind derivation of the joint optimum before the reveal.

Probe: `probes/probe_phase2c5.py`; log `logs/probe_phase2c5_log.txt`.
