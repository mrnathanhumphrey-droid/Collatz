# Probe F1 — the period-9 frontier — **rate 0.984 (not ½) confirmed to k=12; period is NON-integer (all integers 2–9 excluded, exact signs); the oscillation is p-adic-LATTICE-forced (archimedean non-lattice) — Wilson's premise correct**

**Date:** 2026-07-22. Probe `probes/probe_F1_period9.py`. Data: exact `ε_k` through k=8, float `ε_k` k=9–12
(`result_epsilon_11.csv` + `S_12`), so `Λ_r = (ε_{r+1}−ε_r)/2` through r=11. Pins the structure of the critical-only
period-~9 mode that R81 established as the true asymptotic rate of `S_n → 7/15`.

## F1-A — the exact ladder (single reference artifact)
| k | ε_k | \|ε\|·2^k | \|ε\|/0.984^k | sign | ‖ | r | Λ_r | \|Λ\|·2^r | Λ_{r+1}/Λ_r |
|---|---|---|---|---|---|---|---|---|---|
| 4 | −2.45e−3 | 0.039 | 0.0026 | − | | 4 | +6.50e−4 | 0.010 | **0.4927** |
| 5 | −1.15e−3 | 0.037 | 0.0012 | − | | 5 | +3.27e−4 | 0.010 | **0.5028** |
| 6 | −4.98e−4 | 0.032 | 0.0006 | − | | 6 | −3.39e−4 | 0.022 | −1.036 |
| 8 | −7.46e−4 | **0.191** | 0.0008 | − | | 8 | +3.69e−4 | 0.094 | 0.987 |
| 10 | +7.21e−4 | **0.738** | 0.0008 | + | | 10 | +3.91e−4 | 0.400 | 0.989 |
| 12 | +2.27e−3 | **9.317** | 0.0028 | + | | 11 | +3.86e−4 | 0.791 | — |

**Rate 0.984, not ½, confirmed to k=12:** `|ε_k|·2^k` explodes (0.038 → 9.32), so `ε` decays far slower than `2^{−k}`;
`|ε_k|/0.984^k` stays **bounded and oscillating** (~0.001–0.003). And `Λ_r` is essentially **flat** at r=7–11
(~3.7×10⁻⁴, ratios hovering at 1: 0.987, 1.073, 0.989) while `|Λ_r|·2^r` **doubles** each step — the ½ region
(r=4,5) is over by r=6. The corpus's `ρ≈0.984` is the real rate.

## F1-B — INTEGER-vs-IRRATIONAL period: **every integer period 2–9 is EXCLUDED** (exact signs)
`Λ` sign sequence r=1..11: **`− − + + + − + + + + +`** (exact-robust). Testing `sign(Λ_{r+p}) = sign(Λ_r)`:

| p | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 12, 18 |
|---|---|---|---|---|---|---|---|---|---|
| | EXCL | EXCL | EXCL | EXCL | EXCL | EXCL | EXCL | **EXCL** | untestable (need ≥13,19 terms) |

**No integer period 2–9 is consistent** with the 11-term exact sign sequence. An **irrational** period produces no
exact integer sign-period — so this is exactly consistent with the lattice-theoretic `2π/log2 = 9.0647` (irrational)
and with any non-integer near 9. **The period is not an integer** — the decisive structural read, immune to the
<1-period data limitation.

## F1-C — where the mode lives: the **intermediate-j lift tower** (not the rank-1 η=1 part)
`M_n(1+3^j)/S_n` (float π, n=2..7): the highest lift `j=n−1` is `−½` exact (the rank-1 η=1 part, W2-A); `j=0 → ½`,
`j=1 → ¼` **stabilize**; but the **intermediate** entries **oscillate** — e.g. `j=2`: 0.236, 0.056, 0.063, −0.011
(n=4..7); `j=3`: −0.105, −0.016, 0.038. The mode lives in these oscillating middle-tower entries, **not** in the
η=1 triple (rank-1) and **not** in the stabilizing low-j entries. Successive-difference ratios are noisy (0.10,
−0.08, 0.48, 2.05 for j=0), confirming it is **not a clean finite eigenvalue** — the growing-tower/no-finite-operator
character (R29) plus the oscillation. The mode is located, not diagonalized.

## F1-D — the 0.970/0.984 relation (numerology cell, unearned)
`√0.97008 = 0.984926` vs the second-moment envelope ~0.984 (`|Δ| = 9.3×10⁻⁴`). So `√(first-moment shifted mode) ≈
second-moment rate` to ~10⁻³ — **cannot be confirmed or denied at current precision.** Flagged as the numerology
class that killed five leads: a probe cell, not a finding, until a mechanism relates the two.

## F1-E — LATTICE check: archimedean **non-lattice** (Wilson correct); oscillation is **p-adic-lattice-forced**
- **Archimedean step `log M = log3 − v·log2` is NON-LATTICE.** Feller-lattice needs a `d>0` with all values in `dℤ`
  ⟺ `log3, log2` both multiples of `d` ⟺ `log3/log2 ∈ ℚ`. But `log3/log2 = 1.5849625…` is **irrational**
  (elementary: `3^q = 2^p` is impossible by unique factorization). **No common `d` ⟹ non-lattice** — Wilson's premise
  is **correct**. (The "coset `log3 + log2·ℤ`" reading is *not* Feller-lattice; it needs the origin, which
  incommensurability blocks.) Feller's renewal theorem: a non-lattice step gives **no** archimedean log-periodic
  oscillation.
- **Yet the period-~9 oscillation is measured ⟹ it is not archimedean in origin.** It comes from the **3-adic level
  structure** — integer levels `r`, trivially lattice. **Lapidus–Hùng: every p-adic self-similar string is lattice**
  (ζ rational, poles periodic), so the oscillation is **structurally forced by the p-adic setting and forbidden in
  the archimedean analogue.** Its period is set by the level-transfer operator's `arg(λ₂)`, **not** by `2π/log2`
  (that = 9.0647 is suggestive vs measured ~9.2–9.5 but **unearned numerology** like F1-D until a mechanism ties
  `arg(λ₂)` to it).

## Status
**F1: the period-~9 mode is real, dominant (0.984 not ½), non-integer, and p-adic-lattice-forced.** **A** — the exact
ladder confirms rate 0.984 to k=12 (`|ε|·2^k → 9.3`, `|ε|/0.984^k` bounded, `Λ_r` flat at r≥7). **B** — every integer
period 2–9 is excluded by the exact 11-term sign sequence ⟹ **the period is not an integer**, consistent with an
irrational value near 9. **C** — the mode lives in the oscillating intermediate lift-tower entries, not the rank-1
η=1 part, and is not a clean finite eigenvalue. **D** — `√0.97008 = 0.98493 ≈ 0.984` is unearned numerology. **E** —
the archimedean renewal is non-lattice (`log3/log2` irrational, Wilson correct), so the oscillation is **forbidden
archimedean-ly and forced p-adically** (Lapidus–Hùng): period-9 is a **signature of the correct (p-adic) category**,
not an empirical nuisance.

**Consequence for the crux (owed to the pen).** The frontier is now correctly named: the theorem's asymptotic rate is
a **p-adic-lattice-forced log-periodic mode of the level-transfer operator**, non-integer period near 9, dominant
over the subcritical ½. The remaining analytic step (`Σ|A_r(m)| < ∞`) only needs rate `< 1` — which 0.984 satisfies —
so the theorem is unaffected; what is *identified* is the exact character of the dominant mode. Machinery: the
lattice/nonlattice dichotomy is Lapidus–van Frankenhuijsen (Thm 3.6) + Feller §XIII; the p-adic-always-lattice is
Lapidus–Hùng; the large-amplitude log-periodic regime (ours ~0.2–4.6%, vs 10⁻⁵ for Ising on hierarchical lattices)
is Sornette / Luck (arXiv:2403.00432), where the oscillation is dominant not a correction — consistent with 0.984
beating ½. No fitting; exact ladder + exact sign-period exclusions; the 2π/log2 and √0.97008 coincidences flagged as
unearned; F1-E corrected from the coset-misreading to the Feller/Lapidus split.
