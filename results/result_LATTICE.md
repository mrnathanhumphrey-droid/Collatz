# RESULT — LATTICE: the tower is NOT a base-3 Lapidus lattice string; the oscillation is base-2, incommensurate with the base-3 tower (2026-07-27)

**Probe:** `probe_lattice.py`. Is the Syracuse tower a self-similar p-adic fractal string in Lapidus's exact sense
(⟹ rational zeta in 3⁻ˢ, lattice, log-period exactly log3, oscillation bounded/periodic, no crossing by theorem)?
Certified T_i machinery only (T_1..14 recomputed, T_15..20 banked/gated). **NOT** an operator diagonalization (that was
the R29 mush) — the zeta is tested from the scaling ratios.

## L-A — no bare r=1/3 string; no clean finite-rank (rational) zeta at this resolution
- **Λ per-level ratio ≈ 0.88, not 1/3.** `log(0.88)/log3 = −0.114` — not a recognizable rational, and nowhere near
  the `r=1/3` a base-3 string would require. The scaling ratios do **not** close into powers of 1/3 (nor drift toward
  them — they settle at ≈0.865–0.88, away from 1/3).
- **Hankel-rank test (rational zeta ⟺ finite Hankel rank of Λ_i):** normalized singular values on i=8..20 =
  `1.0, 0.060, 0.018, 0.010, 0.0061, 0.0012` — **smooth (≈geometric) decay, no gap/cutoff.** One dominant mode (the ρ₁
  decay) + a graded tail; no clean finite mode set ⟹ **no clean rational-in-3⁻ˢ zeta identifiable** from 13 points.

## L-B — THE DECISIVE FACT: a base-3 lattice mode is INVISIBLE at integer levels; the observed period is base-2
A self-similar string with scale 3 (each level ×3 in the modulus) has complex dimensions `s = D + 2πin/log3`, so its
oscillatory term is `x^{2πin/log3}` with `x = 3ⁱ`, i.e. `e^{2πin·i·log3/log3} = e^{2πin·i} = 1` for **every integer level
i**. **A bare base-3 lattice oscillation is aliased to a constant — period exactly 1 level, structurally invisible.**
Therefore the measured **~9-level** oscillation in Λ_i **cannot be the base-3 lattice fundamental.** It matches
**`2π/log2 = 9.065`** — the base-2 `÷2ᵛ` multiplier — and **`log3/log2 = 1.5850` is irrational**, so the base-2 oscillation
sampled at base-3 level spacing is **quasi-periodic (an irrational rotation), not periodic.** (This reconfirms F1-E's
archimedean non-lattice reading, now directly in the level-index oscillation.)
- ω-fit corroboration (WEAK, 13 pts): full i=8..20 → period 13.8; early i=8..14 → 6.8; late i=14..20 → 7.2. **ω drifts**
  and **none locks to 5.72** (the base-3 lattice value). Per the pre-registered caveat, a *stable* ω would have been strong
  evidence for lattice; the observed *drift* is only weak evidence for non-lattice — so the load-bearing argument is the
  aliasing fact above, not the fit.

## L-C — dimension not cleanly recognizable
`ρ₁ ≈ 0.891/level ⟹ ρ₁ = 3⁻ᴰ ⟹ D = 0.106`. Nearest simple values 1/8=0.125, 2/15=0.133, 1/7=0.143 — none close;
D is not a recognizable rational at current precision (and ρ₁ itself is still drifting down, so D is not yet stable).

## L-D — the boundedness corollary does NOT get the lattice upgrade (but no-crossing survives empirically)
Since the tower is **not** a base-3 lattice, i=20's "no crossing" does **not** upgrade to a lattice theorem via this route.
However the oscillation **is empirically bounded** (detrended log-residual stays in ≈[−0.20, +0.12] across i=8..20), so once
the real ρ₁ decay dominates — which it does, Λ positive and decreasing through i=20 — no sign change occurs; the no-crossing
conclusion holds **empirically**, just not as a lattice theorem. A structural proof of no-crossing would need a different
lever (an amplitude bound on the oscillation, or the C_ρ positivity), not the Lapidus lattice.

## VERDICT — NON-LATTICE by two incommensurate scales; irrational-rotation reading is the live one
The Syracuse tower is **not** a Lapidus base-3 self-similar string. The obstruction is not "approximate self-similarity"
but **two incommensurate scales**: the tower is base-3 (whose lattice mode is invisible at integer levels) while the
oscillation is base-2 (from `÷2ᵛ`), and `log3/log2` is irrational ⟹ the tail is **quasi-periodic, not periodic**, with
rotation number tied to `log3/log2`. So of Wilson's three outcomes, **"non-lattice / drifting ω → the irrational-rotation
reading is live; the 'period' was never a period"** is the one that fires. The clean "log-period = log3, no-crossing by
theorem" does **not** hold.

**Crucial caveat (from the parallel lit hunt, Niven/Lapidus):** non-lattice / quasi-periodic does **NOT** imply S∞ is
irrational. An algebraic subdominant eigenvalue with irrational rotation number sums to a closed form (`½·(3+4i)/5` is the
counterexample). So this refutes the clean log-3 lattice **theorem** but does **not** establish that 7/15 is "wrong in kind."
The rational-vs-irrational verdict is decided by the Pisot/digit structure of ν (= C_ρ), not by the rotation number.

**Not at stake:** S∞ ≈ 0.475 (value stands), the i=20 no-crossing observation (holds empirically), P6D–P6K identities,
S_{i+1}=2T_i. **Newly established:** the log-oscillation is base-2 / incommensurate with the base-3 tower (non-lattice);
Lapidus's clean lattice theorem is not available; the "unpinnable period" is genuinely not a fixed period.
