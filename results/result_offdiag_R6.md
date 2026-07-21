# Probe R6 — the interference ledger — **R6-A GATE PASS: OffDiag₂ = −4/21 derived, not read**

**Date:** 2026-07-21  Exact rationals, minutes. Probe `probes/probe_offdiag_R6.py`. Gates Wilson's off-diagonal
channel formalism: the shell is the primitive-frequency Plancherel mass S_k = Σ_{ξ prim mod 3^k}|μ̂_k(ξ)|² =
𝔼_{X,X′}[c_{3^k}(X−X′)] (c = Ramanujan sum), split into diagonal (v=v′, replicates S_{k−1}) + off-diagonal.

## R6-A — the formalism's first gate: **PASS** (independent of the S-table)
At k=2: X = 2^{−v}(1+3Y) mod 9, v ~ Geom(½) (2^{−v} mod 9 has period 6), Y ~ π₁ (masses 1/3, 2/3 on residues
1,2 — Judge One). Ramanujan c₉ = 6 (9|Δ), −3 (3|Δ, 9∤Δ), 0 (coprime). Split by v=v′ vs v≠v′:

| channel | value | expected |
|---|---|---|
| **diagonal (v=v′)** | **2/3** | S₁ = 2/3 — **REPLICATES** |
| **off-diagonal (v≠v′)** | **−4/21** | −4/21 — **GATE PASS** |
| sum | 10/21 | S₂ ✓ |

The diagonal weight = Σ_j D_j = **1/3 exactly** (P(v=v′) = Σ 4^{−v}, infinite valuations — no folding artifact),
and the primitive lift mod 3 → mod 9 is **3-to-1**, so **3 × ⅓ = 1**: the diagonal channel replicates the previous
shell exactly. **OffDiag₂ = −4/21 is now derived from the v≠v′ character sums, not read off the S-table** — the
new formalism's first gate is green. The boundary's deepest mechanism, literally: (3-adic lift multiplicity 3) ×
(2-adic self-collision weight ⅓) = 1. Two primes, two factors.

## R6-B — the interference ledger (the constant's target)
OffDiag_k = S_k − S_{k−1} (frozen S from `Basic.lean`):

| k | OffDiag_k | float | running total | signs |
|---|---|---|---|---|
| 2 | **−4/21** | −0.190476 | −0.190476 | − |
| 3 | −2980/203889 | −0.014616 | −0.205092 | − |
| 4 | +5699915795296300/2159281421340253987 | +0.002640 | −0.202452 | + |
| 5 | +(…) | +0.001301 | −0.201152 | + |

**Σ_{k≥2} OffDiag = S∞ − S₁ = −1/5** (the constant's derivation target — the total off-diagonal ledger equals
minus one fifth). OffDiag₂ = −4/21 is **95.2%** of it; Σ_{k≥3} = −1/105. Sign pattern **−, −, +, +** — the k≥4
tail is net **positive** (S∞ − S₃ = +5191/1019445), the exact overshoot R5-B measured, now located inside the
interference channel.

## R6-C — channel anatomy by |v−v′| (the two-sign raw material)
Off-diagonal decomposed by valuation gap, exact (Ramanujan c₉ at k=2, c₂₇ at k=3):

| \|v−v′\| | k=2 | k=3 |
|---|---|---|
| **1** | **0** (parity selection) | **0** (parity selection) |
| **2** | **−1/6** | **+2/147** (sign FLIP) |
| **≥3** | −1/42 | −0.02822 |
| total | −4/21 | −2980/203889 |

**Two structural facts for the tail's two-sign law:**
1. **|v−v′|=1 vanishes identically** (both k) — odd valuation gaps flip the mod-3 class (ord₃(2)=2), killing the
   character sum. The interference lives entirely in **even and ≥3 gaps.**
2. **The |v−v′|=2 channel carries the sign flip**: −1/6 at k=2 → **+2/147 at k=3.** As k grows this channel turns
   positive while |v−v′|≥3 stays negative — the competition that produces the overshoot (OffDiag negative for
   k=2,3, then net positive for k≥4). The diagonal channel replicates **exactly at both k=2 (→2/3) and k=3
   (→10/21)** — verified, not assumed.

## Status
**R6-A GATE PASS** — the off-diagonal formalism derives OffDiag₂ = −4/21 directly from the v≠v′ Ramanujan sums,
with the diagonal channel replicating S_{k−1} exactly (3×⅓=1) at k=2 and k=3. **R6-B**: the ledger telescopes to
Σ OffDiag = −1/5 (the constant's target), first term 95.2%. **R6-C**: |v−v′|=1 dead (parity), |v−v′|=2 flips
sign k=2→k=3 (−1/6 → +2/147), |v−v′|≥3 stays negative — the raw material for the pen's two-sign tail law. The
element ⟨1|Ñ|ψ⟩ and the ledger Σ OffDiag = −1/5 are one object; the ledger costume is finitely computable term by
term, forever, and its k=2 term is now derived. No fitting; exact rationals.
