# RESULT — P6: the phase/cascade has a domain source — d₁ is EXACTLY the cross-parity interference (2026-07-26)

**Probe:** `probe_p6.py`. Where in the 2-adic domain does arg π̂ (= the phase-carried cascade, P1LVL) come from?
Conditioning on the first L branch variables = L-fold unrolled SINGLEREC (certified). Per level, no new transport.

## P6-A — GATE passes
Depth-L cylinder sum reproduces π̂_j(ξ) to **1e-16..1e-14** for L=1,2,3 at j=4,6,8 (15 sample primitive ξ). The
unrolling is correctly indexed.

## P6-B — coherence: NOT diffuse (a source exists), = the P_k ⟨2⟩-orbit carrier
Coherence `|Σ_a term_a| / Σ_a|term_a|` over primitive ξ: **mean ~0.55 at every level** (0.555, 0.550, 0.556 for
j=4,6,8; median ~0.55–0.61, min 0.007–0.16, max 0.93–0.98). Argmax modulus is **a=1 for 83–89% of ξ** (small branch
variables dominate). This is *not* the diffuse-zero failure mode — the phase HAS a source — and coherence ~0.55 is
exactly LAMBDA's P_k deterministic carrier (median ~0.57), i.e. the ⟨2⟩-orbit partial alignment. Known structure, not
new.

## P6-C/D — THE MECHANISM: d₁ is exactly the CROSS-PARITY cross-term
Recomputing the shell A_j(1) with the first branch variable a restricted, then bridged:
| subset | j=4 | j=6 | j=8 |
|--------|-----|-----|-----|
| **all** | +4.01e-3 | +3.61e-3 | +2.17e-3 |
| **a even only** | **+1.6e-16** | **+6.9e-17** | **+1.7e-15** |
| **a odd only** | **+1.4e-16** | **−2.1e-15** | **+7.4e-15** |
| **even×odd cross** | **+4.01e-3** | **+3.61e-3** | **+2.17e-3** |

**Same-parity cylinders contribute EXACTLY ZERO to the m=1 channel; the entire shell is the even×odd cross-term** (to
machine precision, all three levels). The parity a mod 2 IS the ⟨2⟩/⟨4⟩ = ℤ/2 coset = the `(−1)^a` quadratic-character
structure that BRIDGE2 identified. **So `d₁ = ⟨cross-parity interference⟩`, exactly, and `d₁ > 0` ⟺ that interference
is positive.**

⚠️ **REASON RETRACTED (see result_P6B):** I originally asserted the reason was "m=1 is odd-lag; lag-1 connects opposite
parities." That is WRONG — the P6B parity-decomposition of A_j(m) shows the rule is **`3∤m → cross-parity, 3|m →
same-parity`** (m=2 is cross, m=3 is same), NOT odd/even lag. d₁ (m=1, 3∤1) is cross-parity so this *finding* stands,
but the odd/even-lag corollary below is FALSE; the correct law is 3∤m, which aligns with the enriched/depleted
dichotomy. The mechanism's derivation is Wilson's pen.

Finer structure (P6-C, a mod 3 — the 3-adic slice of a mod ord): the residue splits into **opposite-sign families**
(j=4: a≡1→−0.0029, a≡2→+0.0025; j=6,8: a≡1→+, a≡2→−, sign flips with j = period-2), a near-cancellation whose net is
the small positive shell. Small a (a≤3) carries net +, large a (a≥4) net small −.

## Pre-registered expectation — CONFIRMED (two primes, via mod 2 and mod 3)
Wilson: structure should live on `a mod 2·3^{j−1} = ord_{3^j}(2)`, not on a. The single-residue-mod-ord test (j=4,
ord=54) is degenerate (each residue too sparse, ~0). But the structure IS on a mod ord — decisively via its **ℤ/2
(parity) factor** (the exact cross-parity mechanism above) and its **3-adic (a mod 3) factor** (the sign-split). The
two primes meet exactly as predicted: the ⟨2⟩ order-2 part carries d₁ (cross-parity), the 3-adic part sets the sign.

## Verdict — P6 closes with a mechanism, not a null
The phase/cascade is NOT diffuse (coherence ~0.55, a source exists) and its m=1 content has an **exact domain-side
form: the cross-parity cylinder interference** (same-parity = 0). This localizes d₁ to the `(−1)^a` = ⟨2⟩/⟨4⟩ = the
two-primes ℤ/2 structure — the same parity BRIDGE2's ℤ/2 flagged and the arc's carry/parity "alternating coherence"
(FOURCELL, CARRYLEMMA) circled without pinning. **Pen target, now exact:** `d₁ = ⟨even×odd cross-term⟩ > 0`; the
governing law (from result_P6B) is `A_j(m) cross-parity ⟺ 3∤m` (= the enriched/depleted dichotomy), this is the
sharpest domain-side handle on the cascade the arc has produced. Sits under Wilson's `7/15 ⟺ Σ_{i≥2}Λ_i = −1/210`.
Not at stake: P1LVL, BRIDGE2, P4, CHANNEL_ID, MEAN1, R1–R30. Cheap (13s).
