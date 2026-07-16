# Result 45 (G1, qx+1 paper) — the q=3 boundary echo/drip decomposition: the CONSTANT-COEFFICIENT recursion (P≡2/3, I≡7/45) is FALSIFIED as stated. What survives is the 7/15 LIMIT and the constant-pressure diagonal tap. The 2/3 as a limiting coefficient is neither confirmed nor cleanly refuted — the residue binning is the wrong split (R19), the exact residual is circular.

**Date:** 2026-07-16. **Gate:** G1 of the L3 Session-One worksheet (desktop agent's), redesigned mid-flight: the total ω(3)=1 is CIRCULAR (banked six ways: R15 slope, R16 D-ratio 1.00092, R22 ρ→1, R25 λ₂→λ₁…), so test the PARTS — P_k (echo) → 2/3 [at risk], I_k (drip) → 7/45 — per level.

**Verdicts: [PART 1, exact arithmetic] constant-coefficient recursion FALSIFIED (residuals drift, super-geometric — R23 stands) / c\*=7/15 SURVIVES as limit / [PART 2, operator] constant-pressure diagonal tap CONFIRMED (D=(1/3)^k exact) / the 2/3 NOT isolated — residue binning ≠ deviation split (R19 subtlety, pre-flagged).**

**Headline: the redesign earned its keep by catching two traps the total-only gate would have hidden. (1) The exact rational residual `c_{m+1} − (2/3)c_m` is NOT the constant 7/45 the hand-derivation needs — it drifts [0.1441, 0.1565, 0.1560], so the constant-coefficient form is only asymptotic (outcome 2, super-geometric, R23 confirmed), NOT exact (outcome 1). (2) That residual landing near 7/45 is CIRCULAR: `c_{m+1}−(2/3)c_m → c*·(1/3) = (7/15)/3 = 7/45` automatically for ANY sequence converging to 7/15, independent of the true propagation coefficient — so the exact test can kill the constant form (it did) but cannot certify the 2/3. The operator provenance, meant to measure the 2/3 directly, used the pair-RESIDUE binning (a==b vs a≠b), whose off-diagonal mass (0.127 at k=2) is NOT the deviation ‖d₂‖²=10/189=0.0529 (2.4× off) — the recursion's 2/3 lives in the DEVIATION split, not the residue one. So the 2/3 is neither confirmed nor cleanly refuted; the binning that would isolate it is the deviation-propagation operator L_k (the incomplete piece of `c_seven_forty_fifth_derivation.py`). What DID survive, exactly: c\*=7/15 as a limit, and the constant-pressure diagonal tap — under binning {a==b ∧ γ=0} the reservoir mass is EXACTLY (1/3)^k at every level, the drip source never depletes.**

Probe: `probe_45_G1_provenance_echo_drip.py`. Log: `result_45_G1_provenance_log.txt`. Runtime: seconds (q=3, L≤3 sparse; exact rationals k≤5).

## PART 1 — exact rational recursion residual (no binning, no operator)

`c_m := ‖d_m‖²·3^m`, with `‖d_m‖² = ‖π_m‖² − (1/3)‖π_{m−1}‖²` (R74). Exact big-int rationals to k=5.

| m | c_m (exact) | c_m decimal | → 7/15 |
|---|---|---|---|
| 2 | 10/21 | 0.476190 | 0.466667 |
| 3 | 31370/67963 | 0.461575 | |
| 4 | 143195649659456490/308468774477179141 | 0.464214 | |
| 5 | …/… | 0.465515 | |

c_m **oscillates and converges toward 7/15 = 0.466667** ⇒ **c\* = 7/15 survives as a limit.**

Constant-coefficient test — is `c_{m+1} − (2/3)c_m` exactly 7/45 = 0.155556?

| m | c_{m+1} − (2/3)c_m | vs 7/45 | constant? |
|---|---|---|---|
| 2 | 88150/611667 = 0.144114 | −1.14e−2 | NO |
| 3 | …/… = 0.156498 | +9.42e−4 | NO |
| 4 | …/… = 0.156039 | +4.83e−4 | NO |

Residuals **[0.144114, 0.156498, 0.156039] DRIFT** — not the flat 7/45 an exact constant-coefficient recursion requires. **Outcome (2): super-geometric, asymptotic-only. R23's exact super-geometric finding stands.**

**⚠️ Circularity caught:** at the fixed point `c_{m+1} − (2/3)c_m → c*(1 − 2/3) = c*/3 = (7/15)/3 = 7/45` for **any** sequence → 7/15, regardless of the true P. So the residual → 7/45 **re-encodes c\*=7/15; it does NOT independently confirm P=2/3.** The exact test kills the constant form but cannot certify the coefficient — exactly the "wrong in two compensating ways" failure mode the redesign was built to block.

## PART 2 — operator provenance (build_M, gate-validated: Σ(M^k v0)=‖π_k‖²)

Echo = O→O flow (already off-diagonal, propagated); drip = D→O flow (diagonal reservoir → fresh off-diagonal). P_k = echo / O_prev.

| binning | k | O_mass | D_mass | (1/3)^k | echo | drip | **P_k** |
|---|---|---|---|---|---|---|---|
| A {a==b} | 2 | 0.126926 | 0.111169 | 0.111111 | 0.052853 | 0.074074 | **0.237840** |
| A {a==b} | 3 | 0.059029 | 0.037432 | 0.037037 | 0.034325 | 0.024704 | **0.270432** |
| B {a==b ∧ γ=0} | 2 | 0.126982 | 0.111113 | 0.111111 | 0.052909 | 0.074074 | **0.238092** |
| B {a==b ∧ γ=0} | 3 | 0.059423 | 0.037038 | 0.037037 | 0.034731 | 0.024691 | **0.273511** |

**(a) Constant-pressure diagonal tap — CONFIRMED.** Under binning B (a==b ∧ γ=0) the diagonal reservoir mass is **exactly (1/3)^k** (0.037038 = (1/3)³ = 0.037037 to rounding; k=2 → 0.1111 = (1/3)²). The drip source is a constant-pressure tap that never depletes — the mechanism by which injection would never slow at q=3. This piece of the drip picture is real and exact.

**(b) The 2/3 is NOT there — binning miss, not carry-error verdict.** The residue off-diagonal `O_mass` (0.127 at k=2) is **not** the deviation `‖d₂‖² = 0.0529` (2.4× off — verified via the log's `O_mass vs ‖d_k‖²` line). The recursion's 2/3 lives in the **deviation** decomposition `‖π_k‖² = (1/3)‖π_{k−1}‖² + ‖d_k‖²`; the operator here measured the **pair-residue** decomposition `‖π_k‖² = (1/3)^k [pure diagonal] + O_mass [rest]`. Different splits. The residue-echo P_k is **0.238 → 0.274, rising, nowhere near 0.667.** Since it is the wrong split, this neither confirms nor kills the 2/3 — it is R19's pre-flagged subtlety ("colliders concentrate on 2-coord-differing pairs; the two-bin story is subtler") made concrete.

## Verdict on G1

- **SURVIVES:** c\* = 7/15 (limit of c_m); the constant-pressure diagonal tap (D = (1/3)^k exact).
- **FALSIFIED as stated:** the constant-coefficient recursion P_k ≡ 2/3, I_k ≡ 7/45 per level — the exact residuals drift (super-geometric; R23 stands). 2/3 and 7/45 are limits at best, not exact per-level coefficients.
- **UNRESOLVED (correctly, not faked):** the 2/3 as a limiting propagation coefficient. The exact residual is circular (auto-follows from c\*=7/15); the residue-binning operator measures the wrong split (O_mass ≠ ‖d‖²). Neither instrument isolates it.
- **The correct instrument, identified:** the deviation-propagation operator **L_k** — propagate the deviation vector d_k through the lifting operator and read *its* echo coefficient. That is the honest home of the 2/3 (the incomplete construction in `c_seven_forty_fifth_derivation.py`). Session-Two / G1b.

Three pre-registered outcomes, mapped: PART 1 → **(2)** (super-geometric, R23 stands) on the recursion FORM; the 2/3 COEFFICIENT itself lands in none of (1)/(2)/(3) cleanly because both instruments were the wrong lens for it — the honest fourth outcome the spec's own §"binning under test" anticipated.

## Not at stake
R1–R44. This tests the echo/drip split of the q=3 boundary; it changes no r_q value and does not touch L3's statement. The L3 target (r_q<1 for d≥3) is untouched; G1 was the d=2 entrance exam, and it reports the boundary's structure honestly rather than confirming a hand value the corpus already owns.

_Reporting discipline: the redesign (test PARTS not the banked total) caught (i) a super-geometric drift the total-only gate hides, and (ii) a circular residual that would have manufactured a false "P=2/3 confirmed." The a==b/a≠b binning's mismatch with ‖d‖² is disclosed (log line `O_mass vs ‖d_k‖²`), and the 2/3 is left UNRESOLVED rather than force-fit to the residue-echo 0.25 or the circular 7/45. Exact rationals throughout PART 1; gate-validated operator (Σ=‖π_k‖²) throughout PART 2. No fit._
