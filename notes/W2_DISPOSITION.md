# W2.E — Disposition: did monotone cumulants close `−1/30`?

**Date:** 2026-05-14
**Task:** Track A wrinkle 2, step 5 (final). One-page disposition of whether
the monotone-cumulant framework supplies the closed-form `−1/30` coefficient,
and what is open.
**Mode E:** verbatim where citing prior work.

---

## 1. Verdict

**The monotone-cumulant framework does NOT supply a closed-form derivation
of `−1/30` from first principles.** It supplies:

- **Mechanism** (κ_2^B subdominant + per-step B-additivity): rigorous, modulo
  Mode-E gap on B-amalgamated lift of HS Thm 3.26 (Wrinkle 1).
- **Rate** `(1/2)^n`: comes from κ_2^B's B-spectral structure, specifically
  from the rate-(1/2) eigenvalue of the bilinear pair operator T_M (R77 §3
  empirical, R76 §6 / R77 §6 open rigorously).
- **Sign**: negative (S_n → 7/15 from below for n ≥ 3), rigorous via R76
  Thm 76.3 `−2` factor and the empirical sign of `R_n − R_∞`.

The **amplitude `1/30`** is **not** a Hasebe combinatorial constant. The
naive expectation that the count of one-2-block monotone partitions on [n]
gives a constant 14 = 2·7 is **incorrect**: the actual count is
`n · H_{n−1} − (n−1) ~ n · ln(n)`, which grows polynomial-times-log in n,
not as a constant.

---

## 2. What W2.A-W2.D established

| Deliverable | Result |
|---|---|
| W2.A | κ_2^B(Off_j, Off_j)|_{(1,4)} has leading bilinear coupling weight 1/8 from (v=1,v'=3)∪(v=3,v'=1). Cross-step κ_2^B = 0 structurally. Closed-form rational amplitude requires T_M (open). |
| W2.B | Monotone partitions on [n] with one 2-block + (n−2) singletons count is `n·H_{n−1} − (n−1)`, polynomial-times-log in n. The combinatorial factor is NOT a constant 14. |
| W2.C | `1/30 = 1/(2·15)` rigorous decomposition: `2` from R76 Thm 76.3 bilinear pair, `15 = 3·5` from R75 Plancherel + R77 T_diag. The `14 = 2·7` decomposition is empirical fingerprint (the `7` cancels through `S_∞ = 7/15`). |
| W2.D | Closed-form structural derivation: `ε_n = −2α·(1/2)^n` with `α = 1/60 = 1/(4·15)`, where `4 = 2² = R64.B class-mass ratio` and `15 = R75+R77 normalization`. The closed-form `α = 1/60` reduces to the T_M eigenvector amplitude being exactly 1 (empirical, not rigorous). |

---

## 3. What's rigorous

- R76 Thm 76.3 `S_n = −2·R_n` (verbatim, rigorous)
- R75 Plancherel S_n = Σ |μ̂_n(ξ)|² (verbatim, rigorous)
- R77 T_diag = (1/5)·[[1,1],[4,4]] eigenstructure on (1, 4) (verbatim, rigorous)
- Hasebe Defn 3.23 monotone-partition counting on [n] (verbatim, rigorous)
- The rate `(1/2)^n` IS the κ_2^B subdominant in the monotone framework
  (mechanism rigorous modulo Wrinkle 1 B-lift)

## 4. What's open

- T_M spectral identification at finite k (R75 §8, R76 §6, R77 §6) — the
  rigorous derivation of λ_2 = 1/2 and the eigenvector amplitude `1/60`.
- B-amalgamated lift of HS 2011 Thm 3.26 (Wrinkle 1 in MONOTONE_CLOSURE_WRITEUP).
- Closed-form rational amplitude of κ_2^B(Off_j)|_{(1,4)} — depends on T_M.

---

## 5. Match to empirical

`|ε_n|·2^n` plateau (n=2..6): 0.0381, 0.0407, 0.0392, 0.0369, 0.0349
(c_seven_forty_fifth.md line 129; tail values from W2.D computation).

Predicted `1/30 ≈ 0.03333`. The empirical plateau decreases toward `1/30`
from above as `n` grows, consistent with `−1/30·(1/2)^n` plus a positive
O((1/4)^n) correction.

**The derived closed-form `−1/30·(1/2)^n` matches the empirical R77 §4 fit
to within the documented O((1/4)^n) tail.**

(k=7, 8 deviations 0.150, 0.191 reflect the multi-spectral transient onset
noted in R77 §6 — this is Wrinkle 3 territory, not W2.)

---

## 6. Headline

**Did monotone cumulants supply the closed-form `−1/30`?**

**No, not rigorously.** The framework supplies:
- the correct rate `(1/2)^n` (via κ_2^B B-spectrum, mechanism rigorous),
- the correct sign (via R76 Thm 76.3 `−2` and the B-spectrum sign),
- the structural decomposition `1/30 = 1/(2·15)` with `2` rigorous from R76
  Thm 76.3 and `15 = 3·5` rigorous from R75 + R77.

But the closed-form AMPLITUDE `α = 1/60` of `R_n − R_∞` (equivalent to the
eigenvector amplitude of T_M at eigenvalue 1/2) is **not** a Hasebe
combinatorial output. It is the same outstanding T_M spectral calculation
flagged in R75 §8, R76 §6, R77 §6.

**Wrinkle 2 is therefore reframed: the monotone-cumulant framework reduces
Wrinkle 2 to the SAME T_M spectral analysis as the rigorous `(1/2)^n` rate
derivation (Wrinkle 2's sibling).** Closing T_M closes both. The
combinatorial 14 = 2·7 is not the right object — the right object is the
T_M eigenvector amplitude `1`.

---

## 7. Implication for the broader closure program

Per MONOTONE_CLOSURE_WRITEUP.md §2 Wrinkle 2, the effort estimate was
"4-8 hours focused" with risk "medium". The actual finding:

- The work CAN be done in 4-8 hours — done, this writeup.
- The CONCLUSION is that the closed-form `1/30` is not derivable from
  monotone-cumulant combinatorics alone; it requires T_M spectral analysis
  (a separate open step).
- This **changes the disposition of Wrinkle 2** from "dig-hard tractable
  combinatorial calculation" to "redirects to the T_M spectral step
  (R75 §8 / R76 §6 / R77 §6 outstanding)".
- The T_M step is **separately tractable** (R75 §8 estimates "another
  session"), so the closure remains within the dig-hard regime, but the
  path is via T_M, not via Hasebe Defn 3.23.

**Updated recommendation.** Reframe Wrinkle 2 as: "close T_M spectral
identification (R75 §8 / R76 §6 / R77 §6)". This is a separate analytical
problem from monotone-cumulant combinatorics, and it is the load-bearing
step for closing the full subdominant `−1/30·(1/2)^n` statement.

---

## 8. Mode-E gaps in W2

1. **Closed-form rational amplitude of κ_2^B(Off_j)|_{(1,4)}**: requires
   T_M (W2.A §6, R76 §6 open).
2. **B-amalgamated lift of HS Thm 3.26**: same as Wrinkle 1 (MONOTONE_CLOSURE_
   WRITEUP §2.1).
3. **Mixed monotone cumulants**: HS 2011 §6 last paragraph references
   "T. Hasebe and H. Saigo, in preparation" but follow-up status not
   resolved in corpus (Mode-E gap #3 in Deliverable A).
4. **κ_2^B B-spectral structure**: the leading-bilinear coupling
   (v=1, v'=3) ∪ (v=3, v'=1) with weight 1/8 is identified, but the full
   B-spectral content of κ_2^B requires T_M (= same T_M as the rate).

---

## 9. Files

- W2_KAPPA2_CALC.md
- W2_PARTITION_COUNT.md
- W2_PLANCHEREL_NORM.md
- W2_CLOSED_FORM.md
- MONOTONE_CLOSURE_WRITEUP.md (parent document; Wrinkle 2 disposition
  should be updated to reflect this redirect to T_M)
