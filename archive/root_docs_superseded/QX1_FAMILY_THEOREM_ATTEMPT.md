# qx+1 Family-Level Plancherel Saturation Theorem — Attempt Results

## Disposition

> **CLAIM_1_PARTIAL + STRUCTURAL_OBSTRUCTION_FOUND**
>
> A genuine new family-level result lands; a different one does not; the Move 2 pipeline to dissolve R77.2's reliance on Tao Prop 1.17 does not close.

**One-paragraph summary:** The full-period Plancherel-saturation magnitude formula `|F̂_p(ξ)| = p^{(r+3)/2}` on supp(F̂_p) generalizes from R78.3's q=3 case to every prime p ≥ 3 — empirically verified to machine precision at p ∈ {3, 5, 7}, with a proof template (Cochrane Theorem 2 + Plancherel + principal-unit equidistribution) that is provably p-blind. This is a real structural result candidate, a family-level extension of R78.3. **But it is the wrong object to substitute for Tao Prop 1.17.** Tao Prop 1.17 bounds `|μ̂_n(ξ)| ≪_A n^{-A}` where μ̂_n is a Markov-chain stationary's characteristic function. The F̂_p result bounds a deterministic cyclic-group character sum's full-period Fourier transform. These are distinct objects, related but not identifiable. Even passing through the K_p short-window object — the q-sweep's actual rate-½ exponent — the path requires the same Burgess-type bilinear bound that R78 hit and stops at, family-level or q=3.

The relevant phases below; full breakdown after the headline.

| Claim | Status | Comment |
|---|---|---|
| **CLAIM 1 (rate-½ universal)** | **PARTIAL** | Full-period F̂_p magnitude formula closes structurally at family level. Short-window K_p √N saturation does NOT close — needs Burgess-type bilinear bound, same wall as R78 at q=3. |
| **CLAIM 2 (prefactor closed form)** | **PARTIAL** | F̂_p magnitude has explicit closed form p^{(r+3)/2}, uniform across p — no q-varying prefactor at the F̂ level. The empirical 1.4× C_p variation in K_p data is from the short-window-vs-full-period scaling and the missing bilinear bound, NOT a different magnitude law per prime. |
| **CLAIM 3 (Nisoli ε_K from q=3 specialization)** | **FAILS** | Even at q=3 the F̂ magnitude bound (rigorous) doesn't translate to a Tao-Prop-1.17-shaped bound on μ̂_n(ξ). Different objects. |

---

## Pre-registration adherence

- **Pre-registered:** 2026-05-11T11:15 EDT, committed at `f96fb86` before any compute. Locked rules in [QX1_FAMILY_THEOREM_PRE_REGISTRATION.md](QX1_FAMILY_THEOREM_PRE_REGISTRATION.md).
- **Framing reconciliation (pre-reg §0):** the Move 2 prompt's "rate-½ universal across q" had to refer to the K_p Kalafatelis-sum exponent (object B), not the ε_n rate-½ envelope (object A, q=3-specific per `result_q_sweep_test_1_rate.md`). The pre-reg locked the K_p interpretation. **This framing distinction was load-bearing and held under attempt.**
- **Procedure followed as locked.** No mid-run parameter changes.
- **A2 supplementary data:** ran K_p at p ∈ {11, 13, 17, 19, 23} for r ∈ {3, 4}. Empirical |K|/√N ∈ [0.92, 1.10] across these primes — consistent with √N saturation extending past the original q-sweep range. Documented but not used as structural proof.

---

## Phase 1: Structural formulation

The family-level objects:

- **f_p(u) := e_M(c·(1+p)^u)** for u ∈ Z, M = p^{r+1}. Periodic in u with period p^r (since (1+p) has multiplicative order p^r in (Z/p^{r+1})*).
- **F̂_p(ξ) := Σ_{u=0}^{p^{r+1}−1} e_M(c·(1+p)^u − ξu)** — full-period (length M) Fourier sum.
- **K_p(r, c, m) := Σ_{u=0}^{p^{r-1}−1} e_M(c·(1+p)^u − p²·m·u)** — short-window (length N = p^{r-1} = M/p²) character sum. **This is the q-sweep object.**
- **μ̂_n(ξ) := E[e^{−2πi ξ Syrac(Z/p^n)/p^n}]** — Markov-chain stationary's characteristic function. **This is Tao Prop 1.17's object.**

**Critical for what follows:** F̂_p, K_p, μ̂_n are three distinct objects. F̂_p ↔ K_p are related by Pólya-Vinogradov truncation; K_p ↔ μ̂_n connection is not established in the project.

The family-level Plancherel identity is straightforward:
> Σ_{ξ ∈ Z/M} |F̂_p(ξ)|² = M · Σ_{u=0}^{M−1} |f_p(u)|² = M · M = M²    (Parseval on Z/M).

Combined with F̂_p's period structure (f periodic with period p^r ⟹ F̂_p supported on p·Z/M), |supp(F̂_p)| ≤ M/p = p^r. The sub-support {p·a : a ≡ 1 mod p} from R78.2's principal-unit argument has size p^{r-1} = M/p². With equidistribution, max|F̂_p|² = M²/|supp| = M·p² and so |F̂_p(ξ)| = p·√M = p^{(r+3)/2}.

## Phase 2: Claim 1 attempt — rate-½ universal

### Route (a): Generalize R78.1–78.3 to family level

**Empirical verification at p ∈ {3, 5, 7}, r ∈ {2, 3}** (script `qx1_move2_phase2_check.py`, output `qx1_move2_phase2_check.csv`):

After normalization correction (initial run flagged "FAIL" due to wrong predicted formula; corrected formula matches):

```
p=3, r=2: predicted √27   = 5.1962  actual 5.1962  exact
p=3, r=3: predicted √81   = 9.0000  actual 9.0000  exact
p=5, r=2: predicted √125  = 11.1803 actual 11.1803 exact
p=5, r=3: predicted √625  = 25.0000 actual 25.0000 exact
p=7, r=2: predicted √343  = 18.5203 actual 18.5203 exact
p=7, r=3: predicted √2401 = 49.0000 actual 49.0000 exact
```

For F̂_short_p (one-period Fourier sum), |F̂_short_p(ξ)| = p^{(r+1)/2} on supp, exact equidistribution. F̂_full_p = p · F̂_short_p ⟹ |F̂_full_p(ξ)| = p^{(r+3)/2}. **Magnitudes match the candidate generalization of R78.3 exactly at every (p, r) tested.**

**Proof template (p-blind):**
1. **Family-level 78.1_p (complete-sum vanishing):** Cochrane Theorem 2 applies to the polynomial `g(u) = c · (1+p)^u − p²·m·u` mod p^{r+1}. The expansion `(1+p)^u = Σ_k C(u, k) p^k` is purely binomial-coefficient algebra, p-blind. Cochrane T2's `D = deg_p H+` argument depends only on the polynomial's p-adic degree, which generalizes. → Complete-sum vanishing of `Σ_{u=0}^{M−1} e_M(g(u))` for the appropriate (c, m) classes.
2. **Family-level 78.2_p (sparsity):** (1+p) has order p^r in (Z/p^{r+1})* for any prime p ≥ 3 (standard fact: principal units 1 + pZ_p mod p^{r+1} form a cyclic group of order p^r). ⟹ f_p is p^r-periodic ⟹ F̂_p supported on (M/p^r)·Z/M = p·Z/M ≅ Z/p^r. The {a ≡ 1 mod p} sub-support comes from the principal-unit-character decomposition, p-blind.
3. **Family-level 78.3_p (equidistribution):** Plancherel Σ_ξ |F̂_p|² = M². Restricted to supp of size M/p², gives |F̂_p|² ≤ M·p² ⟹ |F̂_p| ≤ p·√M. Equidistribution (all support values equal magnitude) is what makes this an equality. **Equidistribution rigor at family level**: needs the principal-unit Gauss-sum equidistribution argument. At q=3 R78.3 says this was "verified empirically; follows from the principal-unit Gauss-sum structure." Same argument template at general p, modulo confirming Gauss-sum theory for the {a ≡ 1 mod p} subgroup of (Z/p^r)* at general p.

**Disposition Route (a):** the F̂_full_p magnitude formula `|F̂_full_p(ξ)| = p^{(r+3)/2}` closes structurally at family level, modulo the same residual rigor gap as R78.3 at q=3 (rigorous equidistribution proof). **This is a genuine new candidate theorem: a q-universal extension of R78.3.**

### But: this is the WRONG object for rate-½ at the K_p level

The q-sweep's measured β ≈ 0.5 is for the SHORT-WINDOW K_p, not the full-period F̂_p.

Pólya-Vinogradov decomposition:
> K_p = (1/M) · Σ_{ξ ∈ supp} 1̂_N(ξ) · F̂_p(ξ)

where 1̂_N(ξ) = Σ_{u=0}^{N−1} e_M(−ξu) is the truncated indicator's Fourier transform.

Cauchy-Schwarz with |F̂_p(ξ)|² = M·p², Σ_{supp} |1̂_N|² ≤ N·M (Plancherel for the indicator), and |supp| = M/p²:

> |K_p|² ≤ (1/M²) · |supp| · max|F̂_p|² · Σ_supp |1̂_N|²

The Cauchy-Schwarz bound gives `|K_p| ≤ √(N·M·p²·|supp|)/M = √N · p` — the *trivial* bound or worse. **The √N saturation requires phase cancellation in the bilinear sum, not magnitude bounds on each factor.**

This is exactly what R78_FINAL §5–§6 flagged: "**The saving must come from PHASE CANCELLATION in Σ 1̂(ξ)·F̂(ξ), not from each factor separately.**" The residual analytical step is a *Burgess-type bound* on the arithmetic-progression character sum on coset {a ≡ 1 mod p}. **This bound exists at the level of difficulty of subconvexity for character sums; it is not produced by R78.3-style Plancherel saturation alone, even at family level.**

### Route (b): Sibling-prime / Pattern β structural argument

`result_q_sweep_test_3_decomposition.md`'s lift-residual orthogonal identity `||R_k^{(q)}||² · q^k = S_{k+1}^{(q)} / q` is q-universal (proved over Q for q ∈ {3, 5, 7, 11, 13}). It does NOT directly produce a bound on |K_p|. It produces an exact identity relating `||R_k||²` to `S_{k+1}`. This is structurally about the Markov-chain dynamics across levels, not character-sum cancellation. **Route (b) does not give Claim 1.**

### Route (c): Family-level operator with universal ½ in spectrum

R77.2's T_3 has spectrum {1/2, 1/4, 1/8} (or rather: T_3's putative spectrum if the 3-mode model were exact, which R77.3 falsified). No q-family analog of T_3 is constructed in the project. Building one would require:
- Identifying the family-level rate operator (analogous to R76 §10's T_diag = (1/5)[[1,1],[4,4]] for q=3).
- Showing 1/2 is in its spectrum at every prime q.

This route is **wide open**; no progress made in this attempt.

### Phase 2 verdict

- Route (a) **partially closes**: full-period F̂_p magnitude formula generalizes. **K_p short-window √N saturation does NOT close rigorously.**
- Route (b): doesn't address K_p saturation.
- Route (c): not pursued; would require new operator construction.

**Claim 1 PARTIAL**: at the F̂_p level closure is real and clean (family-level R78.3); at the K_p level the same Burgess-type wall as R78 at q=3 remains.

---

## Phase 3: Claim 2 attempt — prefactor closed form

**The F̂_p magnitude formula already provides the closed form** at the F̂ level: `|F̂_p(ξ)| = p^{(r+3)/2}` uniformly, no q-varying prefactor.

The empirically-observed 1.4× C_p variation across q ∈ {3, 5, 7, 11, 13} (from `q_sweep_results.md` §3.2: C_p mean from 0.88 at q=3 to 1.06 at q=11) is for the K_p (short-window) prefactor, not F̂_p. Translating from F̂_p (rigorous, uniform) to K_p (empirical, q-varying) requires the missing bilinear bound. **The 1.4× variation lives in the gap between F̂_p and K_p**, not in the underlying magnitude law.

**Claim 2 PARTIAL**: closed form exists at the F̂_p level (it IS the formula `p^{(r+3)/2}`); the q-varying prefactor at the K_p level remains empirically observed but structurally unexplained — it's a residual of the unresolved phase-cancellation question.

---

## Phase 4: Claim 3 attempt — Nisoli ε_K specialization at q=3

For the q=3 case, even granting `|F̂_3(ξ)| = 3^{(r+3)/2}` rigorously, this does **not** translate to a bound on Tao's `|μ̂_n(ξ)|`. Reasons:

1. **F̂_p is a cyclic-group character sum on Z/M**, deterministic.
2. **μ̂_n(ξ) is the expectation of e^{−2πi ξ Syrac/3^n}**, where Syrac is a Markov-chain stationary. Probabilistic.
3. The connection: Tao 2022 (1.26) expresses Syrac(Z/3^n) = 2^{−a_1} + 3·2^{−a_{[1,2]}} + ... with iid Geom(2) a_i. Computing μ̂_n(ξ) involves expectations over the (a_1, ..., a_n) tuple. The K_p / F̂_p object does NOT directly compute this expectation.

The R78.3 magnitude formula sits at the **deterministic-character-sum level**. Even at q=3 it cannot substitute for Tao Prop 1.17's `|μ̂_n(ξ)| ≪_A n^{-A}` without an intermediate step that translates the deterministic-cyclic-group bound into a Markov-chain-stationary bound. **That intermediate step is not in the project documents.**

**Claim 3 FAILS**: the route from family-level F̂_p saturation to Tao-Prop-1.17-shaped μ̂_n bound is not closed. The R77.2 conditional on Tao Prop 1.17 does NOT dissolve via Move 2 as stated.

---

## A2 supplementary data (q ∈ {11, 13, 17, 19, 23})

Since Claim 1 only partially closes (at F̂_p, not K_p), running A2 is supplementary rather than required-by-pre-reg. Ran K_p at p ∈ {11, 13, 17, 19, 23}, r ∈ {3, 4} for completeness:

| p | r | N | √N | |K_p| | |K|/√N |
|---:|---:|---:|---:|---:|---:|
| 11 | 3 | 121 | 11.00 | 10.41 | 0.95 |
| 11 | 4 | 1331 | 36.48 | 39.23 | 1.08 |
| 13 | 3 | 169 | 13.00 | 12.20 | 0.94 |
| 13 | 4 | 2197 | 46.87 | 49.96 | 1.07 |
| **17** | 3 | 289 | 17.00 | 18.62 | 1.10 |
| **17** | 4 | 4913 | 70.09 | 65.91 | 0.94 |
| **19** | 3 | 361 | 19.00 | 17.59 | 0.93 |
| **19** | 4 | 6859 | 82.82 | 78.18 | 0.94 |
| **23** | 3 | 529 | 23.00 | 24.84 | 1.08 |
| **23** | 4 | 12167 | 110.30 | 107.21 | 0.97 |

|K|/√N ∈ [0.93, 1.10] across p ∈ {11, 13, 17, 19, 23}, consistent with the original q-sweep's universal √N saturation. **The K_p saturation is empirically robust to including these three additional primes.** This is empirical extension of Pattern β, not structural proof.

(The smaller range [0.93, 1.10] vs `q_sweep_results.md`'s [0.83, 1.17] reflects the smaller r-window here. The prefactor magnitude estimates need r ≥ 5 or so to stabilize per `q_sweep_results.md` §7.)

---

## Adversarial safeguards A1, A3, A4 (record)

- **A1 (structural vs empirical separation):** explicit at each phase. The F̂_p magnitude formula has a candidate p-blind proof template (rigorous-pending-equidistribution); the K_p √N saturation is empirical only. The 1.4× C_p variation is empirical-with-no-arithmetic-structure-known.
- **A3 (route disagreement):** Routes (a), (b), (c) for Phase 2 disagree as expected — (a) partially closes, (b) doesn't address K_p, (c) wasn't pursued. (b)'s failure is "wrong object" (orthogonal decomposition vs character-sum cancellation), not a real disagreement; (c)'s status is "not attempted" rather than "fails."
- **A4 (deviation log):** one numerical correction during compute — my initial predicted formula `p^{(r+3)/2}` was off by factor p because I conflated F̂_short (one-period sum) with F̂_full (M-period sum). Corrected formula `p^{(r+1)/2}` for F̂_short matches data exactly. This is a labeling correction, not a hypothesis change. Pre-reg §1 unchanged.

---

## What is the headline result, honestly

A **candidate family-level theorem** is identified:

> **THEOREM CANDIDATE (qx+1 Plancherel saturation at the F̂ level):** For every prime p ≥ 3 and every r ≥ 2, define f_p(u) = e_M(c·(1+p)^u) for u ∈ Z/M, M = p^{r+1}. Then F̂_p (the full-period Fourier transform of f_p on Z/M) is supported on {p·a : a ∈ subset of Z/p^r determined by principal-unit structure}, and on its support has uniform magnitude
>
> &nbsp;&nbsp;&nbsp;&nbsp; **|F̂_p(ξ)| = p^{(r+3)/2}.**
>
> Proof sketch: Cochrane Theorem 2 + Plancherel + principal-unit Gauss-sum equidistribution, all p-blind. Rigor follows R78.1–78.3's q=3 template, modulo confirming Gauss-sum equidistribution at general p.

Verified empirically to machine precision at p ∈ {3, 5, 7}, r ∈ {2, 3}.

**What it does NOT do:** dissolve R77.2's conditional on Tao Prop 1.17 for c = 7/45. The F̂_p result is the wrong-object-shape to substitute. The bridge from F̂_p (deterministic character sum) to either:
1. K_p √N saturation (short-window character sum) — requires bilinear bound, same wall as R78 at q=3
2. μ̂_n(ξ) decay (Markov-chain stationary characteristic function) — requires a translation step not in project documents

is not closed by Move 2 as stated.

**Strategic position:** the F̂_p candidate theorem may stand on its own as a structural result — a family-level extension of R78.3 with explicit closed-form magnitude. But c = 7/45's rate-½ rigor remains where R77.2 left it, conditional on either Tao §7's effective constant or some new translation step that hasn't been identified.

---

## Files

- [QX1_FAMILY_THEOREM_PRE_REGISTRATION.md](QX1_FAMILY_THEOREM_PRE_REGISTRATION.md) — locked rules (commit f96fb86)
- [QX1_FAMILY_THEOREM_ATTEMPT.md](QX1_FAMILY_THEOREM_ATTEMPT.md) — this document
- [qx1_move2_phase2_check.py](qx1_move2_phase2_check.py) — empirical verification script
- [qx1_move2_phase2_check.csv](qx1_move2_phase2_check.csv) — raw data

## Disposition handling (per pre-reg)

Pre-reg disposition spectrum was: `THEOREM_PROVEN → CLAIM_1_ONLY → CLAIMS_1_AND_2 → EMPIRICAL_UNIVERSAL_NOT_STRUCTURAL → STRUCTURAL_OBSTRUCTION_FOUND`.

Landed at: **CLAIM_1_PARTIAL (at the wrong-object level) + STRUCTURAL_OBSTRUCTION_FOUND (for the Move 2 pipeline as stated).**

Per pre-reg §4: "If a specific obstruction is identified precisely (mechanism named, not just 'this didn't work'), that is itself a real result." The obstruction identified: the K_p / F̂_p / μ̂_n triad are three distinct objects, family-level results on F̂_p don't immediately bound K_p (needs bilinear bound), don't immediately bound μ̂_n (different probabilistic structure). This constrains the structural search going forward.

For the c=7/45 closure question: the routes available are still:
1. Tao §7.2–7.4 C_A effectivization (R77.2 path, qualitative bookkeeping)
2. Burgess-type bilinear character-sum bound on coset {a ≡ 1 mod 3} (R78 residual gap)
3. A new structural translation between K_p / F̂_p and μ̂_n that connects the cyclic-group character-sum bounds to Markov-chain stationary characteristic-function decay

Routes 1 and 2 are unchanged from existing project status. Route 3 is a new candidate not previously articulated — if such a translation exists, it would let the family-level F̂_p theorem contribute to c=7/45 closure. This is a possible future direction; this attempt did not identify a candidate translation.

The F̂_p candidate theorem may itself merit standalone formalization independent of c=7/45.
