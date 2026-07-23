# Probe W2 — the tailor on the second moment — **does NOT derive ½; EXPOSES ½ as a transient. The η=1 triple is rank-1; the asymptotic second-moment rate is the period-9 ~0.984, not ½ — reconciling Thread-3 with the corpus's May refutation.**

**Date:** 2026-07-22. Probe `probes/probe_W2_tailor.py` (exact-rational ε/Λ ladder from the k=8 exact table +
Thm 76.3; first-moment recheck via W13). W2 asked: does the 76.2 self-inverse/pair structure, on the
second-moment channel, produce `|λ₂| = ½`? **It does not — and the way it fails is the finding.**

## W2-A — the η=1 lift-triple is **rank-1** (degenerate), exact
`triple_n = (M_n(1), M_n(1+3^{n−1}), M_n(1+2·3^{n−1}))`. By **Thm 76.3** (level m=n−1) `S_n = −2 M_n(1+3^{n−1})`,
and by **Lemma 76.0** (M-real) the two pair entries are equal. So exactly, for n=2..8:
> `triple_n = S_n · (1, −½, −½)` — **RANK-1.**

A 3×3 on a rank-1 family is degenerate: the single nontrivial eigenvalue is the S-ratio `S_{n+1}/S_n → 1`
(leading), and `M(1+3^n)` carries the **same** rate as `S_{n+1}` (they are proportional by −½). **There is no
independent ½ subdominant in the η=1 triple** — it tracks only the magnitude `S_n`, not the deviation mode.
[Wilson's pre-registered guardrail / W2-E outcome-3: the triple is too coarse.]

## W2-B — the 76.2 `−2` is exact, but it is a **within-level magnitude**, not a rate factor
`M(1) = −2·M(1+3^{n−1})` holds EXACT n=2..8 (it *is* Thm 76.3). But it relates `S_n` to its pair **at fixed n**;
every entry of the triple has the **same** rate (all `= const·S_n`). So "the −2 connects the core rate to ½" is
**false as posed** — the −2 is a magnitude relation at one level, and the rate lives in how `S_n` itself evolves,
which the triple does not resolve.

## W2-C / W2-E — the exact rate ladder: **½ is a transient (r ≤ 5); it breaks at r = 6–7**
`ε_k = S_k − 7/15`, `Λ_r = (ε_{r+1}−ε_r)/2` (= OffDiag/2), exact through the k=8 table. Ratios, **no fit**:

| k | \|ε_k\|·2^k | | r | Λ_{r+1}/Λ_r | \|Λ_r\|·2^r |
|---|---|---|---|---|---|
| 3 | 0.0407 | | 3 | −0.181 | 0.0106 |
| 4 | 0.0392 | | 4 | **0.4927** | 0.0104 |
| 5 | 0.0369 | | 5 | **0.5028** | 0.0105 |
| 6 | 0.0319 | | 6 | **−1.036** | 0.0217 |
| **7** | **0.1504** (4.7× jump) | | 7 | −0.634 | 0.0275 |
| **8** | **0.1909** | | | | |

- `|ε_k|·2^k` is **flat ~0.037 through k=6, then jumps 4.7×** to 0.150 (k=7), 0.191 (k=8). `|Λ_r|·2^r` **grows**
  (0.010 → 0.022 → 0.028). The Λ ratios are ½ **only at r=4,5** (0.4927, 0.5028 — the two R29-D read), then
  scatter and flip sign (−1.036, −0.634).
- This is **exactly the `|ε_7|·2^7 = 0.150` jump by which the corpus refuted rate-½ (Conjecture 77.2) in May 2026**
  and established the asymptotic `ρ ≈ 0.984` (from `|ε_k|^{1/k}` at k=13..16). So `ε_k` decays **slower** than
  `2^{−k}` past k=6 ⟹ the asymptotic rate is **> ½**, consistent with 0.984.

**Reconciliation (no contradiction with R26).** ½ is the mode that continues **subcritically** — R26's
`|λ₂|/ρ = 2λ²`, derived at λ = ½+ε where the period-9 oscillation is **absent** (it is q=3-critical only). **At
criticality a slower period-9 mode (~0.984) appears and dominates**, so the asymptotic approach `S_n → 7/15` is
governed by 0.984, and **½ is a faster transient sub-mode**. Two different modes; `0.984 > ½` wins asymptotically.
This resolves the long-standing ½-vs-0.984 tension and the R18-A three-numbers puzzle: they were never the same
object, and the **dominant** subdominant (true `|λ₂|`) is the period-9 ~0.984, not ½.

## W2-D — period-9 not cleanly reproduced in the first-moment fixed-level operator
The shifted operator's slow mode is `|λ₂| = 0.970 @ 60°` (period **6**, stable r=3..6); the core-units subdominant
gives arg → 0 with period 27, 81, 243 = 3^r-scaling (fixed-level artifacts). So the fixed-level first-moment model
does **not** cleanly carry period-9 — its exact home remains open (the true operator is growing-spaces; the proxy
gives period-6/3^r artifacts, not 9).

## Status
**W2: the tailor does NOT deliver ½; it exposes ½ as a transient and relocates the asymptotic rate to the period-9
~0.984.** **A** — the η=1 lift-triple is rank-1 (`S_n·(1,−½,−½)` exact), degenerate, no independent ½ subdominant.
**B** — the 76.2 `−2` is exact but a within-level magnitude, not a rate factor. **C/E** — the exact ladder shows ½
holds only r≤5 and breaks at r=6–7 (`|ε|·2^k` jumps 4.7×, `|Λ|·2^r` grows, ratios flip sign), matching the
corpus's May refutation of rate-½ and its asymptotic `ρ≈0.984`. **D** — period-9's exact home is not the
fixed-level first-moment operator (period-6/3^r artifacts).

**Consequence for the crux (the reframe, owed to the pen).** Thread-3's central `|λ₂| = ½` (R18-A, R26, R27, R29-D)
is **real but mislabeled**: ½ is the **subcritical-continuation / fast-transient** mode, valid where R26 derived it
(λ≠½), but it is **not the asymptotic rate** of `S_n → 7/15`. At criticality the slower **period-9 (~0.984)** mode —
critical-only, which the subcritical extrapolation structurally cannot see — is the **dominant** subdominant and
sets the true convergence rate. So the "derive `|λ₂|=½`" program (R25–R30 gate-hunting, the W2 tailor) was aiming
at the **faster** mode; the analytic step the theorem actually needs is control of the **period-9 0.984 object**,
which is where the corpus (Conj 77.2 refutation, T_lead spectrum {43/45, 0}) had already localized it. No fitting;
exact rank-1 triple, exact `−2`, exact ladder ratios; the ½-transient exposure and the 0.984 relocation reported
plainly as the (large) honest negative, reconciling Thread-3 with the R74–R77 corpus rather than contradicting it.
The subcritical ½ (R26) and M-reality (R80) stand; what is retracted is **½ as the asymptotic `|λ₂|`**.
