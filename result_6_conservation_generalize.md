# Result 6 (qx+1 paper) — does R76 generalize to `(Z/q^k)*`? Conservation YES (free), leading-mode NO (q=3 is structurally unique)

**Date:** 2026-07-15. **Verdicts: H_CONS CONFIRMED (verbatim port). H_EQUAL REFUTED (decisively). H_LEAD_BREAKS CONFIRMED.**

**Headline: the qx+1 paper's stated "one line to a full theorem" is mischaracterized — the identity it names generalizes for free, and generalizing it does not close the gap.**

Probe: `probe_6_conservation_generalize.py`. Log: `result_6_conservation_log.txt`. Scopes the open step flagged in `result_5_universal_rate.md` §45 / `QX1_UNIVERSAL_RATE_WRITEUP_2026_07_14.md` §32. Compute: seconds (N ≤ 343).

## What R5 claimed was left

R5 (H_PROVED-at-mechanism) flagged exactly one open step:

> "…that the *sub-leading* character contributions to `‖π_k‖²` do not perturb the leading `(1/3)^k` rate — i.e. that the diagonal self-overlap dominates uniformly in q. **This is exactly the R76-style conservation identity generalized to `(Z/q^k)*`**; … the clean algebraic proof of uniform domination is the remaining line."

This probe takes that at its word and tries to port R76. **The framing is wrong twice over.**

## Validation gate — PASS

Reproduces `result_76_conservation_law.md` §3 exactly before any new claim:

| quantity | got | R76 target |
|---|---|---|
| `M_2(1) = S_2` | 0.4761904762 | 0.4761904762 (=10/21) ✓ |
| `M_2(1+3) = M_2(4)` | −0.2380952381 | −0.2380952381 ✓ |
| `S_2 = −2·M_2(4)` | 0.4761904762 = 0.4761904762 | ✓ |

## H_CONS — CONFIRMED. Conservation ports verbatim, and it is FREE.

`Σ_{j=0}^{q-1} M_{k+1}(η_0 + j·q^k) = 0` for every `η_0 ∈ (Z/q^k)*`:

| q | k→k+1 | N | max\|Σ_j M\| |
|---|---|---|---|
| 3 | 1→2 | 9 | 5.6e−17 |
| 3 | 2→3 | 27 | 2.8e−17 |
| 5 | 1→2 | 25 | 1.2e−16 |
| 5 | 2→3 | 125 | 5.2e−16 |
| 7 | 1→2 | 49 | 1.6e−15 |
| 7 | 2→3 | 343 | 8.9e−16 |
| 11 | 1→2 | 121 | 2.9e−15 |
| **9** (composite) | 1→2 | 81 | 7.8e−16 |
| **15** (composite) | 1→2 | 225 | 2.2e−15 |

**And it is provable by inspection — no new work needed.** R76's Thm 76.1 proof is a complete-character-sum vanishing argument. The inner sum is

> `Σ_{j=0}^{q-1} e^{2πi·rξj/q}` = `q` if `q | rξ`, else `0`.

The chain is supported on `gcd(r,q)=1`, so `q | rξ ⟺ q | ξ`; restricting the ξ-sum to `q ∤ ξ` kills every term. **Nothing in that argument uses q=3** — not primality (composites q=9,15 confirm), not the modulus, only `gcd(r,q)=1`. Theorem 76.1 generalizes to every odd q with the proof unchanged.

**Consequence: conservation cannot be the missing step.** If it were, the paper would already be finished — the port costs one line and it is done above.

## Lift structure — the pairing generalizes; the COUNT does not

Thm 76.3 (`S_{n+1} = −2·M_{n+1}(1+3^n)`) rests on **Lemma 76.2**: among the lifts of `η_0`, one is self-inverse and the rest are mutual inverses. Measured at `η_0 = 1`:

| q | lifts of 1 | self-inverse | mutual-inverse pairs `(q−1)/2` | conservation yields |
|---|---|---|---|---|
| **3** | 3 | {1} | **1** | 1 eqn, **1 unknown → SOLVED** |
| 5 | 5 | {1} | 2 | 1 eqn, 2 unknowns |
| 7 | 7 | {1} | 3 | 1 eqn, 3 unknowns |
| 11 | 11 | {1} | 5 | 1 eqn, 5 unknowns |

The **pairing itself ports fine** — `(1+j·q^k)(1+(q−j)·q^k) = 1 + q·q^k + j(q−j)q^{2k} ≡ 1 mod q^{k+1}`, so `j ↔ q−j` are mutual inverses and `M(η)=M(η^{−1})` makes M equal on each pair. This is **directly visible as a palindrome** in the measured values (below). `−1` is never a lift of 1 (for `q^k > 2`), so exactly one lift (namely 1) is self-inverse for every odd q.

**What breaks is arithmetic, not structure:** Thm 76.3 collapses `S_{n+1}` onto a *single* mode only because `(q−1)/2 = 1`, and

> **`(q−1)/2 = 1  ⟺  q = 3`.**

## H_EQUAL — REFUTED (the branch that would have saved the port)

Pre-registered hopeful branch: if `M_{k+1}(1+j·q^k)` were independent of `j`, conservation would give `S_{k+1} = −(q−1)·M_{k+1}(1+q^k)`, reproducing 76.3's `−2` as `−(q−1)` and porting cleanly. **Dead.** Threshold was 1e−10; measured spread 0.97–5.71.

| q | k | `M(1+j·q^k)`, j=1..q−1 | spread | S | `−(q−1)·M(1+q^k)` | match |
|---|---|---|---|---|---|---|
| 3 | 1 | −0.23809524, −0.23809524 | 5.6e−17 | +0.47619048 | +0.47619048 | **YES** |
| 3 | 2 | −0.23078734, −0.23078734 | 2.8e−17 | +0.46157468 | +0.46157468 | **YES** |
| 5 | 1 | +0.14118367, −0.82480986, −0.82480986, +0.14118367 | 9.66e−01 | +1.36725238 | −0.56473468 | NO |
| 5 | 2 | −1.05191607, −0.08149309, −0.08149309, −1.05191607 | 9.70e−01 | +2.26681831 | +4.20766426 | NO |
| 7 | 1 | −1.08486180, +0.05446389, −1.09761910, ×2 (palindrome) | 1.15e+00 | +4.25603402 | +6.50917078 | NO |
| 7 | 2 | −0.60494839, −2.31455690, −2.05134538, ×2 (palindrome) | 1.71e+00 | +9.94170134 | +3.62969035 | NO |
| 11 | 1 | −0.629, −2.855, +2.837, −1.374, −2.878, ×2 (palindrome) | 5.71e+00 | +9.79842301 | +6.29463950 | NO |

Note the palindromic structure at every q≥5 — that is Lemma 76.2's pairing surviving intact. The identity fails not because the pairing breaks but because there are **too many pairs to pin with one equation.**

## Disposition — what this does and does not touch

**Does NOT touch R5's result.** `S_k^(q) ~ (q/3)^k` remains H_PROVED-at-mechanism: the derivation (`1/3 = Σ4^{−v} = E[2^{−v}]`, q-blind) never routed through 76.3, and the adversarial-q falsifier is unaffected. Pillars 2 (`c̃_q=(q−3)/q`) and 3 (`δ_q≈0.82/ord_q(2)`) untouched.

**DOES change the roadmap.** The paper's stated remaining step is not one line and is not that line:
1. The named identity (conservation) generalizes **for free** — proof unchanged, confirmed to 1e−15 including composite q.
2. Generalizing it is **insufficient** — the q=3 leading-mode collapse it feeds has no analogue at q≥5.
⇒ **Uniform domination needs a genuinely different argument** — a direct decay/spectral bound on sub-leading characters, not a conservation collapse. The R76 route is BLOCKED, not incomplete.

**Consolation prize (real, and paper-usable):** this is a **second, structural sense in which q=3 is critical** — it is the unique odd q for which conservation determines the leading mode. The paper currently justifies "q=3 is the critical case" only via `S_k^(3)→7/15`. This adds a phase-boundary reason with a one-line proof (`(q−1)/2=1 ⟺ q=3`), which is exactly the register the Bernoulli-convolution template (Erdős/Solomyak/Hochman) wants.

**Not at stake:** THEOREM_C_745 (c=7/45), Th 78.1–78.3, R81b, ε_k. This probe touches only the standalone qx+1 paper's rigor step.

## ⚠️ FLAGGED — object-identification gap (next probe)

`M(1)` measured here at q≥5 grows (1.37, 4.26, 9.80 at q=5,7,11). **R76's `S_n` (coprime-restricted Plancherel of the deviation `d_k`) and the q-sweep's `S_k^(q) = q^k‖π_k‖²` (stationary L² mass) share a name and coincide at q=3 — but this probe did NOT verify they are the same object off q=3.** Every cross-q claim built on R76 machinery is contingent on that identification. **Pin it before building further.**

_Reporting discipline: the falsifier (H_EQUAL) was the branch that would have saved the port, was pre-registered to be tested hardest, and lost. The gate reproduced R76's published table to 10 digits before any new claim was made. The result is reported as a blocked route + a reframing, not dressed up as progress; the one genuine byproduct (q=3 structural uniqueness) is stated as a consolation prize, not a headline. Author's structural priors: now 1-for-6 this arc (H_CONS/H_LEAD_BREAKS called correctly; H_EQUAL's hopeful branch lost, as pre-stated)._
