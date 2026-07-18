# Perturbation probe — c=7/45 structural stability under (3+ε) multiplier

**Date:** 2026-05-06.
**Brief:** test whether S_k(ε) — the perturbed-multiplier analog of S_k under
(3+ε)x+1 dynamics with methodology (a) rounded perturbation — converges to 7/15
(stable), shifts smoothly with ε (smooth landscape), or breaks structurally.

## Verdict — Outcome D: methodology (a) fails

> Rounded perturbation **cannot answer the structural-stability question** for
> c = 7/45 at the tested ε range and k=5. Two regimes exist, both unhelpful:
>
> - **|ε| ≤ 0.002:** perturbation is below the rounding threshold for every
>   coprime state at k=5 (max coprime r = 242, threshold |ε| ≥ 0.5/242 ≈ 0.00207).
>   K(ε) = K(0) exactly → S_k(ε) = S_k(0) exactly → centered derivative dS/dε = 0
>   trivially. No signal.
>
> - **|ε| ≥ 0.005:** S_k(ε) jumps chaotically. Values across the tested grid
>   span [−18.62, +2.53] versus the unperturbed S_5 = 0.466. Centered numerical
>   derivatives at h = 0.005, 0.01, 0.02, 0.05 are −87, +170, −506, −25 —
>   inconsistent in both magnitude and sign. The structural-stability question
>   is masked by methodology artifacts.
>
> No "smooth derivative dS/dε at ε=0" can be extracted from this methodology.
> Per the brief's walk-back gate ("If S_k(ε) is non-monotone or chaotic in ε..."),
> stopping at k=5 and reporting rather than extending to k=6, 7.

## Sanity check (passed)

At ε = 0, S_5 = X_5 − X_4 = 3.5341611514 − 3.0686462316 = **0.4655149198**.
Canonical from `c_seven_forty_fifth.md` table: 0.46551492. Diff −2.5×10⁻¹⁰
(matches at float64 precision). The probe's K-builder reproduces the
unperturbed chain correctly.

## Setup recap

For each tested k, build:
```
base_real = (3 + ε) * r + 1
base_int  = round_half_up(base_real)            # round-to-nearest, ties → up
base_mod  = base_int mod 3^k
For each v in {1..M=ord_{3^k}(2)}:
    target = (base_mod * inv2^v) mod 3^k
    if target ≢ 0 mod 3:
        K[r, target] += 2^(-v) / (1 − 2^(-M))
```
Leak handling: rows where `base_int ≡ 0 mod 3` get a self-loop (every target
would land on a non-coprime state); rows with partial per-v leaks are
renormalized to sum to 1.

S_k(ε) := X_k(ε) − X_{k−1}(ε)  where X_k(ε) := 3^k · ‖π_k(ε)‖².

ε grid: {±0.05, ±0.02, ±0.01, ±0.005, ±0.002, ±0.001, 0}.
Built at k=5 and k=4 separately (S_5 needs both levels).

## Results — full table

| ε | X_4(ε) | X_5(ε) | S_5(ε) | ε_5(perturb) | n_pert k=5 | n_pert k=4 | leaks k=5 | leaks k=4 |
|---|---|---|---|---|---|---|---|---|
| −0.0500 | 4.6380 | 5.2740 | +0.6360 | +0.169 | 155/162 | 47/54 | 52 | 20 |
| −0.0200 | 2.5781 | 4.1879 | +1.6098 | +1.143 | 145 | 37 | 66 | 33 |
| −0.0100 | 4.9056 | 4.0395 | −0.8661 | −1.333 | 128 | 20 | 66 | 20 |
| −0.0050 | 3.0686 | 2.8185 | −0.2502 | −0.717 | 95 | 0 | 95 | 0 |
| −0.0020 | 3.0686 | 3.5342 | **+0.4655** | −0.00115 | **0** | **0** | **0** | **0** |
| −0.0010 | 3.0686 | 3.5342 | **+0.4655** | −0.00115 | **0** | **0** | **0** | **0** |
| **0.0000** | **3.0686** | **3.5342** | **+0.4655** | **−0.00115** | **0** | **0** | **0** | **0** |
| +0.0010 | 3.0686 | 3.5342 | **+0.4655** | −0.00115 | **0** | **0** | **0** | **0** |
| +0.0020 | 3.0686 | 3.5342 | **+0.4655** | −0.00115 | **0** | **0** | **0** | **0** |
| +0.0050 | 3.0686 | 1.9450 | −1.1237 | −1.590 | 96 | 0 | 0 | 0 |
| +0.0100 | 1.9820 | 4.5146 | +2.5327 | +2.066 | 129 | 21 | 62 | 0 |
| +0.0200 | **24.7003** | 6.0828 | −18.6175 | −19.084 | 146 | 38 | 45 | 4 |
| +0.0500 | 7.1335 | 5.2949 | −1.8386 | −2.305 | 156 | 48 | 52 | 13 |

(Headers: n_pert = states with `base_int ≠ 3r+1`; leaks = states with full row leak.
 ε_5(perturb) = S_5(ε) − 7/15.)

The four ε values within ±0.002 give bit-identical results to ε=0 — no rounding
shift triggers anywhere. The remaining 9 values produce S_5 values scattered
across two orders of magnitude with sign flips, no monotone structure.

## Why methodology (a) fails — structural diagnosis

Three interacting issues:

### 1. The rounding floor

Rounding shift for state r is `δ(r, ε) = round((3+ε)r+1) − (3r+1)`. δ is 0
unless `|εr| ≥ 0.5`. At k=5 (max r = 242):

- **|ε| < 1/(2·242) ≈ 0.00207** → δ ≡ 0 for every state → K(ε) = K(0) exactly.

So any |ε| below ~0.002 is invisible to the discrete chain. The probe's brief
ε set {±0.001, ±0.002} sits squarely in this dead zone, which is why centered
derivatives at h ≤ 0.002 are exactly zero.

### 2. Discontinuous transition fans at threshold crossings

When `δ` shifts from 0 to ±1, the state's full v-weighted transition fan
gets *replaced* by a different fan — the new `base_mod` lands on a different
residue class, and `inv2^v · base_mod mod 3^k` for v ∈ {1..M} traces a
completely different orbit through coprime states. There's no continuity in
the transition kernel: a 1-state, 1-shift change can move ~M weight units to
arbitrary other states.

The base-shift histograms at increasing |ε| show the staircase pattern:

| ε | shift histogram |
|---|---|
| ±0.001..±0.002 | {0:162} (all states unperturbed) |
| ±0.005 | {0:67, ±1:95} (~half states shifted by 1) |
| ±0.01 | {0:34, ±1:66, ±2:62} (most states shifted by 1 or 2) |
| ±0.02 | {0:17, shifts ±1..±5: 12-34 states each} |
| ±0.05 | shifts 0..±12 distributed across all 162 states |

Each shift step transforms the row's transition pattern wholesale — there's
no per-row gradient with respect to ε, so the assembled K(ε) is a step
function of ε with O(n_states) discrete plateaus.

### 3. Full-row leaks at "wrong" rounding directions

For some states r and shifts δ, `(3r+1+δ) mod 3` lands at 0 instead of 1,
giving full-row leaks (every target ≡ 0 mod 3, all transitions go to non-coprime
states). The methodology converts these to self-loops, which act as absorbing
mass in π. At ε = −0.005, 95 of 162 states (59%) full-leak; π(ε=−0.005)
concentrates strongly on the self-loops.

This effect is *asymmetric in ε sign*: at ε=+0.005, all shifts are +1, and
(3r+1+1) mod 3 = 2 ≠ 0 (no leaks). At ε=−0.005, all shifts are −1, and
(3r+1−1) mod 3 = 0 (every shifted state leaks). This explains why
S_5(−0.005) = −0.250 and S_5(+0.005) = −1.124 are not symmetric — they
correspond to qualitatively different perturbed chains, not just sign-flipped
versions of one perturbation.

### 4. Cross-level leak asymmetry

S_k = X_k − X_{k−1} subtracts levels. At ε = +0.0200, the level-4 chain has
*4 full leaks* (out of 54 states), pushing X_4(+0.02) to 24.70 (versus 3.07
unperturbed). The corresponding level-5 chain has *45 full leaks* but
spreads them across 162 states, giving X_5 = 6.08. The subtraction
S_5 = 6.08 − 24.70 = −18.62 isn't tracking "leading-mode plateau" anymore;
it's tracking arithmetic difference of two leak-dominated quantities.

## Pre-registered outcomes — what got returned

| Outcome | Status |
|---|---|
| **A: smooth landscape** | NO — S_k(ε) is not smooth at any tested h |
| **B: rigid stability** | NO — S_k(ε) ≠ S_k(0) for |ε| ≥ 0.005 |
| **C: discontinuous** | PARTIALLY — S_k(ε) is chaotic above the rounding threshold, but the "chaos" is methodological (rounding artifacts), not a phase transition in the underlying dynamics |
| **D: methodology fails** | **PRIMARY** |

The "stepwise constant below threshold" + "chaotic above threshold" pattern
matches the brief's walk-back gate "If the rounding artifact is large enough
to dominate the ε signal at small ε..." → halt, report, methodology (b) needed.

## What we *can* conclude despite outcome D

1. **The c = 7/45 result is invisible to small-ε rounded perturbation.**
   Below the rounding threshold, K(ε) ≡ K(0), so the probe degenerates to
   the unperturbed chain. This is a methodological non-result, not evidence
   of rigidity (outcome B).
2. **Methodology (a) has structural flaws that make it unsuitable for
   perturbation analysis at any k.** Higher k *worsens* the issue: at k=6,
   threshold is 0.5/728 ≈ 0.000687 (tighter dead zone) and leak count grows.
3. **The structural-stability question remains open.** A probe answering it
   needs methodology (b) (continuous-state framework on Ẑ_3 or T = R/Z) or
   methodology (a') (probabilistic rounding to smooth out the step
   function). Both are real builds, deferred.

## Suggested follow-ups

Ranked by tractability:

1. **Methodology (a'): probabilistic rounding.** For each state, weight the
   floor and ceil targets by `(1 − frac)` and `frac` where `frac = (3+ε)r + 1
   − floor((3+ε)r + 1)`. This gives a stochastic transition that depends
   smoothly on ε. K(ε) then has well-defined ε-derivative and dS/dε is
   numerical-derivative-friendly. Cost: same as current probe, runs in seconds
   at k=5. Would directly resolve A vs B vs C.
2. **Methodology (b): continuous-state on T = R/Z.** Define the perturbed
   chain on the circle (or Z_3) with smooth transition kernel. Build by
   discretizing the kernel at fine resolution (n_bins ≫ 3^k) and finding
   stationary. More principled but requires architectural rebuild. Cost:
   day-scale.
3. **Methodology (c): keep multiplier integer, perturb v-distribution.**
   Replace the truncated-Geom(1/2) with truncated-Geom(1/2 + ε). Stays in
   the integer-arithmetic Markov chain framework but tests sensitivity to
   the v-distribution assumption. Different question — closer to "how
   sensitive is 7/15 to the Geom(1/2) ansatz" than "how sensitive is it
   to the multiplier 3" — but easier to execute and answers a closely
   related stability question.
4. **Run methodology (a) at much larger ε** (e.g., ε ∈ {0.5, 1, 2}). At
   ε=2, multiplier becomes 5, which is the q=5 chain — and we already have
   q=5 results from the q-sweep. So this just reproduces the qx+1 family
   (which is what the q-spectrum probe already covered, with `c̃_5 ≈ 0.487`
   and `c̃_q = (q−3)/q` at q ≥ 11). Not a new probe.

The cleanest cheap follow-up is **methodology (a')** — it preserves the
discrete chain structure, makes K depend smoothly on ε, and runs in seconds.
Recommended next pass if the user wants to actually answer the
structural-stability question.

## Files

- [perturbation_probe.py](perturbation_probe.py) — script
- [perturbation_probe.log](perturbation_probe.log) — full stdout
- [result_perturbation_S_k_curve.csv](result_perturbation_S_k_curve.csv) — S_k(ε) per ε
- [result_perturbation_derivative.csv](result_perturbation_derivative.csv) — numerical dS/dε estimates
- [perturbation_findings.md](perturbation_findings.md) — this writeup

## Honest framing

The brief asked a real structural question (is c=7/45 perturbation-stable?).
The probe's methodology (a) was the natural starting point — discrete chain
preserved, ε perturbation injected via rounding. But the methodology has a
fundamental impedance mismatch with the question: the discrete chain only
"sees" perturbations that cross integer rounding thresholds, so any ε small
enough to be a "small perturbation" is invisible, and any ε large enough to
register triggers wholesale row replacements rather than infinitesimal nudges.

This is the kind of finding that's useful to document because someone (perhaps
future-self) might otherwise repeat the methodology choice. Outcome D is a
methodology verdict, not a structural one — c = 7/45's stability question
stays open, requiring methodology (a') or (b) to actually answer.
