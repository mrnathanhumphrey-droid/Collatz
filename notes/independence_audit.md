# Independence Audit — Validation Task 2

**Date:** 2026-05-03. Audits 7 cross-result independence claims in the
framework synthesis. Verdicts: **(α)** independent confirmation,
**(β)** consistent characterization (shared underlying object), **(γ)** circular.

**Headline:** mostly **(β)**. Of 7 claims, 0 are independent confirmations,
6 are consistent characterizations, 1 is propose-not-demonstrate.

| # | Claim | Verdict |
|---|---|---|
| 1 | M_closed independent of Chang | (β) — same dynamics, different transfer operators |
| 2 | Forward-orbit and inverse-tree are independent data | (β) — same dynamics, different sampling ensembles |
| 3 | R58 and R60 independently reproduce D_avg | **(β) empirically: cross-Pearson +0.92** |
| 4 | Bridge equation Constants 1+2+4 mutually confirm | (β) — internal consistency check, not independent |
| 5 | Bridge equation predicts Tao's K_h | **(γ-input) — K_h is INPUT not OUTPUT** |
| 6 | R22 σ-quartile machinery would close R58's residual | proposed-not-demonstrated |
| 7 | Multiple methods locate same H-dim | **(β) — no algebraic identity; ~10% agreement window** |

## Headline reframings

These are the textual changes the audit recommends for any external-facing
synthesis (paper, Lagarias/Tao/Chang correspondence, UMD application):

- "Three independent confirmations" → "three methodologically distinct
  characterizations of the same underlying object (the trajectory measure
  on Z₂)"
- "Result 23's λ_max gives 0.6755 ≈ Chang's 0.68" → "Result 23's branching
  dim 0.338 doubles to 0.675; Chang's exact 0.694; our spatial info dim
  D₁ = 0.608. Five distinct values in a 0.07-wide window; no algebraic
  identity holds. Consistency at the ~10% level, not sharp identification."
- "Bridge equation matches Tao's K_h" → "Bridge equation uses Tao's K_h as
  structural input; the slope-1.000 verification confirms internal
  consistency at K_h, not independent derivation of K_h"
- "R22 machinery closes R58's residual" → "R22 σ-quartile Esscher-tilt
  closure is a candidate not yet demonstrated"

---

## Claim 1: M_closed and Chang's P operate independently

### Setup

- **Chang's kernel P** (`chang_kernel_pi.py`): forward Syracuse map
  m → (3m+1)/2^v on residues mod 64 with depth-13 lift averaging.
  P[r → r'] = Pr[next residue r' | current residue r] under
  Lebesgue-uniform integers in [0, 2^13) above each cylinder.
- **R23 M_closed** (`experiments/80_dim_h_validation.py`): inverse Syracuse
  rules at modulus 2^k. Doubling rule m → 2m always; inverse-3 rule
  ((m-1)·3⁻¹) when m ≡ 4 mod 6. Pure algebraic on residues.

### Trace

Both built from the same map (Collatz/Syracuse), but as DIFFERENT operators:
- Chang's P: forward transfer operator, mod 64, with cylinder-13 averaging
- M_closed: backward transfer operator (inverse map), mod 2^k, no averaging

Neither uses output from the other. The branching-dim 0.338 = log(λ_max)/log 2
from M_closed and Chang's exact dim_H = log(φ)/log 2 = 0.694 are
genuinely distinct invariants.

### Verdict

**(β)** Methodologically distinct, but rooted in the same Collatz dynamics.
NOT "independent confirmation" — both ultimately encode different facets
of the same dynamics. Their numerical proximity (0.338×2 ≈ 0.694) is a
near-miss within 0.02, **already walked back per dim_h_validation outcome
(γ)**.

---

## Claim 2: Forward-orbit and inverse-tree data are independent

### Setup

- **Forward-orbit data**: random m₀ ∈ [1, 2^N], walk until m=1. Generates
  D_avg, sigma cache, all R50/R51/R53 results, R60 kernel.
- **Inverse-tree data**: build tree from m=1 backward via doubling and
  inverse-3 rules, value-truncated at max_value. Generates R52/R58 weights,
  R23 BFS density, R59 zadic mass, R61 spatial Z_q.

### Trace

Both probe the same Collatz reality. Every forward orbit reaching m=1 traces
a path that IS in the inverse tree (if max_value is large enough). However,
the SAMPLING distribution differs:

- Forward at N=2^32: orbit-starts uniform in [1, 2^32]; visits
  intermediate m ∈ [1, peak] weighted by orbit time-average
- Inverse-tree at max_value=2^22: nodes m ≤ 2^22 weighted by
  inverse-tree multiplicity (variant a = subtree size, etc.)

These produce DIFFERENT measures on different ensembles, but probe the same
underlying dynamics.

R58 explicitly noted this: "value-truncated inverse tree from m=1 (truncated
to m ≤ N = 2^22) — matches the integer-uniform sampling regime" of the
forward orbits. By design, the value-truncation aligns the two samplings.

### Verdict

**(β)** Same dynamics, different sampling ensembles. The two data sources
are NOT independent — they're two views of the Collatz forward map.
The fact that both reproduce the same target (D_avg) is consistent
characterization, not independent confirmation.

---

## Claim 3: R58 and R60 independently reproduce D_avg

### Empirical check (PRIMARY for this audit)

```
Pearson(R58 variant_a, R60 D_pred)        = +0.9243
Pearson(R58 variant_a, R58 empirical D)   = +0.8665
Pearson(R60 D_pred,    R60 empirical D)   = +0.8033
Pearson(R58 variant_a, R60 empirical)     = +0.8677
Pearson(R58 empirical, R60 empirical)     = +0.9184

MAE(R58 - D_avg)  = 0.205
MAE(R60 - D_avg)  = 0.212
MAE(R58 - R60)    = 0.082    <-- much smaller than each MAE-vs-target
```

### Interpretation

R58 and R60 produce predictions that are **92% Pearson-correlated** with
each other and only 0.08 MAE apart — much closer to each other than either
is to the empirical target (~0.21). They are **the same identification
computed via two methodologies**, not independent confirmations.

Even the empirical TARGETS differ slightly (R58 used D at t=90 from 2M
orbits at N=2^32; R60 used late-t average D_avg from 10M orbits at smaller
N), and those targets cross-Pearson at +0.918. So even the "ground truth"
D_emp comparators are correlated, not independent observations.

### Verdict

**(β) DECISIVE.** R58 and R60 are highly redundant predictions of the same
underlying object. The "two different methods both work" framing should be
"the same identification recovered by two methodologies, with internal
consistency at Pearson 0.92."

The substantive content remains: the trajectory measure has a real
identification that multiple methods can capture. But "independent
confirmation" is wrong; "robust under methodological variation" is right.

---

## Claim 4: Bridge equation Constants 1+2+4 mutually confirm

### Setup (per `tao_bridge_findings.md`)

`s_mean(r) ≈ α_det(r) + K_h · log(N/f(N)) + ε`

- **Constant 1: ⟨α_det⟩** = E[prefix_steps] − K_h · ⟨descent during prefix⟩
  — **algebraic invariant of prefix algebra at modulus 2^k**, computed
  exactly to +6.23 across all k. NO empirical fit.
- **Constant 2: K_h = 3/log(4/3)** — taken from Tao 2022 (input).
- **Constant 4: ε** structural correction — measured EMPIRICALLY from
  σ data at multiple N; ≈ −2.45 with finite-N drift ~0.01.

### Trace

Constant 1 is purely algebraic (no data). Constant 2 is cited from external
theory. Constant 4 is empirical.

The "internal consistency" check is: assuming K_h, does
σ̄ − K_h · log N − ⟨α_det⟩ = ε converge to a stable constant across N?
Empirically yes (−2.4468 at 2²⁵ → −2.4574 at 2³², drift 0.01 << K_h · log 2).

This is a bridge-equation **internal consistency check**, not three
mutually confirming derivations.

### Verdict

**(β)** Internal consistency between assumed K_h, algebraic ⟨α_det⟩, and
measured ε. Not "independent" because Constants 1 and 2 are not derived
from data — they're algebra and citation. Constant 4 is the only empirical
piece, and the test is whether it's stable, not whether it agrees with
something else.

---

## Claim 5: Bridge equation predicts Tao's K_h

### Trace

`tao_bridge_findings.md` line 5: "K_h = 3/log(4/3), slope at K_h = 1.000
± 0.005 across 40 verification cells".

K_h is the SLOPE coefficient assumed in the linear regression. The
verification fits `s_mean(r) ≈ K · log(N/f(N)) + const` and reports
the slope K. K = 1.000 ± 0.005 means: **assuming K_h = 3/log(4/3) is the
correct unit, the regression slope is 1.000**. Equivalently:
the regression slope of σ̄ vs log N is 3/log(4/3) ± 0.5%.

So the bridge framework MEASURES the slope and finds it ≈ Tao's K_h. That
IS empirical confirmation that the empirical σ̄ scales with rate K_h.

### Refinement

K_h = 3/log(4/3) is Tao's PARAMETERIZATION of "fictional mean Collatz
descent rate per step." The regression doesn't derive this rate from
first principles inside the bridge framework — it assumes σ scales
linearly with log N (and verifies that empirically) and reports the
fitted slope. Comparing the fitted slope to 3/log(4/3) is a check that
matches Tao's prediction.

So the framing is mixed: K_h's NUMERICAL VALUE is empirical (slope of σ̄
vs log N regression). K_h's CLOSED-FORM (3/log(4/3)) is from Tao's external
heuristic. The bridge framework **verifies** rather than derives K_h.

### Verdict

**(γ-input)** K_h is taken as input/parameterization, not derived from
first principles inside the bridge framework. The "matches Tao's K_h" is
a verification that empirical slope = 3/log(4/3) numerically, which is
real cross-validation of Tao's heuristic but not a new derivation.

Reframe: "Bridge equation empirically verifies that σ̄ scales with slope
3/log(4/3) per Tao's heuristic prediction" — clean and accurate.

---

## Claim 6: R22 σ-quartile machinery would close R58's residual

### Trace

R58 (`inverse_tree_weighting.md` line 122, 142, 154, 171):
- "30% residual at r=5, r=23, r=13 reflects Esscher tilt under survivor
  conditioning beyond what time-averaged visit count captures (Result 22's
  σ-quartile tilt machinery)"
- "Open: variant (a) ↔ D_emp gap closure via Esscher tilt"
- "Concrete next moves: 1. Esscher tilt closure ... Test whether Pearson
  lifts from +0.86 to >+0.95"

### Verdict

**Proposed not demonstrated.** R58 IDENTIFIES R22's machinery as the
candidate closure mechanism but has not RUN the closure test. The
Pearson +0.95 claim is conjectural; no Esscher-tilt-closed variant has
been computed.

Reframe: "R58 leaves a 30% residual at QSD extremes that may close via
R22's σ-quartile Esscher-tilt machinery (untested)."

---

## Claim 7: Multiple methods locate same H-dim

### Empirical check

```
Quantity                                Value     Source
log(lam_max)/log(2)                     0.337726  R23 (Furstenberg branching)
2 * log(lam_max)/log(2)                 0.675452  R23 derived
Chang exact = log(phi)/log(2)           0.694242  Chang dim_H of survivor set
R61 spatial info dim D_1                0.607772  trajectory measure
Chang heuristic                         0.6800    Chang's reported approximation
R59 dim_q2(k=12)                        0.67      Chang-aligned q=2 mass dim
R59 dim_q2(k=7)                         0.83      coarse-scale q=2
R59 dim_q2(k=15)                        0.54      fine-scale q=2

Pairwise differences:
D_1 vs 2*log(lam_max)/log(2):           0.0677    (NO algebraic identity)
D_1 vs log(phi)/log(2):                 0.0865    (distinct)
2*log(lam_max) vs log(phi):             0.0188    (smallest, still distinct)
0.68 vs D_1:                            0.0722    (within ~10%)
```

### Interpretation

There are **five distinct numbers** all sitting in the [0.6, 0.7] window:
- 0.608 (info dim D₁ of trajectory measure)
- 0.675 (R23's 2 × branching dim)
- 0.680 (Chang heuristic, rounded)
- 0.694 (Chang exact log(φ)/log(2))
- 0.67 (R59's q=2 mass dim at k=12, scale-dependent)

**No algebraic identity holds among them.** The pairwise differences range
0.019 to 0.087. The "agreement" is at ~10% relative level, far from sharp.

R59's dim_q2(k) sweep 0.83 → 0.67 → 0.54 across k ∈ {7, 12, 15} is the
standard multifractal-crossover signature of a finite-resolution
correlation-dimension estimate. Catching the value 0.67 at k=12 is a
SCALE-DEPENDENT artifact, not a sharp invariant.

The dim_h_validation walk-back already addressed the M_closed ↔ Chang
proximity. R61 sharpens this by placing all the candidates on the f(α)
spectrum: they're all near the q=1 (information dim) region but no two
land on exactly the same point.

### Verdict

**(β)** Multiple methods produce dim-like quantities that all live in the
[0.6, 0.7] region but are demonstrably distinct at 4-decimal precision.
"Multiple methods agree on H-dim" should be reframed as "several methods
produce information-dimension-region values consistent within ~10%, all
locating in the q ∈ [0.5, 1] portion of f(α)."

---

## Cross-cutting observation

ALL the results in this body of work derive from one underlying object:
**Collatz dynamics on Z**. Every method (forward orbits, inverse tree,
M_closed transfer, Chang's averaged kernel, size-stratified Markov, MF
Z_q, bridge equation σ scaling) is a different reduction or summarization
of that same dynamics.

This is FINE for the substantive science: characterizing the same object
multiple ways is informative robustness. But it is NOT "independent
confirmation" in the methodological sense, because there's only one
underlying reality being probed.

The honest framing for the body of work: **methodologically distinct
characterizations of the trajectory measure on Z₂, all consistent with
each other within their respective error bars, with the size-stratified
Markov framework (R60) and spatial multifractal Z_q (R61) providing the
strongest formal identifications**.

## Reframing recommendations for v3.6 / external correspondence

For each cross-result claim in the framework synthesis chapter, replace:

1. "Result 23 and Chang independently confirm dim ≈ 0.68" →
   "Result 23's branching dim doubled is 0.675; Chang's exact is 0.694;
   these are within 0.02 but not algebraically identified. Both lie in
   the q≈1 (information dim) region of the trajectory measure's f(α)
   spectrum."

2. "Three measures on shared Λ" →
   "Three methodologically distinct projections of the same Collatz
   forward map: Chang's cylinder-averaged kernel, Result 23's M_closed
   inverse operator, and the empirical D_avg. Their consistency is
   consistency of the same dynamics, not independence."

3. "Result 60 and Result 58 both reproduce D_avg, providing independent
   confirmation" →
   "Result 60 (size-stratified Markov on forward orbits) and Result 58
   (value-truncated inverse-tree subtree-size) both predict D_avg; their
   predictions are themselves Pearson-correlated at +0.92 — they recover
   the same identification by methodologically distinct routes, not
   independent confirmations."

4. "Bridge equation matches Tao's K_h" →
   "Bridge equation's empirical regression of σ̄ on log N gives slope
   1.000 ± 0.005 in units of K_h = 3/log(4/3), verifying Tao's heuristic
   prediction at this empirical level. K_h is structural input; bridge
   equation does not derive it."

5. "R22 machinery closes R58's residual" →
   "R22's σ-quartile Esscher-tilt machinery is the proposed closure for
   R58's 30% residual at QSD extremes; the closure test has not been
   executed."

These reframings preserve substantive content while honestly representing
evidentiary structure. They make the body of work MORE defensible to
serious external readers, not less interesting.

## Files

- `independence_audit.md` — this document
- `independence_audit_compute.py` — empirical Pearson and dim-relation checks
- `independence_audit_compute.txt` — output log
- `result_dependency_graph.md` — directed graph of result dependencies
- `claim_strength_summary.csv` — per-claim verdict summary
