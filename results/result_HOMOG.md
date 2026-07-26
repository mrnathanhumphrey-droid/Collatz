# RESULT — HOMOG: the ln|π̂| shape is NOT forced (skew drifts 20σ); homogeneous across nodes; tail diffuse (2026-07-26)

**Probe:** `probe_homog.py`. Measurement only, π̂ (additive) throughout — no ρ/channel/U² quantity enters. Tests whether
the stationary log-spectrum shape is FORCED by ultrametric homogeneity + R8 self-similarity (⟹ exactly stationary,
theorem) or merely EMERGENT (settling with decaying drift) or NEITHER (drift persists, no settling). L = ln|π̂_k(a)|
over `{1≤a<3^k : 3∤a}` (= units mod 3^k, count 2·3^{k−1}), k=8..16.

## A — GATE passes, with visible drift
| k | count | μ | var | skew | exkurt |
|---|-------|---|-----|------|--------|
| 8 | 4374 | −4.93 | 0.436 | −0.719 | 1.48 |
| 12 | 354294 | −7.25 | 0.540 | −0.596 | 1.47 |
| 16 | 28697814 | −9.55 | 0.621 | −0.468 | 1.21 |

Banked −0.65 / +1.4 reproduced (skew −0.72→−0.47, exkurt 1.48→1.21 across k=8..16). **But skew drifts monotonically
toward 0 and exkurt down** — visible, not noise.

## B — DECIDER: NOT forced; and NEITHER (drift resolved, settling not)
Standardized cumulants vs k, deparitied (period-2 guard):
| cumulant | const (scatter) | lin slope | **dep \|slope/SE\|** | geom c∞ / ρ / R² |
|----------|-----------------|-----------|------------------|------------------|
| **κ₃ skew** | −0.599 (0.085) | +0.0325 | **20.67** | +1.01 / 0.980 / 0.979 |
| κ₄ exkurt | +1.399 (0.131) | −0.0307 | 2.04 | −0.09 / 0.980 / 0.354 |
| κ₅ | −2.75 (0.59) | +0.028 | 0.02 | — / — / 0.013 |
| κ₆ | +8.88 (3.49) | +0.448 | 1.48 | — / 0.690 / 0.180 |

- **FORCED is dead:** κ₃ drifts at **20.67σ** (deparitied — survives the period-2 guard, so a real drift not a wobble;
  slope +0.032/k). The shape is **not exactly stationary.** (Pre-reg forced needed |slope/SE|<2 for *every* cumulant.)
- **EMERGENT not supported:** the geometric fit pins **ρ=0.980** (near 1), and the skew increments are ~constant
  (+0.03/k, no deceleration through k=16) — geometric does not beat linear. **ρ<1 is not resolved.**
- Per the pre-registered rubric (drift resolved + ρ not resolved below 1) ⟹ **NEITHER.** The "stationary shape" reading
  is wrong: the shape is still moving at k=16 and where it settles (if it does) is not visible in k≤16. (κ₄ variance
  increments DO decelerate — 0.033→0.018 — so *variance* mildly settles, but the standardized *shape* via κ₃ does not.)

## C — HOMOGENEITY holds: the shape is class-independent
At k=12, standardized shape per 3-adic class:
- **by a mod 9:** skew −0.590..−0.600, exkurt 1.43..1.49 — **essentially identical across all 6 unit residues.**
- **by v₂(a):** skew −0.596..−0.662, exkurt 1.07..1.68 — tight; **high-v₂ does NOT carry a different shape** (contra the
  pre-registered expectation that the sup-region v₂ would stand out).
- **by a mod 27:** mean skew −0.611, spread 0.127 over 18 classes — modest.

So the shape **is homogeneous** (same at every node within a level) — consistent with ultrametric homogeneity — **but
that homogeneous shape drifts with k (the depth).** A homogeneous *non-stationary* flow down the tree: within-level
node-independence + across-level drift are not in tension (k is the depth).

## D — the heavy left tail is DIFFUSE, not structured
Bottom 1% of |π̂(a)| at k=12 and k=16 (both identical):
- **v₂(a):** 50 / 29.7 / 10.4 / 5.4 / 2.2 / 1.2 % — **≈ the uniform baseline** (50/25/12.5/…). Not enriched in high-v₂.
- **a mod 9:** 16–17% each — uniform.
- **a/N:** near-0 10%, mid(.4–.6) 19%, near-1 10%, median 0.499 — **spread across the whole range**, not concentrated
  near 3^k/2 or on the ⟨2⟩-orbit (⟨2⟩ = all units anyway, 2 primitive root — that partition is degenerate).

**The suppressed frequencies are a diffuse set, not a structured one — no object/handle there** (contra the pre-reg
hope that the tail would localize).

## Verdict
- **The spectral shape is NOT forced/exactly stationary** — κ₃ drifts at 20σ. "Depth is provably spent via a forced
  stationary shape" is FALSE.
- **NEITHER emergent nor forced** by the rubric — drift is resolved but settling (ρ<1) is not, over k≤16. The shape is
  a homogeneous but still-moving flow; its limit is not visible.
- **Homogeneity itself holds** (class-independent shape) — the one structural positive, consistent with the ultrametric
  tree.
- **The left tail is diffuse** — no structured suppressed set to grab.

Net: HOMOG's best cases (convert "depth spent" into a theorem; localize the tail) both **fail** — the shape drifts and
the tail is diffuse. Honest negatives; nothing moves S_∞ (as pre-stated). What's real: a 20σ skew drift with unresolved
settling, and confirmed within-level homogeneity. Not at stake: LOGNORMAL, VALPROFILE (denom theorem), NORMCHECK,
CHANNEL_ID, R1–R30. Cost 614s (k≤16).
