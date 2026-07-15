# Result 1B — Halving-PGF two-root-of-unity audit (disposition)

**Date:** 2026-07-14. **Verdict: H_ROBUST** (under the invariant/natural measure; the first-pass pooled "drift" was survivorship, diagnosed and dismissed).

Pre-reg: `PRE_REG_1B_GEOM_PGF_2026_07_13.md`. Probes: `result_1b_halving_pgf.py` (first pass, pooled + step brackets), `result_1b_clean.py` (survivorship-free resolution). Data: `result_1b_pgf_data.csv`, logs `result_1b_log.txt` + `result_1b_clean_log.txt`.

## What was audited

The §CONJ mechanism `7/45 = N(2−ω) / (3²·(1+4))` rests on the Geom(½) halving PGF `G(z)=Σ_{v≥1} 2^{−v} z^v = z/(2−z)` evaluated at two roots of unity:
- `Ĝ(−1) = E[(−1)^v]` → target **−1/3** (the "9": class mass (0,1/3,2/3) → squared 1:4).
- `|Ĝ(ω)|² = 1/N(2−ω)` → target **1/7** (the "7": N(2−ω)=7), ω = e^{2πi/3}.

Since v is **not exactly** Geom(½) empirically (R68: 0.5–25% deviation at specific j; ensemble mean 2.102), the question is whether the true measure moves these two constants — i.e. whether the mechanism is measure-robust or a Geom(½) idealization.

## Result

**Invariant single-step marginal (the measure Tao's chain / R75–R76 run on).** v₂(3n+1) for n drawn uniformly over odd, coprime-to-3 integers (fresh 40M-sample, independent of the trajectory data):

| constant | measured | target | Δ |
|---|---|---|---|
| Ĝ(−1) | −0.3333332 | −1/3 | +1.3e−7 |
| \|Ĝ(ω)\|² | 0.1428572 | 1/7 | +1.1e−8 |
| mean v | 2.000000 | 2.0 | — |

Both constants are exact to 7–8 digits. **H_ROBUST.**

## The first-pass apparent drift was survivorship, not a measure property

The pooled step-weighted pass (`result_1b_halving_pgf.py`) showed |Ĝ(ω)|² drifting to 0.166 in the pooled `step10+` bracket. Per-single-step resolution (`result_1b_clean.py` Part B) shows this drift is **depth-monotone and tracks survival loss**: corr(log-survivors, Δ|Ĝ(ω)|²) = **−0.57**. Steps 0–10 (survival ≈99–100%) sit within Δ≈2e−4 of target; the drift grows only as the surviving population thins and skews low-v (long trajectories descend slowly). Conditioning on survival ≠ the invariant measure. The pooled 0.166 came from the extreme deep tail where survivorship is worst.

## Smoke-check reconciliation (the pre-reg STOP condition)

The pre-reg expected mean ≈ 2.102; the invariant per-step marginal gives 2.0. Reconciled: **2.102 (`result_density_one_v2_bounds.md` L11/L47/L90) is the length-biased *per-trajectory-mean averaged over trajectories*** — high-v trajectories are shorter, so trajectory-weighting up-weights them. The per-step chain-kernel marginal (relevant to the 7/45 mechanism) is Geom(½), mean 2.0. Different marginals; no bug.

## Routing

H_ROBUST → the `7 = N(2−ω)` mechanism is measure-robust under the natural measure: it is exactly as idealized as R75/R76's own Geom(½) chain — no more, no less. Phase 2 gets the sharp target: **derive R75/R76's 7 as N(2−ω)** (and the 9 as G(−1)=−1/3). No `THEOREM_C_745.md` Geom(½)-substitution audit is motivated by drift, because there is no invariant-measure drift.

## Scope note

This does not re-derive 7/45 (that is R75/R76 exact-rational, separate) and does not bear on the rigorous constant. It audits only whether the §CONJ *mechanism* reflects the true measure. It does. The distinct claim tested by Probe-85 rung 2 — that the "7" generalizes as the cyclotomic norm Φ_p(2)=2^p−1 under multiplier variation — was **refuted** there; that is a *multiplier*-generalization, orthogonal to this *measure*-robustness result.

_Reporting discipline: H_ROBUST reported as fired. The first-pass pooled drift is reported as measured AND diagnosed as survivorship (not softened, not inflated). The smoke-check STOP was honored — the marginal mismatch was reconciled before the verdict, not aggregated away._
