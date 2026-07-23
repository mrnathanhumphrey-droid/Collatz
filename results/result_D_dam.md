# Probe D — the dam (transient deflation) — **cross-validation gate FAILS: the pure bore region spans only ~½ period, so at k=12 the bore's period is still INFERRED, not measured. Rate ~0.98 reconfirmed; the k=12 hope corrected.**

**Date:** 2026-07-22. Probe `probes/probe_D_dam.py`. Deflate `B_r = Λ_r − T_r`, `T_r = A·(½)^r`, to expose the bore
(ρ≈0.984, period ~9) across r=3..11. Exact rational `Λ_r` for r≤7 (ε exact k≤8); 15-digit float r=8..11
(`result_epsilon_11.csv` + `S_12`). **The pre-registered D-B cross-validation gate fails — reported per spec.**

## D-A — characterize the transient (gate): no clean region, μ deviates from ½ by the bore
- `A_r = Λ_r·2^r`: **A_3 = 0.010559, A_4 = 0.010404, A_5 = 0.010461** — spread **1.49%**, agree to ~2 digits, **not
  the pre-registered ≥4**.
- Best single rate `μ = (Λ_5/Λ_3)^{1/2} = 0.497687`, i.e. `μ − ½ = −2.3×10⁻³` (the ratios `Λ_4/Λ_3 = 0.4927`,
  `Λ_5/Λ_4 = 0.5028` straddle ½). This exceeds the pre-registered `1×10⁻⁴` stop-threshold.

**Read:** by the size ledger the bore is already **12–47%** of the transient at r=3,4,5, so there is **no
fully-clean region**; the μ-deviation *measures* that contamination, it does **not** refute the exact-½ transient
(R26 derived ½ subcritically). LSQ transient (μ=½ fixed) on r=3,4,5: `A* = 0.010525`.

## D-C — the deflated ladder (artifact; bore part under-resolved)
| r | Λ_r | T_r | B_r=Λ−T | \|B\|/0.984^r | sign B |
|---|---|---|---|---|---|
| 3E | 1.320e−3 | 1.316e−3 | +4.27e−6 | ~0 | + |
| 4E | 6.50e−4 | 6.58e−4 | −7.54e−6 | ~0 | − |
| 5E | 3.27e−4 | 3.29e−4 | −1.98e−6 | ~0 | − |
| 6E | −3.39e−4 | 1.64e−4 | **−5.03e−4** | 0.00055 | − |
| 7E | 2.15e−4 | 8.22e−5 | +1.33e−4 | 0.00015 | + |
| 8f | 3.69e−4 | 4.11e−5 | +3.28e−4 | 0.00037 | + |
| 9f | 3.64e−4 | 2.06e−5 | +3.44e−4 | 0.00040 | + |
| 10f | 3.91e−4 | 1.03e−5 | +3.80e−4 | 0.00045 | + |
| 11f | 3.86e−4 | 5.14e−6 | +3.81e−4 | 0.00046 | + |

`B_3,4,5 ≈ 0` by construction (the LSQ fit zeroed them — circular in the fit region). The bore is only exposed at
r=6..11: a `−` trough at r=6, then `+` rising r=7..11 — **locally monotone, not a resolved oscillation.**

## D-B — CROSS-VALIDATION GATE: **FAIL** (the decisive item)
Fit the bore recurrence `B_{r+1}=a B_r + b B_{r−1}` on the pure region r=8..11: `a=−0.850, b=+2.051` ⟹ **real roots,
`ρ=1.92`** — the four points are locally monotone (`3.28, 3.44, 3.80, 3.81 ×10⁻⁴`, ~0.44 periods of 9), so **no
damped oscillation can be fit to them.** Backward-extrapolation to r=3,4,5 then gives ~2.3–2.7×10⁻⁴ against the
deflated `B_3,4,5 ≈ 4×10⁻⁶` — rel err **5000–13000%**. Gate **fails**.

**Diagnosis (per pre-reg: "third mode OR transient not exactly geometric — report which"):** neither, in fact — it
is **data insufficiency**. The pure region spans `<½ period` and is locally monotone, so the 4-parameter oscillation
is **under-determined**, and the early terms are circular (fit-zeroed). This is not a proven third mode; it is that
**k=12 does not yet expose ≥1 clean bore period.** Per the guardrail, D-C..F are void as a clean measurement.

## D-E — rate/period (windows disagree ⟹ period unmeasurable)
| window | steps (periods) | ρ | period |
|---|---|---|---|
| r=5..9 | 4 (~0.44) | 0.652 | 4.46 (θ=80.7°) |
| r=7..11 | 4 (~0.44) | 1.052 | real roots — none resolved |
| r=6..11 holdout | 5 (~0.55) | **0.981** | real roots — none resolved |

The windows **disagree** (ρ = 0.65 / 1.05 / 0.98); no period is pinnable. The holdout `ρ = 0.981` is at least
consistent with the 0.984 rate (F1/R81), but "real roots" means the oscillation is **not resolved in-window.** Free
check: amplitude cap `|Λ_r| ≤ 7/45 = 0.1556` holds at **every** r ✓.

## D-D — sign sequence on deflated data
`sign(B_r)` r=3..11: `+ − − − + + + + +` (but B_3,4,5 signs are fit-residual noise ~0). One genuine change (r=6→7).
Period 6 excluded; 9/12/18 untestable (need ≥10/13/19 deflated terms; have 9). Under-resolved.

## Status
**Probe D: the dam does not hold at k=12 — the cross-validation gate fails on data insufficiency, and the bore's
period remains inferred, not measured.** **A** — no clean transient region (A spread 1.5%, μ off ½ by 2.3e−3 = the
bore's 12–47% contamination; ½ not refuted). **B GATE FAIL** — the pure region r=8..11 is ~0.44 periods and locally
monotone, so the bore fit is under-determined (spurious real-root ρ=1.92) and the deflated early terms are circular;
the two-mode single-bore model is **not validated** — but the cause is resolution, not a proven third mode. **C** —
the deflated ladder is the artifact, but its bore part (r=6..11, 6 terms) spans <1 period. **E** — windows disagree
(ρ=0.65/1.05/0.98), period unpinnable; holdout ρ=0.981 consistent with 0.984 (rate OK). **D/F** — void
(under-resolved / oscillation not resolved).

**Consequence for the crux (owed to the pen).** The optimistic read — "k=12 ≈ 1.1 periods ⟹ the bore can finally be
measured" — is **corrected**: the *usable* clean region (transient ≲5%) is only r=8..11 ≈ **½ period** and is locally
monotone, and deflation is circular in the transient-fit region, so the period is still **not measurable**. What is
reconfirmed: the **rate ~0.98** (holdout ρ=0.981, consistent with R81's 0.984) and the amplitude cap `|Λ|≤7/45`. To
actually *measure* the bore's period (and run the D-F two-vs-three-mode / beat test) needs the pure region to span
≥1 period — roughly **ε to k≈17–20** — which requires exact (or high-precision) `S_k` at those levels, presently out
of reach. Honest negative, pre-registered: the gate was built to catch exactly this, and it did. No fitting beyond
the labeled windows; the D-B failure and window-disagreement reported plainly, not smoothed into a false period.
