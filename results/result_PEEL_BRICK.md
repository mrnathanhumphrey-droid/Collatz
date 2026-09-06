# RESULT — PEEL-ONE-BRICK: is S_∞'s transcendence reducible to log₂3? **OUTCOME C (NULL) — below floor.** (2026-09-06)

**Probe:** `probes/probe_peel_brick.py`. Tests the **dimension** of the ∞-unipotent kernel (is it secretly
1-dimensional, all transcendence on the single generator `log₂3`?), NOT the value of S_∞ (MAHLER forbids).
Kill-first, pre-registered A (one brick) / B (∞-dim) / C (null). Guardrails: de-rotate the period-9.06
(`ω=ln2`) modulation BEFORE any algebraicity test (the trap — `log₂3` is baked into the rotation);
high-precision PSLQ stable 50→100; mandatory negative control; report the floor.

## Inputs (banked, not recomputed)
- **Exact `Λ_r`, r=1..7** (arbitrary precision): `Λ_1=−2/21`, `Λ_2=−1490/203889`, `Λ_3≈0.001319864`, …,
  `Λ_7≈0.00021484523` (from exact `S_k`, k=1..8). Sign pattern `−,−,+,+,+,−,+`.
- **Float tower `Λ_16..20`** (from `T_15..20`, ~12 sig digits after differencing): `2.34256e-4, 2.06651e-4,
  1.79568e-4, 1.56700e-4, 1.35538e-4`.
- `ω_lp = ln2 = 0.6931…` ⟹ modulation **period = 9.065 levels**; `log₂3 = 1.58496…`, `log₂5 = 2.32193…`.

## Step 0 — de-rotation (mandatory) is UNDERDETERMINED
The exact contiguous block r=1..7 spans **0.66 periods** of the 9.06-level modulation. Fitting
`log|Λ_r| = c₀ + c₁r + A·cos(ωr) + B·sin(ωr)` (4 params) on 7 points leaves **3 dof over < 1 period** —
the periodic (rotation) component is **not separable** from the geometric trend (periodic amplitude fit
`1.78` is meaningless at <1 period; residual is fit-noise `≤2.2e−1`, not signal). **The de-rotation the
guardrail demands cannot be performed on the only high-precision (exact) data available.** [de-rotation floor]

## Step 1 — PSLQ at achievable precision: the negative control fits identically
Float ratios `ρ_r (r=16..) = 0.882, 0.869, 0.873, 0.865`; Aitken limit proxy `L ≈ 0.8701` — **but only
~8 digits** (float-tower differencing), and **not de-rotated** (Step 0 floor). Running PSLQ on
`{1, L, g, L·g, g², L²}`:

| tol | `g = log₂3` | `g = log₂5` (neg control) |
|---|---|---|
| 8 digits (data floor) | `[0,−13,−2,24,−5,−8]` (spurious) | `[345,343,−221,−230,55,50]` (spurious) |
| 20 digits | huge-coeff junk | huge-coeff junk |
| 50 digits (guardrail) | **None** | **None** |

**The negative control `log₂5` is treated identically to `log₂3`** — any "relation" at the achievable ~8
digits is a numerical coincidence (the R9-D failure mode), and none survives to the 50-digit stability the
guardrail requires (the data has < 50 digits). Guardrail **UNMET**.

## Step 2 — lobe-ledger cross-check inherits the same floor
The restratification `Λ_r = L_r − M_r` (result_P6J) is exact only where `Λ_r` is exact (r≤7), so the lobe
cross-check spans the same 0.66 periods — no independent lift.

## Floor and verdict
- **Needed** for the A-vs-B distinction: de-rotation over **≥ ~2.5 periods (≥ ~23 contiguous layers)** AND a
  de-rotated residual limit to **≥ 50 digits**, stable 50→100.
- **Achievable:** 0.66 periods exact (**7 layers** — more is blocked: exact `S_k` needs the
  doubly-exponentially-denominatored stationary chain, infeasible past k≈8); ~8 digits (float tower).
- **Deficit: ~1.8 periods of length AND ~42 digits of precision.**

**⟹ OUTCOME C (NULL).** Whether the kernel is 1-dimensional (transcendence reducible to `log₂3`, Outcome A)
or genuinely ∞-dimensional (Outcome B) is **undecidable with the banked data** — the distinction sits far
below the achievable-k/precision floor. Not A, not B. The trap was dodged (a raw ratio test *would* have
"found" `log₂3`; the negative control shows any such fit is an artifact at this precision).

**Capstone unchanged:** this probe can neither sharpen it (A) nor is it needed to defend it (B). No S_∞
value estimated or claimed. Not at stake: MAHLER, the value floor `2·T_20=0.473177`, R1–R30.
