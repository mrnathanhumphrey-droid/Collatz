# Result — PROBE R: R1 hit the STOP condition (closure NOT a clean union; modular (e,γ) gauge broken). R2 (L=4) deferred (local wall + reduction basis now uncertain).

**Date:** 2026-07-16. Forward work for the invariant hunt. **R1 STOPPED as instructed** — the (e,γ) class closure is not a clean union of classes under the modular-carry gauge. **R2 not fired** (heavy-compute-greenlight rule; and its reduction basis is what R1 just invalidated). Probe `probes/probe_phase2b_R1.py`, log `logs/probe_phase2b_R1_log.txt`, tables `outputs/class_table_L{2,3}.tsv`.

## R1 — STOP: closure is NOT a clean union of classes
Convention (as required to state): **2^e = b·a⁻¹ mod 3^L**, e∈Z/D; γ ∈ Z/3^L, **untransformed** (the modular diagonal gauge (a,b,γ)→(sa,sb,γ)). Move type Δe = δ_a − δ_b mod D; e' = e + Δe.

| L | classes | forward closure of {e=0, γ≠0} | expected (Nathan) | clean union? | partial (gauge-split) |
|---|---|---|---|---|---|
| 2 | 54 | **29** closure / 25 support | 21 / 33 | **NO** | **16** |
| 3 | 486 | **291** closure / 195 support | 171 / 315 | **NO** | **240** |

**Robustly verified — not a bug:**
- State-level forward closure of the diag-carry seed reaches **exactly 126 states** at L=2 — **identical to the LALB BFS from C** (adjacency confirmed correct).
- Those 126 states span **29 classes: 13 full + 16 partial.** A partial class carries some but not all of its gauge orbit. Example (e=3, γ=4): states (2,7,4), (4,5,4), (8,1,4) are reached, but their gauge-partners (1,8,4), (5,4,4) are **not**.
- ⇒ **reachability is not invariant under the diagonal shift** — the (e,γ) projection with modular γ is not a faithful quotient of the dynamics.

**Diagnosis (why, and what it means):**
- This is consistent with **L-B**: only the k=0 co-invariant eigenvector is gauge-invariant; the *full operator is not gauge-equivariant*. (No contradiction with L-B — that is a property of one eigenvector; this is a property of state reachability. Both hold.)
- Relabeling e (generator direction) cannot fix it — the split is in the **carry** action, not the e-labeling.
- **This is exactly the constraint Nathan extracted from the S/J refutations:** *any invariance must act on the carry as an INTEGER map, not a modular one.* The modular-γ gauge breaks; the clean 21/33 requires the integer-carry gauge.

**Action (per instruction (b)): STOPPED and reported loudly. No fishing for the right twist** (the decoration anti-pattern Nathan named). The handoff tables are emitted with a **3-way flag (F=full-closure / P=partial-gauge-split / 0=support)** so the raw (e,γ) digraph structure is preserved for whichever gauge is correct — but the invariant hunt as posed (over 21 clean classes) is **paused pending the exact integer-carry gauge/coordinate.**

**Correction to prior banking:** `result_phase2b_LALB.md` inferred "126 = 21 clean classes" loosely (never verified cleanness). R1 corrects this: the C-closure is **29 classes, 16 partial — not clean.** **L-A (no return to Δ) is UNAFFECTED** (that is about Δ, which stays out of the closure) and **L-B is UNAFFECTED** (the k=0 eigenvector gauge-invariance is a separate, still-true fact). Only the "21 clean classes" reduced-chain premise is retracted under the modular gauge.

## R2 — G (L=4 partner): DEFERRED (local wall + reduction basis invalidated)
**Not fired.** Two independent blockers:
1. **Local wall (sizing).** The full L=4 operator is `D²·q^L = 54·54·81 = 236,196` states with ~`2.3e8` nonzeros. A direct LU factorization of a matrix this size and connectivity fills in catastrophically — infeasible on the local box. Per the instrument law (direct/LU only near the EP) *and* the heavy-compute-→-Lambda rule, I do **not** fire it here. This is the standing wall.
2. **Reduction basis now uncertain.** The "localization-reduced" route would lean on the (e,γ) class reduction — which R1 just showed is **not clean** under the modular gauge. Until the correct (integer-carry) gauge is fixed, the reduced basis for a tractable L=4 solve is not defined.

**Consequences, held open (NO extrapolation, per guard):**
- **Rate-law point 4:** the gap sequence `2.9e-3, 1.0e-4, ?` is **unadjudicated at L=4.** The zero-weight 27^{−L} pattern *predicts* ~3.8e-6 — reported as an **untested prediction, not a result.**
- **Braid point 3:** partner-above-c₀ (L=2) → below (L=3) → **L=4 side unknown.** No pre-registered pick; **not adjudicated.**
- c₀ reference at L=4 = `1/3 + (2/3)·2^{−54}` stands (closed form), but the *partner* is what walls.

**Path forward for R2 (when unblocked):** either (a) resolve R1's gauge → a clean reduced chain → localization basis → small direct solve, or (b) a Lambda-sized direct/shift-invert solve on the full 236k operator (needs greenlight + sizing). No local heavy fire without go.

## Deliverables
- `outputs/class_table_L2.tsv`, `outputs/class_table_L3.tsv` — raw (e,γ) digraph, F/P/0 flagged. **Handoff is emitted but the closure is non-clean; the invariant hunt opens only after the gauge is corrected.**
- This file + STATE entry.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, the H-gates, L-A, L-B (as corrected above), the circulant-family completeness, the J-refutation. No `r_q` value changes.

_Reporting discipline: R1 hit the pre-specified STOP and is reported as a STOP, not massaged toward 21. The non-cleanness is verified robustly (126 states matches LALB; explicit gauge-split example) and diagnosed (modular vs integer carry — Nathan's own constraint), not hand-waved. I did NOT fish for a gauge that yields 21 (Nathan's stated anti-pattern). The prior "21 clean" inference is retracted; L-A/L-B stand. R2 is reported as a wall + a blocked-reduction, with the 27^{−L}=3.8e-6 figure held as an untested prediction and the braid side left unadjudicated — no extrapolation._
