# Result — PROBE C: the (θ = e mod 3, γ) compression is REFUTED; the partner lives in the FULL gauge-invariant (e_ρ, γ)-chain (approximately). Coupling is class-dependent and non-lumpable.

**Date:** 2026-07-16. Compression recon for the hand-derivation (no proof, no rate fit). Direct methods (lumped Galerkin via sparse mat-vec). Probe `probes/probe_phase2b_C.py`, log `logs/probe_phase2b_C_log.txt`, dumps `outputs/compressed_q{3,7}_L*.tsv`.

**Headline: the ℤ/3 premise fails — the partner does NOT compress to (θ = e mod 3, γ) (rel error 0.5–0.7). It requires the FULL gauge-invariant (e_ρ, γ) subspace, and even there the compression is only APPROXIMATE: q=3 rel 1.7e-2 (L=2, near-EP twist content) → 3.2e-3 (L=3); q=7 (gapped) 3.1e-3. It is NOT an exact reduction (C3: non-lumpable, the 0.258 obstruction is present on the partner classes too). The coupling is class-dependent (C2). The derivation target is the (e_ρ, γ)-chain, not a compact θ-tower.**

## ⚠️ Correction to Probe P (found while building C)
Probe P's gauge test (`gauge_resid`) had a **k² bug** (`ω^{k·e_a}` was coded as `ω^{k²·e_a}`), which made k∈{0,6,12} spuriously coincide (all → trivial twist). **P2's "coupling concentrated in k≡0 mod D/3 / ℤ/3 sub-family" was that artifact.** The real signal is k=0 (the partner's left eigenvector is *nearly gauge-invariant*, no twist — residual 9.6e-4 at L=3). **P2's FORK(b) verdict (no exact factorization) STANDS** (k=0 gives 9.6e-4 ≠ 0); only the ℤ/3-sub-family refinement is retracted. **P1 (home = carry tower), P3 (transfer table), P4 (q=7) are UNAFFECTED** (they don't use `gauge_resid`). Code fixed; C1 tests the compression directly and supersedes P2's gauge test.

## C1 — the compression (core deliverable): the θ-reduction is REFUTED
Distinct-partner survival: the compressed spectrum must contain an eigenvalue near the partner **distinct from c₀'s image** (c₀ is kinematic, always in S — it must not masquerade as the partner). Ladder of reductions e_ρ mod {3, 6, 9, D}:

| q | L | partner | (e mod 3, γ) | (e mod 6, γ) | (e mod 9, γ) | full (e_ρ, γ) |
|---|---|---|---|---|---|---|
| 3 | 2 | 0.346827 | rel 4e-1 ✗ | rel 2e-2 ✗ | — | **rel 1.7e-2 ✗** |
| 3 | 3 | 0.333236 | rel 7e-1 ✗ | rel 7e-1 ✗ | rel 5e-1 ✗ | **rel 3.2e-3 ✓** |
| 7 | 2 | 0.158414 | rel 6e-1 ✗ | rel 4e-1 ✗ | rel 2e-1 ✗ | **rel 3.1e-3 ✓** |

- **The (θ = e mod 3, γ) subspace does NOT capture the partner** — off by rel 0.5–0.7. The user's premise (from P2's k² artifact) is **refuted.** Coarser-than-full reductions all fail.
- **The FULL gauge-invariant (e_ρ, γ) subspace captures the partner approximately**, distinct from c₀: q=3 L=3 rel 3.2e-3, q=7 L=2 rel 3.1e-3. **At q=3 L=2 even the full compression is off (rel 1.7e-2)** — near the EP the partner's left eigenvector still carries genuine twist (e_a) content; this shrinks with L (→ 3.2e-3 at L=3), tracking the partner's increasing gauge-invariance.
- **Not an exact reduction** (C3). So the (e_ρ, γ) compression is an *approximate* handoff, improving with L.
- **HANDOFF dumped:** `outputs/compressed_q3_L{2,3}.tsv` (full (e_ρ, γ), dim 54 / 486) + **exact rationals** `outputs/compressed_q3_L2_exact.tsv` (entries like 65/189, 17/567, 8/567 — rational as pre-registered). Labeled by (e_ρ, γ), flagged approximate with the rel error.

## C2 — cascade lemma check: class-DEPENDENT (pre-registration refuted)
Per (θ, γ) class, the new-carry valuation law `P(v₃(γ')=j)` vs the pre-registered `2·3^{−(j+1)}` (tail piled at top):

| L | law | per-class max deviation (min / med / max) | verdict |
|---|---|---|---|
| 2 | [0.667, 0.222, 0.111] | 0.111 / 0.111 / **0.317** | class-DEPENDENT |
| 3 | [0.667, 0.222, 0.074, 0.037] | 0.037 / 0.037 / **0.106** | class-DEPENDENT |

- **The cascade is class-dependent** — the uniform `2·3^{−(j+1)}` law holds only in aggregate (the clean P3 transfer was an average). The phase-coupling lives partly in the cascade.
- **Survival rates are also class-dependent:** average out-weight per class ranges **0.25–0.50** (spread 0.25) at both L. So the coupling lives in *both* the cascade and the survival (T mod 3) rates — not in one place.

## C3 — lumpability diagnostic on (θ, γ)
Within-class spread of class-aggregated outgoing weights: **0.23–0.33** (L=2), **0.25–0.33** (L=3) — the **0.258 obstruction**, present globally and **on the partner-home (γ≠0) classes (max 0.33)**. **Not lumpable** — quantifies why the compression (C1) is only approximate.

## C4 — q=7 control: same reduction structure, different rate
The q=7 tower-partner (0.158414, rank 22) behaves **exactly like q=3**: the θ-reduction fails, the full (e_ρ, γ) subspace captures it (rel 3.1e-3). So the reduction structure is **q-independent** (needs the full e_ρ tower); the q=3-vs-q=7 difference is the **rate** (q=7 gapped, q=3 coalescing). Consistent with Probe P's P4.

## Adjudication / redirect
| probe | verdict |
|---|---|
| C1 | (θ=e mod 3, γ) REFUTED; partner needs full (e_ρ, γ), approximately (rel 3e-3 at L=3 / gapped q=7; 1.7e-2 near-EP L=2). Handoff = the (e_ρ,γ) compression (rational). |
| C2 | cascade class-DEPENDENT (uniform law refuted); survival rates class-dependent 0.25–0.50. |
| C3 | not lumpable (0.258 obstruction) incl partner classes ⇒ compression is approximate, not exact. |
| C4 | q=7 same structure ⇒ reduction is q-independent; difference is the rate. |

**Redirect for the hand-derivation:** the target is the **full gauge-invariant (e_ρ, γ)-chain** (dim D·q^L), with **class-dependent survival rates and a non-lumpable cascade** — NOT a compact 3-state θ-tower. The partner is the (approximate) Perron root of this chain; the finite-L twist correction (the gap between the (e_ρ,γ) compression and the true partner) shrinks with L and is itself part of the object.

## Not at stake
R1–R46, Phases 0/1/2a, Sessions 1–2, D1, F, H, J, L-A, L-B, Real-T1, R1's STOP, Probe P's P1/P3/P4. No `r_q` value changes; **no rate-law fit** (2.9e-3, 1.0e-4 untouched).

_Reporting discipline: the pre-registered ℤ/3 compression is reported REFUTED, not massaged — and its origin (P2's k² bug) is disclosed and the code fixed. The "distinct-partner" survival criterion was tightened after catching that c₀ (kinematic, always in S) sits within the coalescence gap and can masquerade as a surviving partner — the corrected criterion changed several "SURVIVES" to "no". C2's uniform-cascade pre-registration is refuted (class-dependent). C3's non-lumpability is reported as expected (the 0.258 obstruction). Exact-rational dump provided as pre-registered. No rate extrapolation._
