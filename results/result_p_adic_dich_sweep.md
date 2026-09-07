# RESULT — p-ADIC-DICH-SWEEP: the 3-way dissociation across 9 primes. **SURPRISE — dissociation NOT established. Only the ENGINE generalizes (9/9); the FORWARD half does NOT (4/9) — the p=5,7 "forward lift" was small-n noise; the converse-leak has NO clean structure.** (2026-09-06)

**Probe:** `probes/probe_p_adic_dich_sweep.py` (reuses the validated `probe_p_adic_dich` machinery). Sweep
p ∈ {5,7,11,13,17,19,23,29,31}, q=2, chain `X=q^{−v}(1+p·a) mod p^r`. Adaptive exact depth (engine fold on
shallow exact Fractions; deeper γ tables on a float chain), `KMAX=64` so k=p and k=2p are in range for every
prime. Three claims reported as counts /9. White-baseline γ≡1 gate passed for all 9 primes.

## THE THREE COUNTS

| claim | count | prediction | |
|---|---|---|---|
| **(1) ENGINE** `Σ_j ρ_r(t+j·p^{r−1}) = ρ_{r−1}(t)` exact | **9/9** | 9/9 | ✓ **holds** |
| **(2) FORWARD** `p∣k ⟹ γ>1` AND `γ(p)<γ(2p)` | **4/9** | 9/9 | ✗ **FAILS the prediction** |
| **(3) CONVERSE FAILS** (nonempty leak set) | **9/9** | 9/9 | ✓ (converse never holds) |

**DISSOCIATION ESTABLISHED: FALSE.** It required engine 9/9 **AND** forward 9/9 **AND** converse-fails 9/9.
Forward is **4/9**. → **SURPRISE outcome.**

## Per-prime (q=2)

| p | RG | ENGINE | FORWARD | k=p enriched? | leak set (enriched non-multiples, k≤64) |
|---|---|---|---|---|---|
| 5 | 5 | OK | **OK** (5/5, g₅=1.29<g₁₀=2.83) | ✓ 1.288 | {1,4,6,11,14,16,19,21,31,34,36,39,41,44,56,61,64} |
| 7 | 5 | OK | **OK** (9/9, g₇=1.11<g₁₄=1.56) | ✓ 1.113 | {2,9,12,23,26,36,37,40,47,51,58} |
| 11 | 5 | OK | FAIL (mult 4/5) | ✓ 2.569 | {1,4,8,12,14,25,30,36,40,45,58,63} |
| 13 | 4 | OK | FAIL (mult 2/4, g₂₆=0.415 depleted) | ✓ 1.513 | {3,16,21,24,27,31,37,41,50,54,55,57} |
| 17 | 4 | OK | **OK** (3/3, g₁₇=1.19<g₃₄=3.51) | ✓ 1.185 | {2,6,7,13,24,28,30,32,38,44,50,55,62,64} |
| 19 | 4 | OK | FAIL (mult 1/3) | **✗ 0.994** | {17,18,24,33,36,44,45,59,61} |
| 23 | 4 | OK | FAIL (mult **0/2** — both depleted) | **✗ 0.264** | {6,10,17,52,53,60} |
| 29 | 3 | OK | FAIL (mult 1/2, g₅₈=0.553 depleted) | ✓ 3.941 | {5,11,24,25,51,54,56,62,63} |
| 31 | 3 | OK | **OK** (2/2, g₃₁=8.47<g₆₂=33.3) | ✓ 8.465 | {8,23,25,27,48,50} |

## What the sweep KILLS and what it KEEPS

**KEEP — the ENGINE is the one clean general law (9/9, exact, q-robust).** The fiber-tower telescoping
`Σ_j ρ_r(t+j·p^{r−1}) = ρ_{r−1}(t)` holds exactly (Fractions) at every prime, and (secondary, below) at every q
tested. Wilson's telescoping / p-divisible-shift argument is **general and robust**. This is the structural content
worth writing.

**KILL — the FORWARD implication does NOT generalize (4/9), and the p=5,7 "forward lift" was small-n noise.**
- Forward holds cleanly only at p=5,7,17,31.
- Even the *core* "k=p enriched" is only **7/9** — it FAILS at **p=19 (γ(19)=0.994, marginally depleted)** and
  **p=23 (γ(23)=0.264, strongly depleted; γ(46)=0.377 also depleted — BOTH multiples of 23 are depleted)**. At p=23
  the forward implication `p∣k ⟹ enriched` is outright **false**.
- Fiber-depth growth `γ(p)<γ(2p)` is also messy (fails at p=13, 19, 29; "holds" vacuously at p=23 where both
  values are depleted).
- Critically, **p=11 is the DEEPEST prime (RG=5) and still fails forward** — so this is not merely under-resolution.

**KILL — the CONVERSE-LEAK has NO clean description (honest null).** Across all 9 primes:
- The leak-residue-set mod p is **never** a union of ⟨2⟩-orbits (False ×9).
- It is **not** the QR set (e.g. p=7: leak-resid {1,2,5} vs QR {1,2,4}; p=11: {1,3,4,7,8} vs QR {1,3,4,5,9}).
- `|leak|` and `|leak-resid|` show **no law** in p or `(p−1)/ord_p(2)` (p=5:|leak|=17, p=23:6, p=31:6 — not monotone,
  not ∝ anything clean).

Characterization table:

| p | ord_p(2) | (p−1)/ord | \|leak-resid\| | \|leak\| | leak-resid = ⟨2⟩-orbit union? |
|---|---|---|---|---|---|
| 5 | 4 | 1 | 2 | 17 | False |
| 7 | 3 | 2 | 3 | 11 | False |
| 11 | 10 | 1 | 5 | 12 | False |
| 13 | 12 | 1 | 6 | 12 | False |
| 17 | 8 | 2 | 9 | 14 | False |
| 19 | 18 | 1 | 8 | 9 | False |
| 23 | 11 | 2 | 5 | 6 | False |
| 29 | 28 | 1 | 7 | 9 | False |
| 31 | 5 | 6 | 6 | 6 | False |

## Depth caveat (honest)
The larger primes ran at shallow exact-equivalent depth (RG=3 for p=29,31; RG=4 for p=13,17,19,23) with **many
sign-unstable k** (signs flip between the two deepest r — e.g. 7 unstable at p=11, 7 at p=31), so the *exact
leak-sets* are partly depth-noise and should not be over-read element-by-element. **But the headline is depth-robust:**
ENGINE is exact and passes 9/9 regardless of depth; FORWARD fails even at the deepest prime (p=11, RG=5) and fails
*hard* at p=23 (both multiples strongly depleted). Disentangling the *marginal* forward failures (p=19, p=13) from
under-resolution would need deeper exact r at those primes — a targeted follow-up, not run here.

## SECONDARY — q-variation (p=11,13 ; q=3,5)
- **ENGINE q-robust: OK at all four (p,q)** → engine holds at **13/13** (p,q) pairs tested. Strong.
- **FORWARD still fails** at every (p,q): p=11 q=3 (mult 4/5), q=5 (mult 4/5, fiber fails); p=13 q=3 (mult 1/4),
  q=5 (mult 2/4). Not a q artifact — forward is genuinely non-general.
- **Leak-set changes with q** (p=11: q=2→{1,4,8,12,14,25,30,36,40,45,58,63}, q=3→{1,21,23,29,32,54,60},
  q=5→12 elements) with no clean q-law — more evidence the leak is messy, not a fixed orbit structure.

## VERDICT and disposition
**SURPRISE — the dissociation is NOT established.** n>2 did exactly its job: it separated a structural law from
small-n noise, and the FORWARD half turned out to be noise (looked like a general lift at p=5,7; collapses to 4/9,
with a hard failure at p=23 and a failure even at the deepest prime p=11). The enrichment **dichotomy — both the
forward implication and the converse depletion — is 3-specific / non-general**, and the converse-leak has no clean
description.

**The ONE general theorem the sweep certifies is the ENGINE:** the fiber-tower telescoping
`Σ_j ρ_r(t+j·p^{r−1}) = ρ_{r−1}(t)` under a multiplier-equivariant p-divisible affine renewal — exact, 9/9 across
primes, 13/13 across (p,q). **This is what the #2 paper's structural claim should be** — the engine alone, framed
inside the ADD-A-PEAK range of `notes/CHANNEL_ENGINE_LIFT_RECON.md`. The enrichment dichotomy is a p=3 special
feature, **not** general content.

**Retractions (record):** (i) the q²=1+p miss-location mechanism (already retracted in `result_p_adic_dich.md`);
(ii) **now also the "forward-implication lifts" half of the p=5,7 dissociation** — the sweep shows it does not
generalize. What stood at n=2 as "engine + forward lift, converse 3-specific" reduces at n=9 to **"engine lifts;
the whole enrichment dichotomy is 3-specific."** Not at stake: CHANNELDICH p=3 (14/14), the tower/ENGINE identity
(now confirmed 9/9 general), R7–R10, MAHLER, value floor 0.473177.

---

## FOLLOW-UPS #2 (depth study) + #1 (stratified sample) — *is there a pattern to forward pass/fail, and is it depth or real?*

### #2 — DEPTH STUDY (`probes/probe_p_adic_dich_depth.py`): the forward failures are **depth-ROBUST, not artifacts**
Pushed the failing primes deeper and tabulated γ_r at the multiples of p:
- **p=11 (to r=6):** γ_r(11) = 1.370, 1.696, 1.805, 2.569, 1.298 (>1 at **every** r); γ_r(22) >1 at every r — but
  **γ_r(33=3p) = 1.254, 0.611, 0.462, 0.348, 0.352 — stably DEPLETED from r=3 on.** "All multiples enriched" is
  genuinely false: k=3p is depleted at every depth.
- **p=13 (to r=5):** γ_r(13) >1 (settling ~1.52); **γ_r(26=2p) = 1.938, 0.645, 0.415, 0.466 — stably DEPLETED**;
  γ_r(39=3p) also depleted (~0.68). The fiber pair genuinely fails.
- **p=23 (to r=4; r=5 memory-blocked ~1 GB):** **γ_r(23=p) = 0.580, 0.235, 0.264 — below 1 at every r, not rising
  toward 1.** The hard forward failure at k=p is real, not under-resolution.

⟹ **FORWARD is genuinely non-general.** And a finer fact emerges: it is not "no multiples enrich" — *specific*
multiples do (p=11: {p,2p,4p,5p} enriched, 3p depleted; p=13: {p,4p} enriched, {2p,3p} depleted; p=23: p,2p both
depleted). *Which* multiples enrich is itself structured, but messy.

### #1 — STRATIFIED SAMPLE (`probes/probe_p_adic_dich_strat.py`, primes 5..61): **no clean arithmetic law; the "primitive-root" tendency is CONFOUNDED with depth**
The n=9 tendency ("2 primitive root ⟹ forward fails", 7/9) does not survive as a law on 16 primes:

| | forward OK | forward FAIL |
|---|---|---|
| 2 primitive root | 1 (p=5) | 8 |
| 2 NOT prim root | 3 (p=7,17,31) | 4 (p=23,41,43,47) |

- "primroot ⟹ fail" strengthened to **8/9**, BUT "non-primroot ⟹ pass" **collapsed to 3/7** (p=23,41,43,47 also
  fail). So it is **not** a biconditional — many non-primroot primes fail too. Cosets `(p−1)/ord_p(2)` shows only a
  weak monotone hint (cosets=1: 1 OK/8 FAIL; cosets≥2: 3 OK/4 FAIL).
- **Forward holds at only 4 of 16 primes — all small: {5, 7, 17, 31}.** Everything p≥37 "fails."
- **BUT the large-prime failures are depth-CONFOUNDED (noise), not established.** Forced to RG=3 (memory:
  p=37 r=4 = 1.9 M support; r=5 = 69 M, infeasible), the k=p reading at large p is pure under-resolution:
  γ(p) ranges **0.009 (p=59) … 8.747 (p=61)** across adjacent primes — that spread is noise, not signal. So for
  p≥37 the forward test is simply *unresolved*, not cleanly failed.

### Combined verdict on the pattern question
- **Is there a clean law for forward pass/fail?** **No.** Among the *well-resolved* primes (RG≥4: {5,7,11,13,17,19,23})
  forward holds at {5,7,17} and fails at {11,13,19,23} — no primitive-root, coset-count, size, or residue rule
  separates them (p=5 primroot-OK vs p=11,13,19 primroot-FAIL; p=7,17 nonprimroot-OK vs p=23 nonprimroot-FAIL;
  non-monotone in size). The n=9 "primitive-root" pattern was a small-sample coincidence confounded with which
  primes happened to reach adequate depth.
- **Do we need a bigger sample?** **A bigger sample will not settle it** — the obstruction is depth, not n. Reading
  k=p reliably at p≥37 needs RG≥4, which is memory-infeasible (r=5 support runs 10⁷–10⁸). We are at the practical
  depth ceiling; more primes only add RG=3 noise.
- **What IS established:** (a) the ENGINE is general (9/9 exact, 13/13 with q, depth-independent); (b) FORWARD is
  genuinely non-general (depth-robust failures at p=11,13,23); (c) forward holds only for a few small primes with no
  law behind it. **The one clean, general, robust result remains the ENGINE alone.** The enrichment dichotomy —
  forward and converse — is a p=3 special feature with a small-prime tail, not general content.

### ORBIT-SIZE HYPOTHESIS (`probes/probe_p_adic_dich_orbit.py`) — *"does forward fail when the ⟨2⟩-orbit is filled?"* → **REFUTED**
A striking pattern appeared in the strat data: sorting the 16 primes by `ord_p(2)` (the ⟨2⟩-orbit size), the
"all-multiples" forward indicator split with **0 exceptions** — forward OK iff `ord_p(2) ≤ 8` (p=5,7,17,31: ords
4,3,8,5), FAIL iff `ord ≥ 10`. That is exactly the "fails when the orbit is filled (large)" shape. **But proper
verification killed it:**
- **Under a uniform *core* indicator** (just k=p and k=2p enriched + `γ(p)<γ(2p)`, dropping the size-confounded
  "all multiples in [1,64]"), the ord split **already dissolves** — large-ord primes p=37 (ord 36) and p=61
  (ord 60) pass the core, and OK/FAIL ords interleave ({3,4,5,7,8,10,35,36,60} vs {9,11,12,14,18,20,23,28,52,58}).
  The strat "0 exceptions" was an artifact of the all-multiples indicator, whose multiple-count itself scales with
  p (small p has 12 multiples in [1,64]; large p has 2) — i.e. confounded with prime size.
- **The clean decoupling test** — three large primes of matched size (~100) at matched depth r=3, varying only ord:

  | p | ord_p(2) | γ(p) | γ(2p) | forward |
  |---|---|---|---|---|
  | 127 | 7 | 32.4 | 129.3 | OK |
  | 73 | 9 | 0.79 | 0.096 | **FAIL** |
  | 71 | 35 | 3.19 | 11.85 | **OK** |

  **Non-monotone in ord:** the ord-35 prime (largest orbit) *passes* while the ord-9 prime *fails*. If forward
  tracked orbit size, ord-35 should fail hardest and ord-9 should pass — the opposite occurs. **Orbit size does not
  predict forward.**

⟹ **The orbit-size / `ord_p(2)` predictor is REFUTED** (the apparent law was a size-confounded artifact of the
all-multiples indicator). Combined with the earlier nulls (⟨2⟩-orbit residue unions, QR, primitive-root, coset
count), **no clean arithmetic predictor of forward pass/fail exists within reach** — forward is prime-idiosyncratic
(e.g. at matched depth/size, p=73 fails hard while its neighbours p=71, p=127 pass, for no orbit reason). The
decisive tests are also depth-bounded (the trio is only reachable at r≤3, and p=71's γ is still falling with r).
**The ENGINE (9/9 exact, depth-independent) remains the sole general, predictable law.**
