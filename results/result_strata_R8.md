# Probe R8 — uniform kill / ledger weld / sign front / strata — **ALL GATES PASS**

**Date:** 2026-07-21  Exact rationals. Probe `probes/probe_strata_R8.py` (reuses R7's engine verbatim; only a
measure swap + bookkeeping). Ordering load-bearing: R8-A ran first and passed, so R8-C/D are trusted.

## R8-A — THE UNIFORM KILL (hard gate): **PASS**
Replacing μ_{k−1} with the **uniform** measure on ℤ/3^{k−1} (everything else identical), **every C_k(m) = 0
exactly** for k = 2, 3, 4, 5 (all m in each table). No tolerance, exact-rational zero. Confirms Wilson's mechanism:
a ↦ 4^{−m}(1+3a) bijects the ≡1 mod 3 class, so Δ/3 is uniform ⟹ P = 3^{−(k−1)}, Q = 2·3^{−(k−1)}, C = 2P − Q = 0.
**The correlation lives entirely in the non-uniformity of μ_{k−1}.** Walk-back #31 not incurred; the Δ-decomposition
and ĉ normalization are sound. This gate greenlights R8-D and the LTE band derivation.

## R8-B — LEDGER–DEVIATION WELD: **HOLDS**
With d_K = S_K − 7/15 (frozen S-table): **Σ_{k=2}^K OffDiag_k = d_K − 1/5** exact, and **OffDiag_K = d_K − d_{K−1}**
per level, for K = 2, 3, 4, 5. Anchors reproduce: d₁ = 1/5 ✓, d₂ = 1/105 ✓, d₃ = −5191/1019445 ✓. R7-C's running
ledger and the #30 deviation fingerprint are one object, welded at both ends.

| K | Σ_{k≥2} OffDiag | d_K − 1/5 | OffDiag_K | d_K − d_{K−1} | weld |
|---|---|---|---|---|---|
| 2 | −4/21 | −4/21 | −4/21 | −4/21 | ✅ |
| 3 | −41816/203889 | −41816/203889 | −2980/203889 | −2980/203889 | ✅ |
| 4 | (exact) | = | +5699915795296300/… | = | ✅ |
| 5 | (exact) | = | +6958280182844849…/… | = | ✅ |

## R8-C — SIGN-FLIP FRONT (measurement, verbatim; NO fit, NO threshold detector)
**Verdict: there is NO monotone raw-m front m*(k). The sign is organized by stratum j = v₃(m) and by a
palindromic residue pattern within strata.** The signed tables, verbatim:

| k | sign by raw m (m=1…P) |
|---|---|
| 2 | `--+` |
| 3 | `+--++--++` |
| 4 | `++-+---+--++++++--+---+-+++` |
| 5 | `+-----+-+-++----+++++--+---++-+++-+-++++++++-+-+++-++---+--+++++----++-+-+-----++` |

**By stratum (the organizing structure):**
- **DC class (j = k−1, the m≡0 mod P self-orbit): always `+`** — C₂=2, C₃=10/7, C₄=94110/67963, C₅=… (all positive). ✓ anchor.
- **Top bulk stratum (j = k−2, just below DC): uniformly `−`** at every k — k=3 j=1 `--`, k=4 j=2 `--`, k=5 j=3 `--`.
- **Lower strata (j ≤ k−3): mixed, palindromic** — the sign string within each stratum reads the same forwards
  and backwards (consequence of R8-D2's C_k(r)=C_k(P−r)).

Anchors verified: **C₂(1) = −1**, **C₃(1) = +4/49** (the `+` opening k=3), **all DC entries positive**.
Reported verbatim; the pen adjudicates the shape. (No raw-m threshold imposed — R13/R15/R18 trap avoided.)

## R8-D1 — BAND COUNTS (gates the LTE derivation): **PASS**
Direct integer pair-count over (a,a') ∈ (ℤ/3^{k−1})²: **#{v₃(Δ) ≥ k} = 3^{k−1}** and **#{v₃(Δ) = k−1} = 2·3^{k−1}**,
**identical for every m** (observed count-set is a singleton at each k). k=3 → (9, 18); k=4 → (27, 54); k=5 →
(81, 162). Wilson's m-independent affine band is confirmed exactly: the band's pair-count does not depend on m —
the m-dependence of C_k lives purely in *which* pairs land in the band (via the measure weights), not *how many*.

## R8-D2 — PALINDROME: **PASS**
C_k(r) = C_k(3^{k−1} − r) exact for all r ≢ 0, verified k = 4, 5 (k = 2, 3 prior). The m → −m (gap-sign) symmetry.

## R8-D3 — STRATUM WEIGHTS + AVERAGED CORRELATIONS (the pen's raw material; NO fit)
OffDiag_k decomposes cleanly by stratum: **OffDiag_k = (2/3) Σ_j W_j · C̄_k(j)**, verified = frozen OffDiag_k for
every k = 2…5. Two routes for the stratum weight agree exactly (labeled columns, R17 lesson honored):

- **residue-grouped:** W_j = Σ_{r mod P, v₃(r)=j} 4^{−r}/(1−4^{−P})
- **closed form:** W_j = x/(1−x) − x³/(1−x³), x = 4^{−3^j} (derived from the geometric stratum structure)

**Structural fact for the limit law: W_j is k-INDEPENDENT** — a pure geometric weight over the m-strata. W₀ = 20/63,
W₁ = 4160/262143, W₂ = 68719738880/180143…, W₃ = 324518553658426744…, identical across all k rows. **The entire
k-dependence of OffDiag_k — and therefore the whole limit law — lives in the stratum-averaged correlations C̄_k(j).**

| k | C̄_k(0) | C̄_k(1) | C̄_k(2) | DC C̄_k(k−1) | DC weight×value |
|---|---|---|---|---|---|
| 2 | −1 | — | — | +2 | +2/63 |
| 3 | **−324/9709** | −5/7 | — | +10/7 | +5.45e−6 |
| 4 | **+ (≈+0.026)** | −6365212470/229065 | −47055/67963 | +94110/67963 | +7.69e−17 |
| 5 | **+ (≈+0.021)** | + | + | +429586948978369470/… | +2.38e−49 |

**⭐ The overshoot mechanism, refined and located:** the dominant stratum weight is **W₀ = 20/63** (≈0.317, an
order of magnitude above W₁). **C̄_k(0) — the j=0 stratum-averaged correlation — crosses zero between k=3 and k=4**
(−1 → −324/9709 → positive), and that crossing *is* the OffDiag sign flip (−,−,+,+). R7 located the flip in the
single channel C_k(1); R8 shows the responsible object is the whole dominant stratum C̄_k(0) (C_k(1) is one of its
54 members at k=5). The DC/self-orbit class stays positive but its weight collapses as 1/(4^P−1) (Mersenne):
5.45e−6 (k=3) → 7.7e−17 (k=4) → 2.4e−49 (k=5) — numerically dead by k=3, exactly as pre-registered (≲1e−5).

## Status
**R8 ALL GATES PASS.** Uniform kill (A) confirms the correlation is purely a non-uniformity effect and greenlights
the LTE band derivation; the ledger↔deviation weld (B) is exact at both ends; the band pair-count (D1) is
m-independent as Wilson's affine-band claim requires; the palindrome (D2) holds through k=5. The sign front (C) is
**not** a raw-m front — it is stratum-organized (DC always +, top bulk stratum always −, lower strata palindromic).
The strata decomposition (D3) hands the pen its raw material: **W_j is k-independent (closed form x/(1−x)−x³/(1−x³),
x=4^{−3^j}), so the limit law is entirely C̄_∞(j)**; and the −,−,+,+ overshoot is **C̄_k(0) crossing zero between
k=3 and k=4**. **Still owed (pen):** the stationary profile C̄_∞(j) and Σ_k (2/3)Σ_j W_j C̄_k(j) → −1/5 in closed
form (Theorem S's L→∞ limit). No fitting; exact rationals throughout; nothing smoothed.
