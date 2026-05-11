# Obstruction Propagation Test — Results

## Disposition

> **NO_PROPAGATION (NULL CONFIRMED)**

Pre-registered expectation: NULL. **The empirical result matches expectation to machine precision.**

The prefix-decomposition signature S(r) = (a_final, c_final, j, prefix_steps) does **not** propagate as a predictive structure under iteration of the prefix decomposition. The post-prefix re-entry n' = a_final(r) · m + c_final(r) modulo 2^k_target uniformizes signatures over the m-range, independent of starting r. The pushforward distribution at k_target is the marginal distribution of S over all odd residues mod 2^k_target — the same for every starting r.

This is a clean, mathematically grounded null. The empirical observation is corroborated by an explicit arithmetic argument (§4 below).

---

## Pre-registration adherence

- **Pre-registered:** 2026-05-11T10:30 EDT, committed at `49cf2ce` *before* any compute. See [PROP_TEST_PRE_REGISTRATION.md](PROP_TEST_PRE_REGISTRATION.md).
- **Procedure:** ran as locked. k_base = 6, k_target = 12, k_target2 = 18. Full m-enumeration at iteration 1 (2^12 = 4096 m-values per r). Sampled m, m' at iteration 2 (128 × 128 per r). Adversarial safeguards executed first.
- **No deviation** from pre-registered parameters or thresholds. No mid-run adjustments.

### Framing note (also in pre-reg §0)

The prompt's framing — "c = 7/45 obstruction at residue class r mod 2^k" — does not correspond to any object the framework defines. The framework has:

1. **Prefix decomposition** (per-residue mod 2^k, 2-adic): r → (a_final, c_final, j, prefix_steps).
2. **c = 7/45** (global Plancherel-mass asymptote on Z/3^k, 3-adic): a single scalar at level k → ∞.

There is no per-residue-mod-2^k c-value in the framework. The pre-registration substituted the closest framework-defined per-residue invariants (the prefix-decomposition signature) and tested *their* propagation under iteration. This is documented in pre-reg §0–§1; the substitution was locked before compute.

---

## Adversarial safeguards

| ID | Check | Result |
|---|---|---|
| **A1** | Determinism — rerun prefix_decompose on r ∈ {1, 27} twice at k=6 | ✓ identical |
| **A2** | Correctness — direct Collatz from n = 27 (m=0) and n = 91 (m=1) reaches predicted post-prefix value after the predicted step count | ✓ both match (n=27 → 107, n=91 → 350) |
| **A3** | Methodological consistency — same prefix_decompose code at k=6 and k=12, no separate code path | ✓ confirmed |
| **A4** | Deviation log — none required (no deviations from pre-reg) | ✓ |

---

## Falsification cascade

### F1 — Pushforward correspondence (the headline test)

For each odd r ∈ {1, 3, ..., 63} (32 residues), full enumeration of m ∈ {0, ..., 4095}, n'(r, m) = a_final(r) · m + c_final(r), kept the n' that are odd mod 2^12, looked up the signature S(n' mod 4096) at k=12, and computed pushforward summary statistics.

**Result table** (target: predict pushforward statistics from starting-r features):

| Target | Feature | R² | Slope |
|---|---|---:|---:|
| mean j' | j(r) | **0.000** | 0.000 |
| mean j' | prefix_steps(r) | **0.000** | 0.000 |
| mean j' | log(a_final(r)) | **0.000** | 0.000 |
| mean log a_final' | j(r) | **0.000** | — |
| mean log a_final' | log(a_final(r)) | **0.000** | — |

- mean j' across r: mean = 6.50000000, min = 6.50000000, max = 6.50000000, **spread = 0**.
- mean log(a_final') across r: 7.14098000 exactly across all 32 r.
- Full a_final' distribution vector across r-pairs: max L1 deviation = **0.0**.

Every starting r produces literally the same pushforward distribution at k=12. This is exact, not finite-sample-noisy.

**F1 disposition: FAIL** (R² ≈ 0, threshold was ≥ 0.85). The signature does not propagate.

### F3 — Shuffle null

100 random permutations of the (r → S(r)) mapping; rerun F1 on shuffled labels.

- R²(mean_j' ~ shuffled j_r): mean = 0.000, max = 0.000 across 100 trials.
- R²(mean_j' ~ shuffled log_af): mean = 0.000, max = 0.000.

Shuffle produces the same R² as actual. **Consistent with NULL** — there is no real correspondence to disrupt. (If F1 had produced positive R², shuffle would have killed it. Here both are zero.)

**F3 disposition: NULL-CONSISTENT**.

### F2 — Second-iteration test (k_target → k_target2 = 18)

Sampled 128 m-values per starting r, 128 m'-values per resulting n', computed n'' = a_final(n') · m' + c_final(n') mod 2^18, looked up S(n'' mod 2^18).

- mean j_iter2 across r: mean = 9.5048, min = 9.4605, max = 9.5607, spread = **0.100**.
- The full-enumeration limit value is 1 + (k_target2 − 1)/2 = 1 + 17/2 = 9.5 (mean of binomial j-distribution at k=18). All 32 r are within ~0.06 of 9.5 — the spread is sampling noise from finite m, m' samples.
- R²(mean_j_iter2 ~ j_r) = **0.099**. R²(mean_j_iter2 ~ log_af_r) = **0.099**.

Far below the 0.85 threshold. The tiny non-zero R² is finite-sample noise variance, not propagation structure.

**F2 disposition: FAIL** (formally; trivially expected since F1 failed).

### F4 — Holdout cross-validation

Holdout = 7 residues at stride 5 starting r=1: {1, 11, 21, 31, 41, 51, 61}.

In-fold OLS fit on 25 residues: `mean_j' = 6.500000 + 0.000e+00 · j_r`. Slope = 0 because target is constant.

Holdout predictions: all hit `actual = predicted = 6.500000` exactly. **Max error = 0.0**.

**F4 disposition: PASS** (trivially — constant target makes holdout prediction exact). This is "passes by uninformativeness," not by predictive substance.

### Cascade summary

| Check | Pass criterion | Actual | Pass? |
|---|---|---|---|
| F1 | R² ≥ 0.85 | 0.000 | **FAIL** |
| F2 | R² ≥ 0.85 at iter 2 | 0.099 | **FAIL** |
| F3 | shuffle R² < actual − 0.10 | 0 vs 0 (null-consistent) | n/a |
| F4 | holdout matches in-fold residual scale | exact | (vacuous PASS) |

Decision-rule branch: **F1 fails → NO_PROPAGATION**.

---

## §4. The structural characterization — why this is provably NULL

The empirical NULL is corroborated by elementary arithmetic. The pushforward map is

> n'(r, m) = a_final(r) · m + c_final(r),   reduced mod 2^k_target.

Since `a_final(r) = 3^j(r)` is odd, gcd(a_final, 2^k_target) = 1, so multiplication by a_final is a unit modulo 2^k_target. The map `m ↦ a_final · m + c_final mod 2^k_target` is therefore an affine bijection on Z / 2^k_target.

Parity: `a_final · m + c_final ≡ m + c_final (mod 2)` since a_final is odd. So `n'` is odd iff `m + c_final` is odd. For each starting r, there are exactly 2^(k_target − 1) = 2048 valid m-values producing odd n', and the map `m → n' mod 2^k_target` is a bijection from those m onto the 2048 odd residues of Z / 2^k_target.

**Consequence:** for every starting r, the multiset `{n'(r, m) mod 2^k_target : valid m}` is **exactly** the set of all odd residues mod 2^k_target, each hit once. Therefore the empirical distribution of S(n' mod 2^k_target) over the m-range is the marginal distribution of S over odd residues mod 2^k_target — **identical across r**. R² = 0 is forced by the arithmetic.

This holds at every iteration of the prefix decomposition: each successive iteration introduces another affine bijection by a unit (the new a_final), which again uniformizes signatures over the m-range. The structural signature cannot carry forward through iteration of the prefix decomposition.

---

## §5. The "c = 7/45 algebraic consistency" exploratory step

Pre-reg §2 listed this as exploratory: "Compute the global Plancherel mass S_k_target … check whether 7/45 appears as a fixed point under the iteration in any framework-relevant way."

**Result: no role found.** c = 7/45 is a 3-adic Plancherel asymptote on Z/3^k. The propagation question is about 2-adic per-residue structure on Z/2^k. The two objects live in different arithmetics and the framework does not currently connect them at this resolution.

What `result_cycle_obstruction.md` (2026-05-05) already established: "The framework characterizes ergodic asymptotic behavior, not finite-cycle structure." The current result is in the same family — c = 7/45 also does not characterize per-residue propagation under prefix iteration. This is not a contradiction or a surprise; it is consistent with the framework's known scope.

---

## §6. What this result actually says, structurally

The prefix decomposition is **memoryless under re-entry** in the following precise sense:

- The per-residue invariant (a_final, c_final, j, prefix_steps) captures complete information about the *symbolic prefix* from r mod 2^k.
- After the prefix terminates, the post-prefix integer m·a_final + c_final has its residue mod 2^k_target distributed uniformly over the appropriate-parity coset — **completely independent of where it came from**.
- Iterating the prefix decomposition is, in effect, restarting from a fresh uniform draw at the new resolution.

This is consistent with the framework's macro picture: the prefix decomposition explains the per-class intercept α(r) (Result 3, R² = 0.9996) as a one-shot algebraic offset; the post-prefix dynamics is universal random-walk-style behavior (slope β = 10.4282 universal, tail exponential). The propagation null we just confirmed is the precise statement that the prefix's per-class structural fingerprint *does not* persist into post-prefix integer dynamics at finer modular resolutions: once you've accounted for the prefix, the rest is universal.

**This is a strengthening of the existing framework, not a weakening.** It says: the structural per-class effects are exhaustively captured by a single application of the prefix decomposition — there is no hidden multi-scale structural cascade. If you wanted to find per-class structure beyond what the prefix already captures, iterating the prefix decomposition will not surface any.

---

## §7. Files

- [PROP_TEST_PRE_REGISTRATION.md](PROP_TEST_PRE_REGISTRATION.md) — locked rules, committed at 49cf2ce
- [PROP_TEST_RESULTS.md](PROP_TEST_RESULTS.md) — this document
- [prop_test/prefix_signatures_k6.csv](prop_test/prefix_signatures_k6.csv) — 32 odd residues × (a_final, c_final, j, prefix_steps)
- [prop_test/prefix_signatures_k12.csv](prop_test/prefix_signatures_k12.csv) — 2048 odd residues × signatures at k=12
- [prop_test/post_prefix_signatures_k12.csv](prop_test/post_prefix_signatures_k12.csv) — full (r, m) pushforward table, ~65k odd rows
- [prop_test/correspondence_summary.csv](prop_test/correspondence_summary.csv) — per-r pushforward summary statistics
- [prop_test/shuffle_null_results.csv](prop_test/shuffle_null_results.csv) — F3 shuffle trial outputs
- [prop_test/iteration2_results.csv](prop_test/iteration2_results.csv) — F2 iteration-2 sampled outputs
- [prop_test/holdout_validation.csv](prop_test/holdout_validation.csv) — F4 holdout predictions
- [prop_test/src/01_signatures.py](prop_test/src/01_signatures.py) — prefix decomposition + A1/A2 safeguards
- [prop_test/src/02_pushforward.py](prop_test/src/02_pushforward.py) — pushforward enumeration
- [prop_test/src/03_falsification_cascade.py](prop_test/src/03_falsification_cascade.py) — F1, F2, F3, F4

---

## §8. Disposition handling (per pre-reg)

The pre-registration's disposition-handling section said:

> **If NO_PROPAGATION:** cleanly documents that the c = 7/45 framework characterizes single-resolution structure only. Paper 4 stays at its current scope. The hypothesis is closed and the framework's limits are documented honestly.

The result lands here. The framework's **prefix decomposition** also turns out to characterize single-resolution structure only (it doesn't propagate under iteration). c = 7/45 separately remains a global 3-adic asymptote with no per-residue analog. These are now both documented as scope limits — clean honest results, not failures.

The mathematical content of the negative result is concrete: §4 shows that NO_PROPAGATION is *forced by arithmetic* (affine bijection by an odd unit modulo 2^k_target uniformizes signatures over the m-range). This is a theorem-level statement, not just an empirical observation — and stronger than the prompt's pre-registered NULL expectation in that it pins down *why* no propagation can occur.

If Paper 4's scope is structural decomposition + per-class characterization, this result tightens the claim: the structural decomposition captures all the per-class signal at modular resolution k_base, and iterating to higher resolutions does not surface any further structure. The prefix is a *complete* per-class invariant in the propagation-under-iteration sense.

**Compute log:** total elapsed under 5 minutes. Verifies the pre-reg compute budget estimate of "under 1 hour."
