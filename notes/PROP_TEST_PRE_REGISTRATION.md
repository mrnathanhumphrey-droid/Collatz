# Pre-Registration: Obstruction Propagation Test

**Pre-registered:** 2026-05-11T10:30 EDT, before any compute.
**Author:** Claude (acting as designated analyst).
**Hypothesis frame:** prompt provided by user, reconciled with framework definitions below.

---

## §0. Framing reconciliation (critical context)

The prompt asks whether "the c = 7/45 obstruction structure propagates under iteration of the prefix decomposition" and references "the c = 7/45 obstruction at residue class r mod 2^k."

After reading `writeup.md` (Result 3), `c_seven_forty_fifth.md`, `logical_chain_findings_to_c745.md` (Block 11 in particular), and `result_cycle_obstruction.md`, the framework as it stands defines two distinct objects:

1. **Prefix decomposition** (per-residue, 2-adic).
   For odd r mod 2^k, symbolic Collatz iteration on state (a, c) starting at (a = 2^k, c = r) — with rules
   - (a even, c even) → (a/2, c/2)
   - (a even, c odd) → (3a, 3c + 1)
   - terminate when a is odd
   — terminates at (a_final, c_final) with **a_final ∈ {3¹, 3², …, 3^k}** for some j(r) ∈ {1, …, k}. The 2^(k−1) odd residues mod 2^k partition into exactly k distinguishable values of a_final.
   Per-residue invariants produced by this map: **(a_final(r), c_final(r), j(r), prefix_steps(r))**.

2. **c = 7/45** (global, 3-adic).
   c = lim_{k→∞} S_k / 3 where S_k = Σ_{ξ ∈ Z/3^k, 3∤ξ} |μ̂_k(ξ)|² is the Plancherel mass of high-frequency Fourier coefficients of Tao's Syracuse random variable on Z/3^k. Also identified as the fixed point of T_diag = (1/5)·[[1, 1], [4, 4]] (R76.1).
   This is **a single scalar at level k → ∞**. It is not a per-residue-class-mod-2^k quantity. `c_seven_forty_fifth.md` characterizes it as a global asymptotic constant. `result_cycle_obstruction.md` (2026-05-05) concludes explicitly: "the framework characterizes ergodic asymptotic behavior, not finite-cycle structure" — and the same applies to per-class structure: c = 7/45 is not picking out specific residue classes.

**There is no object in the framework called "c = 7/45 obstruction signature at residue class r mod 2^k."** The closest framework-defined per-residue invariants are those produced by the prefix decomposition (#1 above).

**Decision:** rather than fabricate a quantity that doesn't exist, this pre-registration runs the closest tractable test that follows the *intent* of the prompt — does the per-residue structural signature propagate under iteration of the prefix decomposition? — with the *framework-defined* signature.

---

## §1. Reformulated hypothesis

Let **S(r) := (a_final(r), c_final(r), j(r), prefix_steps(r))** denote the prefix-decomposition signature of odd residue r mod 2^k.

**H_PROP:** The post-prefix re-entry n' = a_final(r) · m + c_final(r), reduced mod 2^k_target for k_target > k_base, has a structural signature S(n' mod 2^k_target) that depends on the original S(r) — with a predictable map φ: S(r) → distribution over S' values.

**H_NULL:** S(n' mod 2^k_target) under the pushforward over m is independent of S(r). The pushforward distribution at k_target is the same regardless of starting r.

**Pre-registered expectation: NULL.** Arithmetic preview justifying this prior (also pre-registered before compute):

> n' = 3^j(r) · m + c_final(r). Since gcd(3^j, 2^k_target) = 1, the map m ↦ n' mod 2^k_target is an affine bijection on Z/2^k_target. As m ranges over a full residue system mod 2^k_target with the appropriate parity (so n' is odd), n' mod 2^k_target hits *every* odd residue mod 2^k_target exactly once. Therefore the empirical distribution of S(n' mod 2^k_target) over m is the marginal distribution of S over all 2^(k_target − 1) odd residues mod 2^k_target — independent of r.

If this arithmetic preview holds, the empirical test will confirm NULL trivially. The test runs anyway, both to verify the arithmetic and to expose any residual structure (e.g., from finite-m sampling, parity restrictions, etc.).

---

## §2. Locked procedure

### Parameters (frozen before compute)

- **k_base = 6** — base modular resolution (32 odd residues). Matches the existing framework's primary characterization scale (Result 3, R²=0.9967 at k=6).
- **k_target = 12** — target resolution. Two doublings above k_base. 2^11 = 2048 odd residues at k_target.
- **m sample:** full enumeration m ∈ {0, 1, …, 2^k_target − 1} for each starting r at k_base. (k_target = 12 makes 4096 m-values × 32 r-values = 131,072 (r, m) pairs — well within compute budget.)
- **Signature definition:** S(r) := (a_final(r), c_final(r), j(r), prefix_steps(r)) computed by symbolic Collatz prefix iteration with terminating rule "a becomes odd."
- **Pushforward computation:** n'(r, m) := a_final(r) · m + c_final(r), modulo 2^k_target.
- **Iteration-2 procedure:** apply prefix decomposition again at k_target to n' to obtain S'(n' mod 2^k_target), then push forward once more to k_target2 = 18 via n''(r, m, m') := a_final(n') · m' + c_final(n') mod 2^k_target2.
- **Holdout for F4:** 20% of the 32 base residues (6 residues) selected by deterministic stride (every 5th odd residue starting from r=1) — held out from any "map fitting" before F4 prediction; locked seed: 20260511 for any random ops.

### Falsification cascade

- **F1 (structural correspondence):** for each starting r, compute the empirical distribution of S(n' mod 2^k_target) over the m-range. Test whether starting r predicts this distribution.
  Quantitative criterion: R² ≥ 0.85 on a regression of one informative scalar of the pushforward distribution (specifically, the mean of j(n' mod 2^k_target) across m) against scalar features of S(r). If R² < 0.85 → F1 fails.
- **F2 (second iteration):** repeat the propagation analysis one more level (n' → n''). If F1 passed, check whether the map composes with itself: does S(r) at k_base predict S(n'' mod 2^k_target2) at k_target2 with R² ≥ 0.85 under composed φ?
- **F3 (shuffle null):** randomly permute the assignment of S(r) to starting residues r (seed 20260511) and rerun F1 on the shuffled labels. If shuffle R² is comparable to actual R² → apparent correspondence is artifact.
  Quantitative criterion for "comparable": |R²(actual) − R²(shuffled)| < 0.10 means shuffled is comparable; ≥ 0.10 means shuffle breaks the correspondence.
- **F4 (cross-validation):** with the 6-residue holdout, fit any candidate map φ on the 26 in-fold residues, then predict pushforward distributions for the held-out 6. Verify predictions hold to within their in-fold residual scale.

### Adversarial safeguards

- **A1 (determinism):** rerun prefix decomposition on r = 1 and r = 27 twice; verify identical signatures both runs.
- **A2 (correctness):** for n = 27 (well-known long Collatz trajectory), spot-check the post-prefix re-entry against direct Collatz iteration. Verify residue tracking matches.
- **A3 (methodological consistency):** k_target computation uses the same symbolic prefix iteration code as k_base, with k_target substituted for k_base. No new code path.
- **A4 (deviation log):** any deviation from this pre-registration (parameter change, code change, scope change) must be documented in writing before continuing, with explicit reason.

### Decision rules (conservative; null-favored)

- **PROPAGATES_CLEAN:** F1 ∧ F2 ∧ F3 ∧ F4 all pass.
- **PROPAGATES_ONE_LEVEL:** F1 ∧ F3 ∧ F4 pass, F2 fails.
- **PROPAGATES_WEAK:** F1 R² ∈ [0.50, 0.85] AND F3 shuffle R² is comparable (< 0.10 gap) → artifact, not real propagation.
- **NO_PROPAGATION:** F1 fails (R² < 0.85). Document post-prefix distribution structure regardless of "interestingness."
- **INCONCLUSIVE:** any adversarial safeguard A1–A3 fails.

### What "the same as c = 7/45" would look like, if anything

Since c = 7/45 has no per-residue framework definition, there is no direct algebraic-consistency check (procedure step 4b). The closest possible test — and the only one I'll run for that step — is:

- Compute the global Plancherel mass S_k_target on Z/3^k_target via the existing characterization code, both before and after considering the pushforward distribution. Check whether the value 7/45 appears as a fixed point under the iteration in any framework-relevant way. If no such appearance exists → cleanly report that c = 7/45 plays no role at the per-residue propagation level.

This part of the test is exploratory and may produce a "no role found" outcome. That's fine.

---

## §3. Deliverables

- `PROP_TEST_PRE_REGISTRATION.md` — this document, committed before compute
- `PROP_TEST_RESULTS.md` — disposition at top, full report after compute
- `prop_test/prefix_signatures_k6.csv` — baseline S(r) for the 32 odd residues at k=6
- `prop_test/post_prefix_signatures_k12.csv` — S(n' mod 2^12) for all (r, m) pairs
- `prop_test/correspondence_summary.csv` — pushforward distribution summary per r
- `prop_test/shuffle_null_results.csv` — F3 shuffle test outputs
- `prop_test/iteration2_results.csv` — F2 second-iteration analysis
- `prop_test/holdout_validation.csv` — F4 cross-validation outputs
- `prop_test/src/` — source scripts (Python, all reproducible)

---

## §4. Compute budget

- Prefix decomposition at k=6 for 32 residues: ~ seconds
- Pushforward enumeration (131k (r,m) pairs) at k=12: ~ seconds
- Prefix decomposition at k=12 for 2048 odd residues: ~ minutes
- Falsification cascade: ~ minutes
- Total estimate: under 1 hour

Pre-registration locked.
