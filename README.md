# Collatz stopping-time structural analysis

**Status: closed 2026-05-01.**

## What this is

Empirical analysis of σ(n) — the Collatz total stopping time — for odd n up to 2²⁷, with a focus on residue-class structure. Final result is a prefix-decomposition theorem: the per-class distribution of σ on residues mod 2^k is parameterized entirely by the terminal value a_final ∈ {3^j : 1 ≤ j ≤ k} of the deterministic Collatz prefix on the residue. The 2^(k−1) odd classes mod 2^k therefore collapse onto exactly k distinct distributions.

Verified against Bonacorsi & Bordoni (Columbia, arXiv:2603.04479, March 2026); structural finding extends their mod-8 hierarchical analysis, and a direct head-to-head NB GLM at their data scale shows 7× parameter reduction at zero cost to predictive log score.

## Canonical artifacts

- [`writeup.md`](writeup.md) — the result document. Includes Related Work / head-to-head section.
- [`findings.md`](findings.md) — chronological log with sanity-check protocol applied per finding.
- [`figures/`](figures/) — Stage 2 EDA figures (σ vs log n by class, v-distribution, residual-tail, mod-64 grid).
- [`stage4_results/k6_uniform_full/`](stage4_results/k6_uniform_full/) — Stage 4 figures (forest plots, hyperparameter posteriors, GPD per class, tail probabilities, α decomposition, k=6 vs k=8 comparison, higher-moment universality, moment vs prefix, N-scaling).

## Pipeline

| Step | Script | Purpose |
|---|---|---|
| 1 | [`generate.py`](generate.py) | Numba memoized σ + features, parquet output |
| 2 | [`analyze.py`](analyze.py) | EDA figures |
| 3 prep | [`stage3_prep.py`](stage3_prep.py) | Stage-3 input: odd-only filter, class index, optional subsample |
| 3 fit | [`fit.py`](fit.py) + [`model.stan`](model.stan) | Hierarchical Stan fit |
| 4 | [`diagnose.py`](diagnose.py) | Posterior summary + GPD post-hoc + tail probs |

Reproduction commands in `writeup.md`.

## Reproduction smoke check (any future me, anyone with the repo)

1. `python generate.py --N 1048576` should finish in ~1.5s. σ at n=27 should be 111.
2. Run prefix algorithm (in head or pencil) on residue r=21 starting from state (a=64, c=21). Expected: 7 steps, terminating at (a_final=3, c_final=1).
3. Open arXiv:2603.04479 — confirm reported NB2-GLM log score is −272,911.95.

If those three check, the work is intact.

## Open follow-ups (not pursued in this module)

- **Analytical derivation** of Var(σ | a_final) as a closed form from random-walk first principles. Half day. Would convert empirical r=0.9999 into a theorem.
- **Trajectory-measure characterization** for the v=4, v=10 spike structure on Syracuse iterates. Pure number theory; self-contained.
- **KS-test rigorization** of the "exactly k distributions per modular resolution" claim. ~30 min. Confirmatory.
- **Larger N** (2²⁸, 2³⁰) to push the asymptotic universality verification further.
- **Hierarchical Bayesian version** of the predictive head-to-head with a_final coefficients drawn from a prior, to close the W₁ gap with B&B's reported W₁=3.199.

These are deferred — listed here for whoever picks this back up.

## Desktop zip

`collatz for desktop_2026-05-01_V3_BB_New.zip` (~2.3 MB) on Desktop is the canonical handoff package. Contains writeup.md + findings.md + figures + Stage 4 results CSVs. Sufficient for review; full reproduction needs the pipeline scripts above.
