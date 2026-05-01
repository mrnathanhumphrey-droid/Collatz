"""
Experiment 04 — Head-to-head Negative Binomial GLM (frequentist)

Replicates the Bonacorsi-Bordoni setup at N=10^7 with NB GLM via
statsmodels (frequentist MLE, NOT Bayesian; see caveats below) and
compares 5 specifications:

    M0 baseline:        log(n)
    M1 B&B-style:       log(n) + factor(n mod 8)
    M2 B&B-extended:    log(n) + factor(n mod 64)
    M3 ours k=3:        log(n) + factor(a_final at k=3)
    M4 ours k=6:        log(n) + factor(a_final at k=6)

For each: out-of-sample log score and W1 distance on a 50,000-observation
held-out test set (random selection, seed=42 for reproducibility).

CAVEATS — read before interpreting numbers:
  (1) statsmodels.discrete.NegativeBinomial is frequentist MLE, not
      Bayesian NB2-GLM with posterior predictive scoring. The B&B paper
      uses Bayesian inference. Posterior predictive log score is typically
      higher (less negative) than frequentist point-predictive log score
      because it integrates over posterior uncertainty in (beta, alpha).
  (2) W1 here is computed by drawing one NB sample per test point from the
      MLE-fit predictive; B&B compute W1 from posterior predictive samples
      averaged over draws. Different sampling protocol; W1 numbers are not
      directly comparable.
  (3) M3 and M4 use a_final defined via the prefix algorithm extended to
      all residues (including even r). For r=0 mod 2^k the algorithm
      halves all the way to a=1 then stops, giving a_final=1. For other
      even r it halves until c becomes odd then proceeds normally.

These results should be read as: "within the frequentist NB GLM family,
do a_final-derived covariates yield equivalent predictive accuracy with
fewer parameters than raw mod-K covariates?" Cross-framework comparison
to B&B's Bayesian results requires a separate Bayesian implementation
(see 06_bayesian_nb2_head_to_head.py if/when written).

Usage:
    python 04_head_to_head_nb_glm.py --N 10000000 --seed 42
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import nbinom, wasserstein_distance
from statsmodels.discrete.discrete_model import NegativeBinomial

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")


def deterministic_prefix(r, a0, max_steps=400):
    a, c = a0, r; steps = 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3*c + 1
        steps += 1
    return a


def lookup_a_final(k_pow):
    M = 1 << k_pow
    arr = np.zeros(M, dtype=np.int64)
    for r in range(M):
        arr[r] = deterministic_prefix(r, M)
    return arr


def build_X(features, n_rows):
    cols = [np.ones((n_rows, 1))]
    for label, vals, kind in features:
        if kind == "continuous":
            cols.append(vals.reshape(-1, 1).astype(np.float64))
        else:
            uniq = np.unique(vals)
            for u in uniq[1:]:  # K-1 dummies, baseline = first
                cols.append((vals == u).astype(np.float64).reshape(-1, 1))
    return np.hstack(cols)


def fit_nb_and_score(X_train, y_train, X_test, y_test, label, seed=0):
    t0 = time.perf_counter()
    print(f"  fitting {label} (X_train {X_train.shape})", flush=True)
    m = NegativeBinomial(y_train, X_train).fit(disp=0, maxiter=200)
    beta = m.params[:-1]
    alpha = max(m.params[-1], 1e-10)
    mu_test = np.clip(np.exp(X_test @ beta), 1e-6, 1e6)
    n_param = 1.0 / alpha
    p_param = n_param / (n_param + mu_test)
    log_pmf = nbinom.logpmf(y_test, n_param, p_param)
    log_score = float(log_pmf.sum())
    rng = np.random.default_rng(seed)
    nb_samples = rng.negative_binomial(n_param, p_param)
    w1 = float(wasserstein_distance(y_test, nb_samples))
    dt = time.perf_counter() - t0
    print(f"    -> log_score={log_score:,.1f}  W1={w1:.3f}  alpha={alpha:.4f}  ({dt:.1f}s)",
          flush=True)
    return dict(label=label, log_score=log_score, w1=w1, alpha=alpha,
                n_params=X_train.shape[1] + 1, fit_time=dt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n_test", type=int, default=50_000)
    ap.add_argument("--data", type=str, default=None)
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"

    lookup_k3 = lookup_a_final(3)
    lookup_k6 = lookup_a_final(6)
    print(f"a_final unique values: k=3 {sorted(set(lookup_k3.tolist()))}, "
          f"k=6 {sorted(set(lookup_k6.tolist()))}")

    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    sigma = df["sigma"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    mod8 = (n % 8).astype(np.int32)
    mod64 = (n % 64).astype(np.int32)
    a_final_k3 = lookup_k3[mod8].astype(np.int32)
    a_final_k6 = lookup_k6[mod64].astype(np.int32)
    print(f"data: {len(n):,} rows  (sigma min={sigma.min()}, max={sigma.max()})")

    rng = np.random.default_rng(args.seed)
    test_idx = rng.choice(len(n), args.n_test, replace=False)
    test_mask = np.zeros(len(n), dtype=bool); test_mask[test_idx] = True
    print(f"train: {(~test_mask).sum():,}  test: {test_mask.sum():,}  (seed={args.seed})")

    y_train = sigma[~test_mask]; y_test = sigma[test_mask]
    models = [
        ("M0 baseline: log(n)",                 [("log_n", log_n, "continuous")]),
        ("M1 B&B-style: log(n) + n mod 8",       [("log_n", log_n, "continuous"), ("mod8", mod8, "cat")]),
        ("M2 B&B-extended: log(n) + n mod 64",   [("log_n", log_n, "continuous"), ("mod64", mod64, "cat")]),
        ("M3 ours: log(n) + a_final(k=3)",       [("log_n", log_n, "continuous"), ("afk3", a_final_k3, "cat")]),
        ("M4 ours: log(n) + a_final(k=6)",       [("log_n", log_n, "continuous"), ("afk6", a_final_k6, "cat")]),
    ]
    results = []
    for label, feats in models:
        X = build_X(feats, len(n))
        results.append(fit_nb_and_score(X[~test_mask], y_train, X[test_mask], y_test, label,
                                        seed=args.seed))

    print()
    print("=== Final comparison ===")
    print(f"{'Model':<42} {'#p':>4} {'log_score':>14} {'W1':>7}")
    for r in results:
        print(f"{r['label']:<42} {r['n_params']:>4} {r['log_score']:>14,.1f} {r['w1']:>7.3f}")
    print()
    print("B&B reported (Bayesian NB2-GLM, mod-8 hierarchical, N=10^7, N_test=50K):")
    print(f"  log_score = -272,911.95   W1 = 3.199")
    print()
    print("CAVEAT: this is a frequentist comparison. B&B used Bayesian inference with")
    print("posterior predictive scoring. Cross-framework comparison requires Bayesian")
    print("implementation (Stan/PyMC) for parity.")


if __name__ == "__main__":
    main()
