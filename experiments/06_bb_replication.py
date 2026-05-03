"""
Experiment 06 — Bonacorsi-Bordoni framework apples-to-apples replication

Implements the same Bayesian NB2-GLM with log link as B&B's paper, fits
five specifications inside that framework, and scores all of them with
posterior predictive log score and W1 distance on a held-out test set.

Specifications:
    B0: log(n) only                                                  (baseline)
    B1: log(n) + hierarchical RE on (n mod 8), u_j ~ N(0, sigma_u)   (their winning model)
    B2: log(n) + fixed effects on a_final at k=3                     (ours, equivalent resolution)
    B3: log(n) + fixed effects on a_final at k=6                     (ours, finer resolution)
    B4: log(n) + hierarchical RE on (n mod 8) + fixed effects on a_final at k=6
        (B1 + ours; if a_final captures structure, sigma_u in B4 should shrink toward 0)

Scoring:
    Posterior predictive log score:   sum_i log p(y_test_i | data)
    W1 distance:                      between empirical and posterior predictive samples

Stan does the sampling via cmdstanpy. NUTS is slow on N=10^7; we subsample
to N=500,000 (configurable via --N_train) for tractability. Relative model
ordering should be preserved at any reasonable subsample size; absolute
numbers will differ from full-N reported values.

CAVEAT: B1's "hierarchical RE on (n mod 8)" matches B&B's spec to my best
understanding from their paper abstract. If their actual hyperprior on
sigma_u differs, the B1 numbers will shift slightly. The structural
question — whether a_final captures the residue structure better than mod 8
— is robust to this.

Usage:
    python 06_bb_replication.py --N 10000000 --N_train 500000 --N_test 50000
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import nbinom, wasserstein_distance

sys.stdout.reconfigure(encoding="utf-8")


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
    return np.array([deterministic_prefix(r, M) for r in range(M)], dtype=np.int64)


def build_X_and_re(features, re_levels, n):
    """Build design matrix X and re_idx for given features. RE always present
    in data block (set to all-1 if no RE used) per Stan's design."""
    cols = [np.ones((len(n), 1))]
    for label, vals, kind in features:
        if kind == "continuous":
            cols.append(vals.reshape(-1, 1).astype(np.float64))
        else:
            uniq = np.unique(vals)
            for u in uniq[1:]:
                cols.append((vals == u).astype(np.float64).reshape(-1, 1))
    X = np.hstack(cols)
    return X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--N_train", type=int, default=500_000,
                    help="Subsample size for training (Bayesian NUTS is slow at full N).")
    ap.add_argument("--N_test", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--threads-per-chain", type=int, default=4)
    ap.add_argument("--iter-warmup", type=int, default=500)
    ap.add_argument("--iter-sampling", type=int, default=500)
    ap.add_argument("--smoke", action="store_true",
                    help="Smoke test: tiny subsample, 1 chain, 100+100 iter.")
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--only", type=str, default=None,
                    help="Comma-separated list to run only some specs (e.g. B0,B1,B3).")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here.parent / "data"
    out_dir = Path(args.out) if args.out else here.parent / "experiments_output" / "bb_replication"
    out_dir.mkdir(parents=True, exist_ok=True)

    from cmdstanpy import CmdStanModel
    from cmdstanpy.utils import cxx_toolchain_path
    cxx_toolchain_path()

    print(f"[load] N={args.N:,}", flush=True)
    df = pl.read_parquet(data_dir / f"main_N{args.N}.parquet")
    n = df["n"].to_numpy().astype(np.int64)
    sigma = df["sigma"].to_numpy().astype(np.int64)
    log_n = np.log(n.astype(np.float64))
    mod8 = (n % 8).astype(np.int32)
    mod64 = (n % 64).astype(np.int32)
    a_final_k3 = lookup_a_final(3)[mod8].astype(np.int32)
    a_final_k6 = lookup_a_final(6)[mod64].astype(np.int32)
    print(f"        {len(n):,} rows  sigma range [{sigma.min()}, {sigma.max()}]", flush=True)

    # Train/test split: random 50K test, then random N_train from remainder
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(n))
    test_sel = perm[:args.N_test]
    train_sel = perm[args.N_test:args.N_test + args.N_train]
    if args.smoke:
        train_sel = train_sel[:5_000]
        test_sel = test_sel[:1_000]
    print(f"[split] train={len(train_sel):,}  test={len(test_sel):,}  seed={args.seed}", flush=True)

    # Build feature views
    feats = {
        "log_n": log_n, "mod8": mod8, "mod64": mod64,
        "afk3": a_final_k3, "afk6": a_final_k6,
    }

    # Specifications: name -> (X feats list, re vals or None)
    specs = {
        "B0": dict(feats=[("log_n", log_n, "continuous")], re=None),
        "B1": dict(feats=[("log_n", log_n, "continuous")], re=mod8),                 # hierarchical RE on mod 8
        "B2": dict(feats=[("log_n", log_n, "continuous"), ("afk3", a_final_k3, "cat")], re=None),
        "B3": dict(feats=[("log_n", log_n, "continuous"), ("afk6", a_final_k6, "cat")], re=None),
        "B4": dict(feats=[("log_n", log_n, "continuous"), ("afk6", a_final_k6, "cat")], re=mod8),  # B3 + RE on mod 8
    }

    only = set(args.only.split(",")) if args.only else set(specs.keys())

    print(f"[stan] compiling model ...", flush=True)
    t0 = time.perf_counter()
    model = CmdStanModel(stan_file=here / "nb2_glm.stan",
                         cpp_options={"STAN_THREADS": "TRUE"})
    print(f"[stan] compile done in {time.perf_counter()-t0:.1f}s", flush=True)

    chains = 1 if args.smoke else args.chains
    iter_warmup = 100 if args.smoke else args.iter_warmup
    iter_sampling = 100 if args.smoke else args.iter_sampling
    threads = 2 if args.smoke else args.threads_per_chain

    results = []
    for name, spec in specs.items():
        if name not in only:
            continue
        print()
        print(f"=== Spec {name} ===", flush=True)

        X_full = build_X_and_re(spec["feats"], spec["re"], n)
        X_train = X_full[train_sel]
        X_test = X_full[test_sel]
        y_train = sigma[train_sel].astype(int).tolist()
        y_test = sigma[test_sel].astype(int).tolist()

        if spec["re"] is not None:
            re_train = spec["re"][train_sel]
            re_test = spec["re"][test_sel]
            uniq = np.unique(np.concatenate([re_train, re_test]))
            re_map = {v: i+1 for i, v in enumerate(uniq)}
            re_idx_train = np.array([re_map[v] for v in re_train], dtype=np.int32).tolist()
            re_idx_test = np.array([re_map[v] for v in re_test], dtype=np.int32).tolist()
            N_re_levels = int(len(uniq))
            use_re = 1
        else:
            re_idx_train = [1] * len(train_sel)
            re_idx_test = [1] * len(test_sel)
            N_re_levels = 1
            use_re = 0

        N_tr = len(train_sel)
        grainsize = max(64, N_tr // (threads * 64))  # per-chain reduce_sum uses `threads` workers, not chains*threads

        data = {
            "N_train": N_tr,
            "N_test": len(test_sel),
            "P": X_train.shape[1],
            "X_train": X_train.tolist(),
            "X_test": X_test.tolist(),
            "y_train": y_train,
            "y_test": y_test,
            "use_re": use_re,
            "N_re_levels": N_re_levels,
            "re_idx_train": re_idx_train,
            "re_idx_test": re_idx_test,
            "grainsize": grainsize,
        }

        print(f"  X_train={X_train.shape}  use_re={use_re}  N_re_levels={N_re_levels}  "
              f"grainsize={grainsize}", flush=True)
        print(f"  fitting {chains} chains x ({iter_warmup}+{iter_sampling}) iter ...", flush=True)
        t0 = time.perf_counter()
        fit = model.sample(
            data=data,
            chains=chains, parallel_chains=chains,
            threads_per_chain=threads,
            iter_warmup=iter_warmup, iter_sampling=iter_sampling,
            seed=args.seed,
            output_dir=str(out_dir / f"stan_{name}"),
            show_progress=False,
            adapt_delta=0.9, max_treedepth=12,
        )
        fit_dt = time.perf_counter() - t0
        print(f"  fit done in {fit_dt:.1f}s", flush=True)

        # Posterior predictive log score on test
        log_lik = fit.stan_variable("log_lik_test")  # (n_draws, N_test)
        # Per-test-point: log mean(exp(log_lik_draws)) — use logsumexp
        from scipy.special import logsumexp
        n_draws = log_lik.shape[0]
        lpd_per_obs = logsumexp(log_lik, axis=0) - np.log(n_draws)
        log_score = float(lpd_per_obs.sum())

        # Sample W1: for each test obs, sample one y from posterior predictive
        beta_draws = fit.stan_variable("beta")        # (n_draws, P)
        phi_draws = fit.stan_variable("phi")          # (n_draws,)
        if use_re == 1:
            u_draws = fit.stan_variable("u")          # (n_draws, N_re_levels)
        rng_w1 = np.random.default_rng(args.seed + 1)
        # Pick one draw per test obs
        draw_idx = rng_w1.integers(0, n_draws, size=len(test_sel))
        log_mu_test = np.einsum("ij,ij->i", X_test, beta_draws[draw_idx])
        if use_re == 1:
            re_test_arr = np.array(re_idx_test) - 1
            log_mu_test += u_draws[draw_idx, re_test_arr]
        mu_test_draws = np.exp(log_mu_test)
        phi_test_draws = phi_draws[draw_idx]
        # NB2: mean=mu, variance=mu+mu^2/phi → nbinom(n=phi, p=phi/(phi+mu))
        n_param = phi_test_draws
        p_param = phi_test_draws / (phi_test_draws + mu_test_draws)
        nb_samples = rng_w1.negative_binomial(n_param, p_param)
        w1 = float(wasserstein_distance(np.array(y_test), nb_samples))

        # sigma_u summary if applicable
        if use_re == 1:
            sigma_u_post = fit.stan_variable("sigma_u")
            sigma_u_mean = float(sigma_u_post.mean())
            sigma_u_q = np.quantile(sigma_u_post, [0.05, 0.5, 0.95])
        else:
            sigma_u_mean = None
            sigma_u_q = None

        # phi summary
        phi_mean = float(phi_draws.mean())

        n_params = X_train.shape[1] + 1 + (1 if use_re else 0) + (N_re_levels if use_re else 0)
        result = dict(
            spec=name, n_params=n_params, log_score=log_score, w1=w1,
            phi_mean=phi_mean, sigma_u_mean=sigma_u_mean,
            sigma_u_q=sigma_u_q.tolist() if sigma_u_q is not None else None,
            fit_time=fit_dt,
        )
        results.append(result)
        print(f"  log_score={log_score:,.1f}  W1={w1:.3f}  phi={phi_mean:.4f}  "
              f"sigma_u={sigma_u_mean if sigma_u_mean else 'N/A'}", flush=True)

    print()
    print("=== Summary ===")
    print(f"{'Spec':<6} {'#p':>4} {'log_score':>14} {'W1':>7} {'phi':>8} {'sigma_u':>12} {'time(s)':>8}")
    for r in results:
        sigma_u_str = f"{r['sigma_u_mean']:.4f}" if r['sigma_u_mean'] is not None else "  --"
        print(f"{r['spec']:<6} {r['n_params']:>4} {r['log_score']:>14,.1f} {r['w1']:>7.3f} "
              f"{r['phi_mean']:>8.4f} {sigma_u_str:>12} {r['fit_time']:>8.1f}")

    # Save to CSV
    summary = pl.DataFrame([
        {"spec": r["spec"], "n_params": r["n_params"],
         "log_score": r["log_score"], "w1": r["w1"], "phi_mean": r["phi_mean"],
         "sigma_u_mean": r["sigma_u_mean"] if r["sigma_u_mean"] is not None else float("nan"),
         "fit_time_s": r["fit_time"]}
        for r in results
    ])
    summary.write_csv(out_dir / "bb_replication_summary.csv")
    print()
    print(f"[save] {out_dir/'bb_replication_summary.csv'}")


if __name__ == "__main__":
    main()
