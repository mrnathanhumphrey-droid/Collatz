"""
Experiment 06b — Pathfinder version of the BB replication.

Same model + likelihood as 06_bb_replication.py, but uses Stan's Pathfinder
algorithm (Zhang et al. 2022, JMLR) instead of HMC/NUTS. Pathfinder is a
single-pass quasi-Newton variational approximation that avoids HMC warmup
and runs orders of magnitude faster on this scale of data.

Pivoted to this after hitting Stan 2.36+ unified-mode multi-chain lockup
at full N=500K (4 chains in shared-thread-pool process locked all 16 threads
at full CPU with zero sample progress for 60+ min).

Same scoring: posterior predictive log score (LPD) and W1 distance on
held-out test set. Pathfinder log_lik_test computed in Python because
generated_quantities don't run during pathfinder by default.

Usage:
    python 06b_bb_pathfinder.py --N 10000000 --N_train 500000 --N_test 50000
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import nbinom, wasserstein_distance
from scipy.special import logsumexp

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


def build_X(features, n):
    cols = [np.ones((len(n), 1))]
    for label, vals, kind in features:
        if kind == "continuous":
            cols.append(vals.reshape(-1, 1).astype(np.float64))
        else:
            uniq = np.unique(vals)
            for u in uniq[1:]:
                cols.append((vals == u).astype(np.float64).reshape(-1, 1))
    return np.hstack(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10_000_000)
    ap.add_argument("--N_train", type=int, default=500_000)
    ap.add_argument("--N_test", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--num-paths", type=int, default=4)
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--threads-per-path", type=int, default=8)
    ap.add_argument("--data", type=str, default=None)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--only", type=str, default=None)
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
    print(f"        {len(n):,} rows  sigma range [{sigma.min()}, {sigma.max()}]", flush=True)

    a_final_k3_lookup = lookup_a_final(3)
    a_final_k6_lookup = lookup_a_final(6)
    mod8_full = (n % 8).astype(np.int32)
    mod64_full = (n % 64).astype(np.int32)
    a_final_k3 = a_final_k3_lookup[mod8_full].astype(np.int32)
    a_final_k6 = a_final_k6_lookup[mod64_full].astype(np.int32)

    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(len(n))
    test_sel = perm[:args.N_test]
    train_sel = perm[args.N_test:args.N_test + args.N_train]
    print(f"[split] train={len(train_sel):,}  test={len(test_sel):,}  seed={args.seed}", flush=True)

    print(f"[stan] compiling model ...", flush=True)
    t0 = time.perf_counter()
    model = CmdStanModel(stan_file=here / "nb2_glm.stan",
                         cpp_options={"STAN_THREADS": "TRUE"})
    print(f"[stan] compile done in {time.perf_counter()-t0:.1f}s", flush=True)

    specs = {
        "B0": dict(feats=[("log_n", log_n, "continuous")], re=None),
        "B1": dict(feats=[("log_n", log_n, "continuous")], re=mod8),
        "B2": dict(feats=[("log_n", log_n, "continuous"), ("afk3", a_final_k3, "cat")], re=None),
        "B3": dict(feats=[("log_n", log_n, "continuous"), ("afk6", a_final_k6, "cat")], re=None),
        "B4": dict(feats=[("log_n", log_n, "continuous"), ("afk6", a_final_k6, "cat")], re=mod8),
    }
    only = set(args.only.split(",")) if args.only else set(specs.keys())

    results = []
    for name, spec in specs.items():
        if name not in only:
            continue
        print(f"\n=== Spec {name} ===", flush=True)
        feats = spec["feats"]
        re_vals = spec["re"]

        X_full = build_X(feats, n)
        X_train = X_full[train_sel]
        X_test = X_full[test_sel]
        y_train = sigma[train_sel].tolist()
        y_test = sigma[test_sel].tolist()

        if re_vals is not None:
            uniq_re = np.unique(re_vals)
            re_map = {int(u): i + 1 for i, u in enumerate(uniq_re)}
            re_idx_train = [re_map[int(re_vals[i])] for i in train_sel]
            re_idx_test = [re_map[int(re_vals[i])] for i in test_sel]
            N_re_levels = len(uniq_re)
            use_re = 1
        else:
            re_idx_train = [1] * len(train_sel)
            re_idx_test = [1] * len(test_sel)
            N_re_levels = 1
            use_re = 0

        N_tr = len(train_sel)
        grainsize = max(64, N_tr // (args.threads_per_path * 64))

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
        print(f"  X_train={X_train.shape}  use_re={use_re}  N_re_levels={N_re_levels}  grainsize={grainsize}", flush=True)
        print(f"  pathfinder paths={args.num_paths} draws={args.draws} num_threads={args.threads_per_path}", flush=True)

        t0 = time.perf_counter()
        fit = model.pathfinder(
            data=data,
            num_paths=args.num_paths,
            draws=args.draws,
            seed=args.seed,
            output_dir=str(out_dir / f"path_{name}"),
            show_console=False,
            num_threads=args.threads_per_path,
        )
        fit_dt = time.perf_counter() - t0
        print(f"  fit done in {fit_dt:.1f}s", flush=True)

        # Compute log_lik_test in Python (pathfinder doesn't run generated_quantities by default)
        beta_draws = fit.stan_variable("beta")        # (n_draws, P)
        phi_draws = fit.stan_variable("phi")          # (n_draws,)
        if use_re == 1:
            u_draws = fit.stan_variable("u")          # (n_draws, N_re_levels)
        n_draws = beta_draws.shape[0]
        print(f"  pathfinder draws: {n_draws}", flush=True)

        # log_mu_test: (N_test, n_draws)
        log_mu = X_test @ beta_draws.T  # (N_test, n_draws)
        if use_re == 1:
            re_test_arr = np.array(re_idx_test) - 1  # 0-indexed
            log_mu += u_draws[:, re_test_arr].T
        mu = np.exp(log_mu)
        phi_b = phi_draws[None, :]  # (1, n_draws)
        y_test_arr = np.array(y_test)[:, None]
        log_lik = nbinom.logpmf(y_test_arr, n=phi_b, p=phi_b / (phi_b + mu))  # (N_test, n_draws)
        # Posterior predictive LPD: per obs, log mean(exp(log_lik over draws))
        lpd_per_obs = logsumexp(log_lik, axis=1) - np.log(n_draws)
        log_score = float(lpd_per_obs.sum())

        # W1: sample one y_pred per test obs from posterior predictive
        rng_w1 = np.random.default_rng(args.seed + 1)
        draw_idx = rng_w1.integers(0, n_draws, size=len(test_sel))
        log_mu_test_one = np.einsum("ij,ij->i", X_test, beta_draws[draw_idx])
        if use_re == 1:
            log_mu_test_one += u_draws[draw_idx, re_test_arr]
        mu_test_draws = np.exp(log_mu_test_one)
        phi_test_draws = phi_draws[draw_idx]
        n_param = phi_test_draws
        p_param = phi_test_draws / (phi_test_draws + mu_test_draws)
        nb_samples = rng_w1.negative_binomial(n_param, p_param)
        w1 = float(wasserstein_distance(y_test_arr.flatten(), nb_samples))

        # phi / sigma_u summaries
        phi_mean = float(phi_draws.mean())
        if use_re == 1:
            sigma_u_post = fit.stan_variable("sigma_u")
            sigma_u_mean = float(sigma_u_post.mean())
            sigma_u_q = np.quantile(sigma_u_post, [0.05, 0.5, 0.95]).tolist()
        else:
            sigma_u_mean = None
            sigma_u_q = None

        n_params = X_train.shape[1] + 1 + (1 if use_re else 0) + (N_re_levels if use_re else 0)
        result = dict(
            spec=name, n_params=n_params, log_score=log_score, w1=w1,
            phi_mean=phi_mean, sigma_u_mean=sigma_u_mean, sigma_u_q=sigma_u_q,
            fit_time=fit_dt, n_draws=n_draws,
        )
        results.append(result)
        print(f"  log_score={log_score:,.1f}  W1={w1:.3f}  phi={phi_mean:.4f}  "
              f"sigma_u={sigma_u_mean if sigma_u_mean is not None else 'N/A'}", flush=True)

    print()
    print("=== Summary (Pathfinder) ===")
    print(f"{'Spec':<6} {'#p':>4} {'log_score':>14} {'W1':>7} {'phi':>8} {'sigma_u':>12} {'time(s)':>8}")
    for r in results:
        sigma_u_str = f"{r['sigma_u_mean']:.4f}" if r['sigma_u_mean'] is not None else "  --"
        print(f"{r['spec']:<6} {r['n_params']:>4} {r['log_score']:>14,.1f} {r['w1']:>7.3f} "
              f"{r['phi_mean']:>8.4f} {sigma_u_str:>12} {r['fit_time']:>8.1f}")

    summary = pl.DataFrame([
        {"spec": r["spec"], "n_params": r["n_params"],
         "log_score": r["log_score"], "w1": r["w1"], "phi_mean": r["phi_mean"],
         "sigma_u_mean": r["sigma_u_mean"] if r["sigma_u_mean"] is not None else float("nan"),
         "fit_time_s": r["fit_time"], "n_draws": r["n_draws"]}
        for r in results
    ])
    summary.write_csv(out_dir / "bb_pathfinder_summary.csv")
    print(f"[save] {out_dir/'bb_pathfinder_summary.csv'}")


if __name__ == "__main__":
    main()
