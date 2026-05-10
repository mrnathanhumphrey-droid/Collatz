"""Paper 1 §5 head-to-head — HMC validation at N=10^4.

Re-runs B0..B4 from Paper 1 Table (head-to-head) at N=10,000 under HMC
(parallel_chains mode, NOT unified-mode), with 8000 train / 2000 test split.
Also re-runs Pathfinder at the same N for an apples-to-apples Section-3
side-by-side comparison; the §5 Pathfinder numbers were at N=500K so are
not directly comparable in absolute log-score magnitude.

Reuses C:/Collatz/experiments/nb2_glm.stan and the lookup_a_final helper
from 06_bb_replication.py / 06b_bb_pathfinder.py.

Usage:
    python paper1_hmc_n10k_validation.py
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import polars as pl
from scipy.special import logsumexp
sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
COLLATZ = HERE.parent
DATA_DIR = COLLATZ / "data"
OUT_DIR = COLLATZ / "experiments_output" / "paper1_hmc_n10k"
FIG_DIR = COLLATZ / "figures"
DOCS_DIR = COLLATZ / "docs"
PARQUET_OUT = COLLATZ / "data" / "paper1" / "hmc_n10k_results.parquet"
for d in (OUT_DIR, FIG_DIR, DOCS_DIR, PARQUET_OUT.parent):
    d.mkdir(parents=True, exist_ok=True)

SEED_DATA = 20260509
N_DATA = 10_000
N_TRAIN = 8_000
N_TEST = 2_000


def deterministic_prefix(r, a0, max_steps=400):
    a, c, steps = a0, r, 0
    while a % 2 == 0 and steps < max_steps:
        if c % 2 == 0:
            a //= 2; c //= 2
        else:
            a *= 3; c = 3 * c + 1
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


def make_specs(n_arr):
    log_n = np.log(n_arr.astype(np.float64))
    mod8 = (n_arr % 8).astype(np.int32)
    mod64 = (n_arr % 64).astype(np.int32)
    a_final_k3 = lookup_a_final(3)[mod8].astype(np.int32)
    a_final_k6 = lookup_a_final(6)[mod64].astype(np.int32)
    return {
        "B0": dict(feats=[("log_n", log_n, "continuous")], re=None),
        "B1": dict(feats=[("log_n", log_n, "continuous")], re=mod8),
        "B2": dict(feats=[("log_n", log_n, "continuous"),
                          ("afk3", a_final_k3, "cat")], re=None),
        "B3": dict(feats=[("log_n", log_n, "continuous"),
                          ("afk6", a_final_k6, "cat")], re=None),
        "B4": dict(feats=[("log_n", log_n, "continuous"),
                          ("afk6", a_final_k6, "cat")], re=mod8),
    }


def prepare_data_block(spec, n_arr, sigma_arr, train_idx, test_idx):
    X = build_X(spec["feats"], n_arr)
    X_train = X[train_idx]; X_test = X[test_idx]
    y_train = sigma_arr[train_idx].astype(int).tolist()
    y_test = sigma_arr[test_idx].astype(int).tolist()
    if spec["re"] is not None:
        re_train = spec["re"][train_idx]
        re_test = spec["re"][test_idx]
        uniq = np.unique(np.concatenate([re_train, re_test]))
        re_map = {v: i + 1 for i, v in enumerate(uniq)}
        re_idx_train = [re_map[v] for v in re_train]
        re_idx_test = [re_map[v] for v in re_test]
        N_re_levels = len(uniq); use_re = 1
    else:
        re_idx_train = [1] * len(train_idx)
        re_idx_test = [1] * len(test_idx)
        N_re_levels = 1; use_re = 0
    grainsize = max(64, len(train_idx) // 64)
    data = {
        "N_train": int(len(train_idx)),
        "N_test": int(len(test_idx)),
        "P": int(X_train.shape[1]),
        "X_train": X_train.tolist(),
        "X_test": X_test.tolist(),
        "y_train": y_train, "y_test": y_test,
        "use_re": use_re, "N_re_levels": int(N_re_levels),
        "re_idx_train": re_idx_train, "re_idx_test": re_idx_test,
        "grainsize": grainsize,
    }
    return data, X_train.shape[1], use_re, N_re_levels


def log_score_from_log_lik(log_lik):
    """log_lik: (n_draws, N_test). Returns sum_i log mean_d exp(log_lik_di)."""
    n_draws = log_lik.shape[0]
    lpd_per_obs = logsumexp(log_lik, axis=0) - np.log(n_draws)
    return float(lpd_per_obs.sum())


def diagnostics_from_fit(fit):
    """Pull worst-case Rhat / ESS_bulk + divergent / treedepth / EBFMI counts."""
    summary = fit.summary()
    cols_lower = {c.lower(): c for c in summary.columns}
    rhat_col = cols_lower.get("r_hat", cols_lower.get("rhat"))
    essb_col = cols_lower.get("ess_bulk")
    if rhat_col is None or essb_col is None:
        rhat_max = float("nan"); ess_min = float("nan")
    else:
        rhat_max = float(summary[rhat_col].max(skipna=True))
        ess_min = float(summary[essb_col].min(skipna=True))
    md = fit.method_variables()
    div_count = int(md["divergent__"].sum()) if "divergent__" in md else 0
    treedepth = (md["treedepth__"] if "treedepth__" in md else None)
    if treedepth is not None:
        max_td = int(getattr(fit, "metadata", None).cmdstan_config.get("max_treedepth", 10)
                     if hasattr(fit, "metadata") and getattr(fit.metadata, "cmdstan_config", None)
                     else 10)
        td_count = int((treedepth >= max_td).sum())
    else:
        td_count = 0
    try:
        ebfmi = fit.method_variables().get("energy__")
        if ebfmi is not None:
            chains_arr = ebfmi
            ebfmi_per_chain = []
            for c in range(chains_arr.shape[1]):
                e = chains_arr[:, c]
                num = float(np.sum(np.diff(e) ** 2))
                den = float(np.var(e) * len(e))
                ebfmi_per_chain.append(num / den if den > 0 else float("inf"))
            ebfmi_min = float(min(ebfmi_per_chain))
        else:
            ebfmi_min = float("nan")
    except Exception:
        ebfmi_min = float("nan")
    return rhat_max, ess_min, div_count, td_count, ebfmi_min


def run_hmc(model, data, spec_name, seed):
    print(f"  [HMC] {spec_name}: starting...", flush=True)
    t0 = time.perf_counter()
    fit = model.sample(
        data=data,
        chains=4, parallel_chains=4,
        threads_per_chain=1,           # explicit: NOT unified-mode
        iter_warmup=1000, iter_sampling=1000,
        seed=seed,
        output_dir=str(OUT_DIR / f"hmc_{spec_name}"),
        show_progress=False,
        adapt_delta=0.9, max_treedepth=12,
    )
    wall = time.perf_counter() - t0
    rhat, essb, divs, tds, ebfmi = diagnostics_from_fit(fit)
    log_lik = fit.stan_variable("log_lik_test")
    ls = log_score_from_log_lik(log_lik)
    beta = fit.stan_variable("beta")
    beta_log_n_idx = 1
    beta_log_n_mean = float(beta[:, beta_log_n_idx].mean())
    beta_log_n_sd = float(beta[:, beta_log_n_idx].std())
    if data["use_re"] == 1:
        sigma_u = fit.stan_variable("sigma_u")
        su_mean = float(sigma_u.mean())
        su_q025, su_q975 = [float(q) for q in np.quantile(sigma_u, [0.025, 0.975])]
        su_full = sigma_u.copy()
    else:
        su_mean = float("nan"); su_q025 = float("nan"); su_q975 = float("nan")
        su_full = None
    print(f"  [HMC] {spec_name}: log_score={ls:,.2f}  Rhat<{rhat:.4f}  "
          f"ESS_bulk_min={essb:.0f}  divs={divs}  TDsat={tds}  "
          f"EBFMI_min={ebfmi:.3f}  sigma_u={su_mean:.4f}  wall={wall:.1f}s",
          flush=True)
    return {
        "log_score": ls, "rhat_max": rhat, "ess_bulk_min": essb,
        "divergent_count": divs, "treedepth_count": tds, "ebfmi_min": ebfmi,
        "sigma_u_mean": su_mean, "sigma_u_q025": su_q025, "sigma_u_q975": su_q975,
        "sigma_u_full": su_full,
        "beta_log_n_mean": beta_log_n_mean, "beta_log_n_sd": beta_log_n_sd,
        "wall_time_s": wall,
    }


def run_pathfinder(model, data, spec_name, seed):
    print(f"  [PF ] {spec_name}: starting...", flush=True)
    t0 = time.perf_counter()
    fit = model.pathfinder(
        data=data,
        num_paths=4,
        draws=1000,
        seed=seed,
        output_dir=str(OUT_DIR / f"pf_{spec_name}"),
        show_console=False,
    )
    wall = time.perf_counter() - t0
    # Compute log_lik_test in Python (Pathfinder doesn't run gen quants by default)
    beta = fit.stan_variable("beta")               # (D, P)
    phi = fit.stan_variable("phi")                 # (D,)
    if data["use_re"] == 1:
        u = fit.stan_variable("u")                 # (D, N_re_levels)
        sigma_u_arr = fit.stan_variable("sigma_u")
        su_mean = float(sigma_u_arr.mean())
        su_q025, su_q975 = [float(q) for q in np.quantile(sigma_u_arr, [0.025, 0.975])]
        su_full = sigma_u_arr.copy()
    else:
        u = None; su_mean = float("nan"); su_q025 = float("nan"); su_q975 = float("nan")
        su_full = None
    Xt = np.array(data["X_test"]); yt = np.array(data["y_test"])
    re_idx_test = np.array(data["re_idx_test"]) - 1
    log_mu = Xt @ beta.T                           # (N_test, D)
    if data["use_re"] == 1:
        log_mu = log_mu + u.T[re_idx_test, :]      # broadcast: (N_test, D)
    # NB2 log pmf for each (test_i, draw_d):
    # log p(y|mu, phi) using lgamma identity
    from scipy.special import gammaln
    mu = np.exp(log_mu)
    phi_b = phi[None, :]                           # (1, D)
    yt_b = yt[:, None]                             # (N_test, 1)
    log_lik = (gammaln(yt_b + phi_b) - gammaln(yt_b + 1) - gammaln(phi_b)
               + phi_b * (np.log(phi_b) - np.log(phi_b + mu))
               + yt_b * (np.log(mu + 1e-300) - np.log(phi_b + mu)))
    # log_lik shape: (N_test, D); we want (D, N_test) to match HMC convention
    log_lik = log_lik.T
    ls = log_score_from_log_lik(log_lik)
    print(f"  [PF ] {spec_name}: log_score={ls:,.2f}  sigma_u={su_mean:.4f}  "
          f"wall={wall:.1f}s", flush=True)
    return {"log_score": ls, "sigma_u_mean": su_mean,
            "sigma_u_q025": su_q025, "sigma_u_q975": su_q975,
            "sigma_u_full": su_full, "wall_time_s": wall}


def main():
    print(f"[load] N=10M source data", flush=True)
    df = pl.read_parquet(DATA_DIR / "main_N10000000.parquet")
    n_full = df["n"].to_numpy().astype(np.int64)
    sigma_full = df["sigma"].to_numpy().astype(np.int64)
    rng = np.random.default_rng(SEED_DATA)
    sub_idx = rng.choice(len(n_full), size=N_DATA, replace=False)
    n_arr = n_full[sub_idx]
    sigma_arr = sigma_full[sub_idx]
    perm = rng.permutation(N_DATA)
    train_idx = perm[:N_TRAIN]
    test_idx = perm[N_TRAIN:N_TRAIN + N_TEST]
    print(f"        subsample N={N_DATA}, train={N_TRAIN}, test={N_TEST}, "
          f"seed={SEED_DATA}", flush=True)
    print(f"        sigma range [{sigma_arr.min()}, {sigma_arr.max()}]", flush=True)

    from cmdstanpy import CmdStanModel
    print(f"[stan] compiling nb2_glm.stan with STAN_THREADS=TRUE...", flush=True)
    t0 = time.perf_counter()
    model = CmdStanModel(stan_file=HERE / "nb2_glm.stan",
                         cpp_options={"STAN_THREADS": "TRUE"})
    print(f"        compile done in {time.perf_counter() - t0:.1f}s", flush=True)

    specs = make_specs(n_arr)
    rows = []
    sigma_u_draws = {"hmc": {}, "pf": {}}

    for name, spec in specs.items():
        print(f"\n=== {name} ===", flush=True)
        data, P, use_re, N_re = prepare_data_block(spec, n_arr, sigma_arr,
                                                    train_idx, test_idx)
        n_params = P + 1 + (1 if use_re else 0) + (N_re if use_re else 0)

        hmc = run_hmc(model, data, name, seed=SEED_DATA)
        pf = run_pathfinder(model, data, name, seed=SEED_DATA)

        if hmc["sigma_u_full"] is not None:
            sigma_u_draws["hmc"][name] = hmc["sigma_u_full"]
        if pf["sigma_u_full"] is not None:
            sigma_u_draws["pf"][name] = pf["sigma_u_full"]

        rows.append({
            "spec": name, "n_params": n_params,
            "n_train": N_TRAIN, "n_test": N_TEST, "seed": SEED_DATA,
            "log_score_hmc": hmc["log_score"],
            "log_score_pathfinder_n10k": pf["log_score"],
            "sigma_u_mean_hmc": hmc["sigma_u_mean"],
            "sigma_u_q025_hmc": hmc["sigma_u_q025"],
            "sigma_u_q975_hmc": hmc["sigma_u_q975"],
            "sigma_u_mean_pathfinder_n10k": pf["sigma_u_mean"],
            "sigma_u_q025_pathfinder_n10k": pf["sigma_u_q025"],
            "sigma_u_q975_pathfinder_n10k": pf["sigma_u_q975"],
            "beta_log_n_mean": hmc["beta_log_n_mean"],
            "beta_log_n_sd": hmc["beta_log_n_sd"],
            "rhat_max": hmc["rhat_max"],
            "ess_bulk_min": hmc["ess_bulk_min"],
            "divergent_count": hmc["divergent_count"],
            "treedepth_count": hmc["treedepth_count"],
            "ebfmi_min": hmc["ebfmi_min"],
            "wall_time_s_hmc": hmc["wall_time_s"],
            "wall_time_s_pf": pf["wall_time_s"],
        })

    out_df = pl.DataFrame(rows)
    out_df.write_parquet(PARQUET_OUT)
    print(f"\n[save] {PARQUET_OUT}", flush=True)
    print(out_df.select(["spec", "n_params", "log_score_hmc",
                         "log_score_pathfinder_n10k",
                         "sigma_u_mean_hmc",
                         "sigma_u_mean_pathfinder_n10k",
                         "rhat_max", "ess_bulk_min",
                         "divergent_count", "treedepth_count",
                         "ebfmi_min", "wall_time_s_hmc"]))

    # Save sigma_u draws to npz for figure
    np.savez(OUT_DIR / "sigma_u_draws.npz",
             hmc_B1=sigma_u_draws["hmc"].get("B1", np.array([])),
             hmc_B4=sigma_u_draws["hmc"].get("B4", np.array([])),
             pf_B1=sigma_u_draws["pf"].get("B1", np.array([])),
             pf_B4=sigma_u_draws["pf"].get("B4", np.array([])))
    print(f"[save] {OUT_DIR / 'sigma_u_draws.npz'}", flush=True)


if __name__ == "__main__":
    main()
