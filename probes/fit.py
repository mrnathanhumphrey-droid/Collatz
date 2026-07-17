"""
Stage 3 fit: cmdstanpy runner for the hierarchical Normal Collatz model.

Usage:
  # Smoke (compile + 1 chain * 100 iter on 50K rows, ~1-2 min total):
  python fit.py --input data/stage3_input_N33554432_k10.parquet --smoke

  # Full fit (4 chains * (warmup+sampling) iter, threaded):
  python fit.py --input data/stage3_input_N33554432_k10.parquet \
      --chains 4 --threads-per-chain 4 \
      --iter-warmup 500 --iter-sampling 500
"""

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, required=True,
                    help="Stage 3 input parquet (from stage3_prep.py).")
    ap.add_argument("--out", type=str, default=None,
                    help="Output dir for Stan CSVs. Default fits/<tag>/")
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--threads-per-chain", type=int, default=4)
    ap.add_argument("--iter-warmup", type=int, default=500)
    ap.add_argument("--iter-sampling", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--grainsize", type=int, default=0,
                    help="reduce_sum grainsize. 0 -> auto (N // (chains*threads*64)).")
    ap.add_argument("--smoke", action="store_true",
                    help="Smoke test: 50k subsample, 1 chain, 100 warmup+100 sampling.")
    ap.add_argument("--max-treedepth", type=int, default=12,
                    help="Stan max_treedepth (default 12; raise if saturated).")
    ap.add_argument("--adapt-delta", type=float, default=0.9,
                    help="Stan adapt_delta (default 0.9).")
    ap.add_argument("--tag", type=str, default=None,
                    help="Output dir suffix. Default infers from input filename.")
    args = ap.parse_args()

    # Lazy import so script can be inspected without cmdstanpy installed
    from cmdstanpy import CmdStanModel
    from cmdstanpy.utils import cxx_toolchain_path
    cxx_toolchain_path()  # adds bundled RTools40 mingw64+usr/bin to PATH

    here = Path(__file__).resolve().parent

    print(f"[load]  {args.input}", flush=True)
    df = pl.read_parquet(args.input)
    print(f"        {df.height:,} rows, {df.width} cols", flush=True)

    if args.smoke:
        n_smoke = min(50_000, df.height)
        df = df.sample(n=n_smoke, seed=args.seed)
        print(f"[smoke] subsetted to {df.height:,} rows", flush=True)

    K = int(df["class_idx"].max()) + 1
    N = df.height

    chains = 1 if args.smoke else args.chains
    iter_warmup = 100 if args.smoke else args.iter_warmup
    iter_sampling = 100 if args.smoke else args.iter_sampling
    threads_per_chain = 2 if args.smoke else args.threads_per_chain

    grainsize = args.grainsize
    if grainsize <= 0:
        grainsize = max(64, N // (chains * threads_per_chain * 64))

    data = {
        "N": N,
        "K": K,
        "class_idx": (df["class_idx"].to_numpy() + 1).astype(np.int32).tolist(),
        "log_n": df["log_n"].to_numpy().astype(np.float64).tolist(),
        "y": df["sigma"].to_numpy().astype(np.float64).tolist(),
        "grainsize": int(grainsize),
    }

    print(f"[data]  N={N:,}  K={K}  grainsize={grainsize}", flush=True)
    print(f"[fit]   {chains} chains x ({iter_warmup}+{iter_sampling}) iter, "
          f"{threads_per_chain} thr/chain", flush=True)

    print("[stan]  compiling model.stan ...", flush=True)
    t0 = time.perf_counter()
    model = CmdStanModel(
        stan_file=here / "model.stan",
        cpp_options={"STAN_THREADS": "TRUE"},
    )
    print(f"[stan]  compile done in {time.perf_counter() - t0:.1f}s",
          flush=True)

    out_root = Path(args.out) if args.out else here / "fits"
    if args.tag is None:
        # Infer tag from input filename, e.g. stage3_input_N..._k6_uniform.parquet -> k6_uniform
        stem = Path(args.input).stem
        infer = stem.replace("stage3_input_", "").lstrip("N0123456789_")
        tag = f"{infer}_N{N}_{'smoke' if args.smoke else 'full'}"
    else:
        tag = f"{args.tag}_{'smoke' if args.smoke else 'full'}"
    out_dir = out_root / tag
    out_dir.mkdir(parents=True, exist_ok=True)

    # Explicit inits to avoid the occasional inf during default uniform-init
    inits = {
        "mu_alpha": 0.0,
        "mu_beta": 10.4,
        "tau_alpha": 1.0,
        "tau_beta": 0.05,
        "alpha_z": [0.0] * K,
        "beta_z": [0.0] * K,
        "phi": 64.0,
    }

    t0 = time.perf_counter()
    fit = model.sample(
        data=data,
        chains=chains,
        parallel_chains=chains,
        threads_per_chain=threads_per_chain,
        iter_warmup=iter_warmup,
        iter_sampling=iter_sampling,
        seed=args.seed,
        output_dir=str(out_dir),
        show_progress=False,
        inits=inits,
        max_treedepth=args.max_treedepth,
        adapt_delta=args.adapt_delta,
    )
    fit_time = time.perf_counter() - t0
    print(f"[fit]   sampling done in {fit_time:.1f}s", flush=True)

    # Quick diagnostics
    try:
        diag = fit.diagnose()
        print(diag)
    except Exception as e:
        print(f"[diag]  diagnose failed: {e}", flush=True)

    print(f"[save]  CSVs -> {out_dir}", flush=True)


if __name__ == "__main__":
    main()
