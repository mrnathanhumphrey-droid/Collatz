"""
Stage 3 prep: load Stage 1 main parquet, subset to odd n > 1, compute residue
class index mod 2^k (default k=10 = 512 odd classes), emit a Stan-ready
parquet with columns (n, sigma, log_n, class_idx).

Reasoning for choices:
  - Odd-only subset removes the deterministic shift sigma(2n) = sigma(n)+1
    that would otherwise correlate residuals with nu_2(n). Cleaner signal.
  - mod 2^k with odd residues: residues {1, 3, ..., 2^k-1} -> {0, 1, ..., 2^(k-1)-1}.
    For k=10, that's 512 classes. We keep class label as class_idx in [0, K-1]
    and a separate residue column for interpretability.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=1 << 25,
                    help="N from Stage 1 (default 2^25).")
    ap.add_argument("--k", type=int, default=10,
                    help="Residue power: group by n mod 2^k. Default 10.")
    ap.add_argument("--data", type=str, default=None,
                    help="Data dir. Default <script_dir>/data.")
    ap.add_argument("--out", type=str, default=None,
                    help="Output parquet path.")
    ap.add_argument("--per-class-cap", type=int, default=None,
                    help="If set, randomly subsample at most this many rows "
                         "per class (preserving the upper tail via stratified "
                         "sampling). Useful to keep Stan tractable.")
    ap.add_argument("--tail-frac", type=float, default=0.05,
                    help="Per-class tail fraction always kept verbatim when "
                         "subsampling. Default 0.05 = top 5%.")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    data_dir = Path(args.data) if args.data else here / "data"
    src = data_dir / f"main_N{args.N}.parquet"
    print(f"[load]  {src}", flush=True)
    df = pl.read_parquet(src)
    print(f"        full: {df.height:,} rows", flush=True)

    df = df.filter((pl.col("n") % 2 == 1) & (pl.col("n") > 1))
    print(f"        odd>1: {df.height:,} rows", flush=True)

    M = 1 << args.k
    df = df.with_columns([
        (pl.col("n") % M).cast(pl.Int32).alias(f"res_mod_{M}"),
        pl.col("n").cast(pl.Float64).log().alias("log_n"),
    ])
    df = df.with_columns(
        ((pl.col(f"res_mod_{M}") - 1) // 2).cast(pl.Int32).alias("class_idx")
    )

    K = (M // 2)
    classes = df["class_idx"].unique().sort().to_list()
    print(f"[class] {len(classes)} unique classes (expected {K})", flush=True)
    assert len(classes) == K, "missing odd residue classes — data range too small"

    if args.per_class_cap is not None:
        cap = args.per_class_cap
        tf = args.tail_frac
        # heuristic residual = sigma - (2/(log4-log3)) * log(n)
        log_factor = 2.0 / (np.log(4.0) - np.log(3.0))
        df = df.with_columns(
            (pl.col("sigma").cast(pl.Float64) - log_factor * pl.col("log_n"))
            .alias("resid")
        )

        # Per-class: keep top tf-fraction by resid + uniform-random sample
        # filling up to cap rows total.
        keep_parts = []
        rng = np.random.default_rng(42)
        for k in classes:
            sub = df.filter(pl.col("class_idx") == k)
            n_sub = sub.height
            n_tail = max(1, int(round(tf * n_sub)))
            n_uniform = max(0, cap - n_tail)
            sub_sorted = sub.sort("resid", descending=True)
            tail = sub_sorted.head(n_tail)
            rest = sub_sorted.tail(n_sub - n_tail)
            if rest.height > n_uniform:
                # random sample without replacement
                idx = rng.choice(rest.height, size=n_uniform, replace=False)
                rest = rest[idx.tolist()]
            keep_parts.append(pl.concat([tail, rest]))
        df = pl.concat(keep_parts).drop("resid")
        print(f"[cap]   stratified-sampled: {df.height:,} rows "
              f"(cap {cap}/class, top {tf:.0%} retained)", flush=True)

    out_path = (Path(args.out) if args.out
                else data_dir / f"stage3_input_N{args.N}_k{args.k}.parquet")
    df.select(["n", "sigma", "log_n", "class_idx", f"res_mod_{M}"]) \
      .write_parquet(out_path, compression="zstd", compression_level=3)

    sz_mb = out_path.stat().st_size / 1e6
    print(f"[save]  {df.height:,} rows -> {out_path.name} ({sz_mb:.1f} MB)",
          flush=True)


if __name__ == "__main__":
    main()
