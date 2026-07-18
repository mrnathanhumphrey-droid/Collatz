"""
probe_preimage_kalg.py — K_alg variant of preimage structure probe.

Complements probe_preimage_structure.py (which used the brief's M=2^20 integer
lift method, K_emp). This version uses the framework's standard truncated-Geom
chain K_alg matching all prior probes (self-similarity, offset sweep, Atkinson,
R-forcing): v ∈ {1..M_k} with weights 2^{-v}/Z.

K_emp and K_alg are structurally DIFFERENT chains (100% Frobenius diff at k=5):
- K_emp includes v=0 from even integer lifts (each x has its own integer-lift
  identity at v=0, x → (3x+1) mod N as part of the natural-density measure).
- K_alg restricts to v ≥ 1 (one application of /2 forced); this is the "Syracuse
  step on odd integers" chain that drives all framework results (S_k → 7/15 etc.)

Two analyses:
  PART A: full v_eff = M = ord(2 mod 3^k) — proper algebraic chain
  PART B: truncated v_eff = min(M, 60) — matches probe convention

Per column y: |Preimage|_struct (edge count, no underflow), |Preimage|_weighted
(K > 0), max weight, column entropy, column sum.
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

OUTDIR = Path(r"C:\Collatz\probe_preimage_structure")
OUTDIR.mkdir(exist_ok=True)
sys.stdout.reconfigure(encoding="utf-8")


def order_of_2_mod(N):
    o, p = 1, 2 % N
    while p != 1:
        p = (p * 2) % N
        o += 1
    return o


def build_K_alg(k, v_eff_cap=None):
    """K_alg dense build. Returns (K, edge_count, n, v_eff, M, weights)."""
    N = 3 ** k
    M = order_of_2_mod(N)
    v_eff = M if v_eff_cap is None else min(M, v_eff_cap)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(v_eff, dtype=np.int64)
    p = inv2
    for v in range(v_eff):
        powers_inv2[v] = p
        p = (p * inv2) % N
    coprime_mask = np.ones(N, dtype=bool)
    coprime_mask[::3] = False
    coprime_idx_in_N = np.where(coprime_mask)[0]
    n = len(coprime_idx_in_N)
    state_idx = -np.ones(N, dtype=np.int64)
    state_idx[coprime_idx_in_N] = np.arange(n)
    weights = np.array([2.0 ** -(v + 1) for v in range(v_eff)], dtype=np.float64)
    Z = weights.sum()
    if Z > 0:
        weights /= Z
    base = (3 * coprime_idx_in_N + 1) % N
    K = np.zeros((n, n), dtype=np.float64)
    edge_count = np.zeros((n, n), dtype=np.int32)
    row_idx = np.arange(n)
    for v in range(v_eff):
        targets = (base * powers_inv2[v]) % N
        target_states = state_idx[targets]
        # Within fixed v, target_states is a permutation of {0..n-1} (since
        # base is bijective and 2^{-v} is unit), so fancy indexing is safe.
        K[row_idx, target_states] += weights[v]
        edge_count[row_idx, target_states] += 1
    return K, edge_count, n, v_eff, M, weights


def column_stats(K, edge_count):
    n = K.shape[1]
    npre_struct = (edge_count > 0).sum(axis=0)
    npre_weighted = (K > 0).sum(axis=0)
    max_w = K.max(axis=0)
    col_sum = K.sum(axis=0)
    entropy = np.zeros(n, dtype=np.float64)
    for j in range(n):
        col = K[:, j]
        s = col_sum[j]
        if s > 0:
            p = col[col > 0] / s
            entropy[j] = -np.sum(p * np.log(p))
    return npre_struct, npre_weighted, max_w, entropy, col_sum


def linear_fit(x, y):
    A = np.column_stack([np.ones_like(x), x])
    c, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ c
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(c[0]), float(c[1]), float(ss_res), r2


def main():
    print("=" * 78)
    print("K_ALG PREIMAGE PROBE — Tao truncated-Geom chain on (Z/3^k)*, k=5,6,7")
    print("=" * 78)
    print()

    # PART A: full v_eff = M
    print("=" * 78)
    print("PART A: full v_eff = M = ord(2 mod 3^k) (proper Tao chain)")
    print("=" * 78)
    print()
    full = {}
    for k in [5, 6, 7]:
        print(f"--- k = {k} ---")
        t0 = time.time()
        K, ec, n, v_eff, M, weights = build_K_alg(k, v_eff_cap=None)
        t_build = time.time() - t0
        n_underflow = int((weights == 0).sum())
        print(f"  n = {n}, M = {M}, v_eff (math) = {v_eff}, build {t_build:.2f}s")
        if n_underflow > 0:
            print(f"  weights underflow to 0 for {n_underflow} v values "
                  f"(v >= {v_eff - n_underflow}); 2^-v below float64 min "
                  f"~ 5e-324)")
        nz_struct = int((ec > 0).sum())
        nz_w = int((K > 0).sum())
        print(f"  Edges (struct, ec>0): {nz_struct}; expected n·M = {n*M}")
        print(f"  Edges (weighted, K>0): {nz_w}; density {nz_w/(n*n):.4f}")
        nps, npw, mxw, ent, cs = column_stats(K, ec)
        print(f"  |Preimage|_struct: mean={nps.mean():.2f}, "
              f"median={int(np.median(nps))}, max={int(nps.max())}, "
              f"min={int(nps.min())}")
        print(f"  |Preimage|_weighted: mean={npw.mean():.2f}, "
              f"max={int(npw.max())}, min={int(npw.min())}")
        print(f"  Column sum: mean={cs.mean():.6f}, std={cs.std():.4e}, "
              f"min={cs.min():.4e}, max={cs.max():.4e}")
        print(f"  Max weight per column: median={float(np.median(mxw)):.4e}, "
              f"max={mxw.max():.4e}")
        # entropy/log(npw) — concentration ratio
        with np.errstate(divide='ignore', invalid='ignore'):
            ent_ratio = np.where(npw > 1, ent / np.log(npw), 0)
        print(f"  Column entropy: mean={ent.mean():.4f}, "
              f"entropy/log(|P|_w) median = {float(np.median(ent_ratio)):.4f}")
        full[k] = dict(
            n=n, M=M, v_eff=v_eff, n_underflow=n_underflow,
            mean_nps=float(nps.mean()), max_nps=int(nps.max()),
            min_nps=int(nps.min()),
            mean_npw=float(npw.mean()), max_npw=int(npw.max()),
            min_npw=int(npw.min()),
            mean_entropy=float(ent.mean()),
            median_ent_ratio=float(np.median(ent_ratio)),
            mean_col_sum=float(cs.mean()), std_col_sum=float(cs.std()),
            min_col_sum=float(cs.min()), max_col_sum=float(cs.max()),
            median_max_w=float(np.median(mxw)),
        )
        out_csv = OUTDIR / f"preimage_map_k{k}_kalg_full.csv"
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["y_idx", "n_pre_struct", "n_pre_weighted",
                        "max_weight", "col_entropy", "col_sum"])
            for j in range(n):
                w.writerow([j, int(nps[j]), int(npw[j]),
                            f"{float(mxw[j]):.6e}",
                            f"{float(ent[j]):.6f}",
                            f"{float(cs[j]):.6e}"])
        print(f"  saved {out_csv}")
        print()

    # PART B: truncated v_eff = 60
    print("=" * 78)
    print("PART B: truncated v_eff = min(M, 60)")
    print("=" * 78)
    print()
    trunc = {}
    for k in [5, 6, 7]:
        print(f"--- k = {k} ---")
        t0 = time.time()
        K, ec, n, v_eff, M, weights = build_K_alg(k, v_eff_cap=60)
        t_build = time.time() - t0
        print(f"  n = {n}, M = {M}, v_eff = {v_eff}, build {t_build:.2f}s")
        nz_w = int((K > 0).sum())
        print(f"  Edges (weighted, K>0): {nz_w}; density {nz_w/(n*n):.4f}")
        nps, npw, mxw, ent, cs = column_stats(K, ec)
        print(f"  |Preimage|_struct: mean={nps.mean():.2f}, "
              f"max={int(nps.max())}, min={int(nps.min())}")
        print(f"  |Preimage|_weighted: mean={npw.mean():.2f}, "
              f"max={int(npw.max())}, min={int(npw.min())}")
        print(f"  Column sum: mean={cs.mean():.6f}, std={cs.std():.4e}")
        print(f"  Max weight per column: median={float(np.median(mxw)):.4e}")
        with np.errstate(divide='ignore', invalid='ignore'):
            ent_ratio = np.where(npw > 1, ent / np.log(npw), 0)
        print(f"  Column entropy: mean={ent.mean():.4f}, "
              f"entropy/log(|P|_w) median = {float(np.median(ent_ratio)):.4f}")
        trunc[k] = dict(
            n=n, M=M, v_eff=v_eff,
            mean_nps=float(nps.mean()), max_nps=int(nps.max()),
            min_nps=int(nps.min()),
            mean_npw=float(npw.mean()), max_npw=int(npw.max()),
            mean_entropy=float(ent.mean()),
            median_ent_ratio=float(np.median(ent_ratio)),
            mean_col_sum=float(cs.mean()),
            median_max_w=float(np.median(mxw)),
        )
        out_csv = OUTDIR / f"preimage_map_k{k}_kalg_trunc60.csv"
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["y_idx", "n_pre_struct", "n_pre_weighted",
                        "max_weight", "col_entropy", "col_sum"])
            for j in range(n):
                w.writerow([j, int(nps[j]), int(npw[j]),
                            f"{float(mxw[j]):.6e}",
                            f"{float(ent[j]):.6f}",
                            f"{float(cs[j]):.6e}"])
        print(f"  saved {out_csv}")
        print()

    # PART C: scaling fits
    print("=" * 78)
    print("PART C: scaling analysis")
    print("=" * 78)
    print()
    ks = np.array([5, 6, 7], dtype=np.float64)

    nps_full = np.array([full[k]["mean_nps"] for k in [5, 6, 7]])
    npw_full = np.array([full[k]["mean_npw"] for k in [5, 6, 7]])
    nps_trunc = np.array([trunc[k]["mean_nps"] for k in [5, 6, 7]])
    npw_trunc = np.array([trunc[k]["mean_npw"] for k in [5, 6, 7]])
    Ms = np.array([full[k]["M"] for k in [5, 6, 7]], dtype=np.float64)
    ns = np.array([full[k]["n"] for k in [5, 6, 7]], dtype=np.float64)

    print(f"  Full v_eff=M:")
    print(f"    n values:               {ns}")
    print(f"    M values:               {Ms}")
    print(f"    mean |Pre|_struct:      {nps_full}")
    print(f"    mean |Pre|_weighted:    {npw_full}")
    print()

    a, b, ss, r2 = linear_fit(ks, nps_full)
    print(f"  Linear |Pre|_struct = {a:.4f} + {b:.4f}·k, R² = {r2:.4f}")
    a, b, ss, r2 = linear_fit(ks, np.log(nps_full))
    print(f"  log|Pre|_struct = {a:.4f} + {b:.4f}·k → exp(b) = {math.exp(b):.4f}")
    print(f"    Predicted (|Pre|=M=2·3^(k-1)): rate=3 per step, "
          f"intercept=log(2/3)={math.log(2/3):.4f}")
    a, b, ss, r2 = linear_fit(ks, np.log(npw_full))
    print(f"  log|Pre|_weighted = {a:.4f} + {b:.4f}·k → exp(b) = {math.exp(b):.4f}")
    print()
    print(f"  Truncated v_eff=min(M,60):")
    print(f"    mean |Pre|_struct:      {nps_trunc}")
    print(f"    mean |Pre|_weighted:    {npw_trunc}")
    a, b, ss, r2 = linear_fit(ks, nps_trunc)
    print(f"    Linear |Pre|_struct = {a:.4f} + {b:.4f}·k, R² = {r2:.4f}")
    a, b, ss, r2 = linear_fit(ks, npw_trunc)
    print(f"    Linear |Pre|_weighted = {a:.4f} + {b:.4f}·k, R² = {r2:.4f}")
    print()

    # Save scaling CSV
    summary_csv = OUTDIR / "scaling_analysis_kalg.csv"
    with open(summary_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "n", "M", "version", "v_eff", "mean_npre_struct",
                    "max_npre_struct", "min_npre_struct", "mean_npre_weighted",
                    "max_npre_weighted", "mean_entropy", "median_ent_over_log",
                    "mean_col_sum", "median_max_weight"])
        for k in [5, 6, 7]:
            r = full[k]
            w.writerow([k, r["n"], r["M"], "full", r["v_eff"],
                        f"{r['mean_nps']:.4f}", r["max_nps"], r["min_nps"],
                        f"{r['mean_npw']:.4f}", r["max_npw"],
                        f"{r['mean_entropy']:.6f}",
                        f"{r['median_ent_ratio']:.6f}",
                        f"{r['mean_col_sum']:.6f}",
                        f"{r['median_max_w']:.6e}"])
            r = trunc[k]
            w.writerow([k, r["n"], r["M"], "trunc60", r["v_eff"],
                        f"{r['mean_nps']:.4f}", r["max_nps"], r["min_nps"],
                        f"{r['mean_npw']:.4f}", r["max_npw"],
                        f"{r['mean_entropy']:.6f}",
                        f"{r['median_ent_ratio']:.6f}",
                        f"{r['mean_col_sum']:.6f}",
                        f"{r['median_max_w']:.6e}"])
    print(f"  saved {summary_csv}")
    print()

    # PART D: row-collision analysis (3:1 collapse)
    print("=" * 78)
    print("PART D: row-collision structure (3:1 collapse from 3x+1 mod 3^(k-1))")
    print("=" * 78)
    print()
    for k in [5, 6, 7]:
        K, ec, n, v_eff, M, weights = build_K_alg(k, v_eff_cap=None)
        # Group rows by (3x+1) mod 3^(k-1) — should give n/3 groups of 3 identical rows
        N = 3 ** k
        coprime = np.array([r for r in range(N) if r % 3 != 0], dtype=np.int64)
        base_high = (3 * coprime + 1) % N  # full 3x+1 mod 3^k
        base_low = base_high % (3 ** (k - 1))  # mod 3^(k-1)
        unique_classes = np.unique(base_low)
        n_classes = len(unique_classes)
        # Verify rows in same class are identical (within float tolerance)
        sample_class = unique_classes[0]
        sample_rows = np.where(base_low == sample_class)[0]
        if len(sample_rows) >= 2:
            row_diff = np.linalg.norm(K[sample_rows[0]] - K[sample_rows[1]])
            print(f"  k={k}: {n_classes} row-equivalence-classes "
                  f"(predicted n/3 = {n//3}); "
                  f"sample diff between rows in same class = {row_diff:.2e}")
        # rank of K (numerical)
        s = np.linalg.svd(K, compute_uv=False)
        rank_eps = int((s > 1e-12 * s[0]).sum())
        print(f"      numerical rank (rel tol 1e-12): {rank_eps} "
              f"(predicted n/3 = {n//3})")
    print()

    print(f"All outputs in {OUTDIR}")


if __name__ == "__main__":
    main()
