"""
Cliff-mapping follow-up to result_bohr_probe.py.

Restricts to the (a=5, b=4) cell — i.e. joint distribution on
(Z/32)* x (Z/81)*, 16 x 54 = 864-cell histogram — and sweeps depths
densely in the k = 15..30 range to characterize the chi^2 cliff
observed between k=15 and k=20.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

SEED = 20260504
N_TARGET = 10_000_000
DEPTHS = [15, 16, 17, 18, 19, 20, 22, 25, 30]
A = 5
B = 4
CHUNK_SIZE = 1_000_000

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_cliff.csv"
MD_PATH = OUTDIR / "result_bohr_probe_cliff.md"


def coprime_residues_mod_3b(b):
    M = 3 ** b
    return np.array([r for r in range(M) if r % 3 != 0], dtype=np.int64)


def main():
    t0 = time.time()
    n_total = N_TARGET
    print(f"[t={time.time()-t0:6.1f}s] starting cliff probe; N={n_total:,}; "
          f"depths={DEPTHS}; (a,b)=({A},{B})", flush=True)

    M2 = 1 << A
    M3 = 3 ** B
    rows = 1 << (A - 1)
    cols = 2 * (3 ** (B - 1))
    r2_list = np.array([r for r in range(M2) if r % 2 == 1], dtype=np.int64)
    r3_list = coprime_residues_mod_3b(B)
    inv3 = -np.ones(M3, dtype=np.int64)
    cnt = 0
    for r in range(M3):
        if r % 3 != 0:
            inv3[r] = cnt
            cnt += 1

    H = {k: np.zeros((rows, cols), dtype=np.int64) for k in DEPTHS}
    n_alive = {k: 0 for k in DEPTHS}
    n_processed = 0

    rng = np.random.default_rng(seed=SEED)
    n_chunks = (n_total + CHUNK_SIZE - 1) // CHUNK_SIZE
    depth_set = set(DEPTHS)
    depth_max = max(DEPTHS)

    for ci in range(n_chunks):
        c_lo = ci * CHUNK_SIZE
        c_hi = min(c_lo + CHUNK_SIZE, n_total)
        sz = c_hi - c_lo

        x = rng.integers(low=0, high=500_000_000_000, size=sz, dtype=np.int64)
        n64 = (2 * x + 1).astype(np.int64)
        alive = np.ones(sz, dtype=bool)

        for depth in range(1, depth_max + 1):
            sub = np.where(alive)[0]
            if sub.size == 0:
                break
            ns = n64[sub]
            m = 3 * ns + 1
            while True:
                even = (m & 1) == 0
                if not even.any():
                    break
                m[even] >>= 1
            n64[sub] = m
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

            if depth in depth_set:
                alive_idx = np.where(alive)[0]
                vals = n64[alive_idx]
                n_alive[depth] += alive_idx.size
                r2 = (vals & (M2 - 1)).astype(np.int64)
                row = (r2 >> 1)
                r3 = (vals % M3).astype(np.int64)
                col = inv3[r3]
                flat = row * cols + col
                bc = np.bincount(flat, minlength=rows * cols)
                H[depth] += bc.reshape(rows, cols)

        n_processed += sz
        elapsed = time.time() - t0
        rate = n_processed / max(elapsed, 1e-9)
        print(f"[t={elapsed:6.1f}s] chunk {ci+1}/{n_chunks} done; "
              f"rate={rate:,.0f}/s; alive(d{depth_max})={n_alive[depth_max]:,}",
              flush=True)

    print(f"[t={time.time()-t0:6.1f}s] computing chi^2 by depth", flush=True)

    df = (rows - 1) * (cols - 1)
    rows_data = []
    for k in DEPTHS:
        N_eff = n_alive[k]
        Hk = H[k]
        P = Hk.astype(np.float64) / N_eff
        Mr = P.sum(axis=1, keepdims=True)
        Mc = P.sum(axis=0, keepdims=True)
        Q = Mr * Mc
        with np.errstate(divide="ignore", invalid="ignore"):
            contrib = np.where(Q > 0, N_eff * (P - Q) ** 2 / Q, 0.0)
        chi2 = float(contrib.sum())
        chi2_per_df = chi2 / df
        z = (chi2 - df) / np.sqrt(2.0 * df)
        rows_data.append((k, df, N_eff, chi2, chi2_per_df, z))
        print(f"  k={k:>2}  N_eff={N_eff:>10,}  chi2={chi2:>9.2f}  "
              f"chi2/df={chi2_per_df:.4f}  z={z:+.2f}", flush=True)

    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "df", "N_eff", "chi2", "chi2_per_df", "z_score"])
        for r in rows_data:
            w.writerow(r)

    lines = []
    lines.append("# Bohr Probe — Cliff Map at (a=5, b=4)")
    lines.append("")
    lines.append("Follow-up to `result_bohr_probe.md`. The k=15→k=20 jump from "
                 "z=+1.10 to z=+16.50 was a hard cliff. This run scans "
                 "k ∈ {15..20, 22, 25, 30} densely at the same (a=5, b=4) cell "
                 f"with N = {n_total:,} starts.")
    lines.append("")
    lines.append("## chi² trajectory")
    lines.append("")
    lines.append("| k | N_eff | chi² | chi²/df | z |")
    lines.append("|---:|---:|---:|---:|---:|")
    for r in rows_data:
        k, df_, N_eff, chi2, c_per_df, z = r
        lines.append(f"| {k} | {N_eff:,} | {chi2:.2f} | {c_per_df:.4f} | "
                     f"{z:+.2f} |")
    lines.append("")
    zs = [r[5] for r in rows_data]
    z_max = max(zs)
    z_min = min(zs)
    lines.append(f"**Range:** z_min = {z_min:+.2f}, z_max = {z_max:+.2f}")
    lines.append("")
    MD_PATH.write_text("\n".join(lines), encoding="utf-8")

    print(f"[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
