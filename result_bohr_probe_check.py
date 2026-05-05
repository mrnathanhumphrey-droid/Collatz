"""
result_bohr_probe_check.py — three follow-up checks on the Bohr signal.

  CHECK 1 — high-N restriction: starts in [10^9, 10^12] (no small n).
            Compare chi²/df and top 10 cells at k=30 to the unrestricted run.
  CHECK 2 — asymptotic saturation: extend to k = 60, 80, 100. Top 5 per depth.
  CHECK 3 — diagonal mechanism: at k=30, for cells (r, r) with r ∈
            {5, 11, 13, 17, 19, 29}, dump the value distribution of trajectories
            landing there. Quantify how much of the over-mass is from "exact
            small attractor" landings vs larger integers with the right residue.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

SEED = 20260504
N = 1_000_000
DEPTHS = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]
A, B = 5, 4
M2, M3 = 1 << A, 3 ** B
ROWS, COLS = 1 << (A - 1), 2 * 3 ** (B - 1)
DF = (ROWS - 1) * (COLS - 1)
DIAGONAL_RS = [5, 11, 13, 17, 19, 29]
TARGET_DEPTH = 30
LOW, HIGH = 500_000_000, 500_000_000_000  # 2x+1 spans [10^9+1, 10^12-1) odd
INT64_GUARD = 1 << 60

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_check.csv"
TOP_PATH = OUTDIR / "result_bohr_probe_check_top.csv"
VAL_PATH = OUTDIR / "result_bohr_probe_check_values.csv"
MD_PATH = OUTDIR / "result_bohr_probe_check.md"


def coprime_residues_mod_3b(b):
    return np.array([r for r in range(3 ** b) if r % 3 != 0], dtype=np.int64)


def odd_residues_mod_2a(a):
    return np.array([r for r in range(1 << a) if r % 2 == 1], dtype=np.int64)


def build_inv3():
    inv = -np.ones(M3, dtype=np.int64)
    cnt = 0
    for r in range(M3):
        if r % 3 != 0:
            inv[r] = cnt
            cnt += 1
    return inv


def histogram(vals, inv3):
    r2 = (vals & (M2 - 1)).astype(np.int64)
    row = (r2 >> 1)
    r3 = (vals % M3).astype(np.int64)
    col = inv3[r3]
    flat = row * COLS + col
    return np.bincount(flat, minlength=ROWS * COLS).reshape(ROWS, COLS)


def chi2_df_from_H(H, N_eff):
    if N_eff == 0:
        return 0.0, float("nan"), float("nan"), None
    P = H.astype(np.float64) / N_eff
    Mr = P.sum(axis=1, keepdims=True)
    Mc = P.sum(axis=0, keepdims=True)
    Q = Mr * Mc
    with np.errstate(divide="ignore", invalid="ignore"):
        contrib = np.where(Q > 0, N_eff * (P - Q) ** 2 / Q, 0.0)
    chi2 = float(contrib.sum())
    cpdf = chi2 / DF
    z = (chi2 - DF) / np.sqrt(2.0 * DF)
    return chi2, cpdf, z, (P, Q)


def top_cells(H, N_eff, k_count=10):
    _, _, _, PQ = chi2_df_from_H(H, N_eff)
    if PQ is None:
        return []
    P, Q = PQ
    D = P - Q
    r2_list = odd_residues_mod_2a(A)
    r3_list = coprime_residues_mod_3b(B)
    flat = np.argsort(np.abs(D).ravel())[::-1][:k_count]
    out = []
    for fi in flat:
        i = fi // COLS
        j = fi % COLS
        out.append({
            "r2": int(r2_list[i]),
            "r3": int(r3_list[j]),
            "P": float(P[i, j]),
            "Q": float(Q[i, j]),
            "D": float(D[i, j]),
            "z_cell": float(D[i, j] * np.sqrt(N_eff / max(Q[i, j], 1e-30))),
        })
    return out


def main():
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] check probe; N={N:,}; "
          f"starts odd in [{2*LOW+1:,}, {2*HIGH-1:,}]; depths={DEPTHS}",
          flush=True)

    inv3 = build_inv3()
    rng = np.random.default_rng(seed=SEED)
    x = rng.integers(low=LOW, high=HIGH, size=N, dtype=np.int64)
    n_arr = (2 * x + 1).astype(np.int64)
    starts = n_arr.copy()
    alive = np.ones(N, dtype=bool)
    overflow = np.zeros(N, dtype=bool)
    py_vals = {}

    H = {k: np.zeros((ROWS, COLS), dtype=np.int64) for k in DEPTHS}
    n_alive = {k: 0 for k in DEPTHS}
    n_promoted_total = 0

    cell_values = {(r, r): [] for r in DIAGONAL_RS}
    cell_starts = {(r, r): [] for r in DIAGONAL_RS}

    depth_set = set(DEPTHS)
    depth_max = max(DEPTHS)

    for depth in range(1, depth_max + 1):
        # int64 vectorized step
        int64_mask = alive & ~overflow
        sub = np.where(int64_mask)[0]
        if sub.size > 0:
            ns = n_arr[sub]
            m = 3 * ns + 1
            while True:
                even = (m & 1) == 0
                if not even.any():
                    break
                m[even] >>= 1
            big_mask = m > INT64_GUARD
            if big_mask.any():
                big_li = np.where(big_mask)[0]
                for li in big_li:
                    gi = sub[li]
                    overflow[gi] = True
                    py_vals[gi] = int(m[li])
                    n_promoted_total += 1
            n_arr[sub] = m
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

        # Python int step on overflow trajectories
        py_idx = [gi for gi in list(py_vals) if alive[gi]]
        for gi in py_idx:
            v = py_vals[gi]
            m = 3 * v + 1
            while not (m & 1):
                m >>= 1
            if m == 1:
                alive[gi] = False
                del py_vals[gi]
                overflow[gi] = False
                n_arr[gi] = 1
            elif m <= INT64_GUARD:
                n_arr[gi] = m
                del py_vals[gi]
                overflow[gi] = False
            else:
                py_vals[gi] = m

        if depth in depth_set:
            int64_alive_idx = np.where(alive & ~overflow)[0]
            if int64_alive_idx.size > 0:
                H[depth] += histogram(n_arr[int64_alive_idx], inv3)
            py_alive_idx = [gi for gi in py_vals if alive[gi]]
            for gi in py_alive_idx:
                v = py_vals[gi]
                r2 = v & (M2 - 1)
                r3 = v % M3
                if r2 % 2 == 1 and r3 % 3 != 0:
                    H[depth][r2 >> 1, inv3[r3]] += 1
            n_alive[depth] = int(alive.sum())

            if depth == TARGET_DEPTH:
                # CHECK 3: gather actual values for target diagonal cells
                # int64-resident
                vals = n_arr[int64_alive_idx]
                r2v = (vals & (M2 - 1)).astype(np.int64)
                r3v = (vals % M3).astype(np.int64)
                for r in DIAGONAL_RS:
                    cell_mask = (r2v == r) & (r3v == r)
                    hits = int64_alive_idx[cell_mask]
                    for gi in hits:
                        cell_values[(r, r)].append(int(n_arr[gi]))
                        cell_starts[(r, r)].append(int(starts[gi]))
                # python ints
                for gi in py_alive_idx:
                    v = py_vals[gi]
                    r2 = v & (M2 - 1)
                    r3 = v % M3
                    if r2 == r3 and r2 in DIAGONAL_RS:
                        cell_values[(r2, r2)].append(v)
                        cell_starts[(r2, r2)].append(int(starts[gi]))

            print(f"  [t={time.time()-t0:6.1f}s] depth={depth:>3}  "
                  f"alive={n_alive[depth]:>9,} ({100*n_alive[depth]/N:6.2f}%)  "
                  f"overflow_total={n_promoted_total}", flush=True)

    print(f"[t={time.time()-t0:6.1f}s] computing chi2 + top cells", flush=True)

    rows_data = []
    top_per_depth = {}
    for k in DEPTHS:
        chi2, cpdf, z, _ = chi2_df_from_H(H[k], n_alive[k])
        rows_data.append({
            "k": k, "n_alive": n_alive[k], "chi2": chi2, "cpdf": cpdf, "z": z,
        })
        n_top = 10 if k in (30, 60, 80, 100) else 5
        top_per_depth[k] = top_cells(H[k], n_alive[k], k_count=n_top)
        print(f"  k={k:>3}  alive={n_alive[k]:>9,}  "
              f"cpdf={cpdf:>10.3f}  z={z:+11.2f}", flush=True)

    # ---- CHECK 3 analysis: bin cell-values by magnitude ----
    print(f"\n[t={time.time()-t0:6.1f}s] CHECK 3 analysis (target k={TARGET_DEPTH})",
          flush=True)
    val_rows = []
    for r in DIAGONAL_RS:
        vals = cell_values[(r, r)]
        starts_r = cell_starts[(r, r)]
        n_total = len(vals)
        if n_total == 0:
            print(f"  cell ({r},{r}): EMPTY", flush=True)
            continue
        n_eq_r = sum(1 for v in vals if v == r)
        n_le_100 = sum(1 for v in vals if v <= 100)
        n_le_1e4 = sum(1 for v in vals if v <= 10_000)
        n_le_1e8 = sum(1 for v in vals if v <= 100_000_000)
        v_min = min(vals)
        v_med = sorted(vals)[len(vals) // 2]
        v_max = max(vals)
        # split by starting integer size: starts > 10^9 by construction, so
        # we report the start range and how many starts ≤ 10^10 vs > 10^10
        # (subdividing the [10^9, 10^12] band)
        n_st_le_1e10 = sum(1 for s in starts_r if s <= 10_000_000_000)
        n_st_gt_1e11 = sum(1 for s in starts_r if s > 100_000_000_000)
        val_rows.append({
            "r": r,
            "n_total": n_total,
            "n_eq_r": n_eq_r,
            "n_le_100": n_le_100,
            "n_le_1e4": n_le_1e4,
            "n_le_1e8": n_le_1e8,
            "v_min": v_min,
            "v_med": v_med,
            "v_max": v_max,
            "n_st_le_1e10": n_st_le_1e10,
            "n_st_gt_1e11": n_st_gt_1e11,
        })
        print(f"  cell ({r},{r}): n={n_total}  v_eq_r={n_eq_r}  "
              f"v_le_100={n_le_100}  v_le_1e4={n_le_1e4}  "
              f"v_med={v_med}  v_max={v_max}", flush=True)

    # ---- write CSVs ----
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "n_alive", "chi2", "chi2_per_df", "z_score"])
        for r in rows_data:
            w.writerow([r["k"], r["n_alive"], f"{r['chi2']:.4f}",
                        f"{r['cpdf']:.6f}", f"{r['z']:.4f}"])

    with open(TOP_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "rank", "r_2", "r_3", "P", "Q", "P_minus_Q", "z_cell"])
        for k in DEPTHS:
            for rank, c in enumerate(top_per_depth[k], 1):
                w.writerow([k, rank, c["r2"], c["r3"],
                            f"{c['P']:.6e}", f"{c['Q']:.6e}",
                            f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])

    with open(VAL_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["r", "n_total", "n_eq_r", "n_le_100", "n_le_1e4",
                    "n_le_1e8", "v_min", "v_median", "v_max",
                    "n_start_le_1e10", "n_start_gt_1e11"])
        for vr in val_rows:
            w.writerow([vr["r"], vr["n_total"], vr["n_eq_r"], vr["n_le_100"],
                        vr["n_le_1e4"], vr["n_le_1e8"], vr["v_min"],
                        vr["v_med"], vr["v_max"],
                        vr["n_st_le_1e10"], vr["n_st_gt_1e11"]])

    # ---- markdown writeup ----
    # Hardcoded comparison: previous unrestricted run's top 10 at k=30
    # (from result_bohr_probe_verify.md, N=10^6 over [1, 10^12])
    prev_top30 = [
        (1, 5, 5, 24.22), (2, 11, 11, 22.12), (3, 17, 17, 16.34),
        (4, 13, 13, 10.00), (5, 19, 19, 15.74), (6, 29, 29, 13.18),
        (7, 15, 47, 8.50), (8, 5, 37, 8.45), (9, 1, 80, 3.82),
        (10, 1, 65, 8.73),
    ]

    lines = []
    lines.append("# Bohr Probe — Verification Checks "
                 "(restriction, saturation, diagonal mechanism)")
    lines.append("")
    lines.append(f"N = {N:,}; (a, b) = ({A}, {B}); d.f. = {DF}; "
                 f"starts odd in [{2*LOW+1:,}, {2*HIGH-1:,}]")
    lines.append("")

    # CHECK 1
    lines.append("## CHECK 1 — High-N restriction (starts ≥ 10^9)")
    lines.append("")
    chi2_30, cpdf_30, z_30, _ = chi2_df_from_H(H[30], n_alive[30])
    lines.append(f"At k = 30: chi²/df = **{cpdf_30:.3f}**, "
                 f"z = **{z_30:+.2f}** (vs unrestricted 6.367, +107.00)")
    lines.append("")
    lines.append("**Top 10 cells at k=30 (high-N restricted) vs unrestricted:**")
    lines.append("")
    lines.append("| rank | this run (r2, r3) | z_cell | unrestricted (r2, r3) | "
                 "z_cell unrestricted |")
    lines.append("|---:|---|---:|---|---:|")
    cur = top_per_depth[30]
    for i in range(10):
        c = cur[i]
        prev = prev_top30[i]
        lines.append(f"| {i+1} | ({c['r2']}, {c['r3']}) | {c['z_cell']:+.2f} | "
                     f"({prev[1]}, {prev[2]}) | {prev[3]:+.2f} |")
    lines.append("")
    diag_in_top10 = sum(1 for c in cur[:10]
                        if c["r2"] == c["r3"] and c["r2"] in DIAGONAL_RS)
    lines.append(f"**Diagonal cells (r=r ∈ {{5,11,13,17,19,29}}) in top 10:** "
                 f"{diag_in_top10}/6")
    lines.append("")

    # CHECK 2
    lines.append("## CHECK 2 — Asymptotic saturation (k = 60, 80, 100)")
    lines.append("")
    lines.append("| k | N_alive | chi²/df | z |")
    lines.append("|---:|---:|---:|---:|")
    for r in rows_data:
        lines.append(f"| {r['k']} | {r['n_alive']:,} | {r['cpdf']:.3f} | "
                     f"{r['z']:+.2f} |")
    lines.append("")
    cpdf50 = next(r["cpdf"] for r in rows_data if r["k"] == 50)
    cpdf100 = next(r["cpdf"] for r in rows_data if r["k"] == 100)
    sat_ratio = cpdf100 / cpdf50 if cpdf50 > 0 else float("nan")
    lines.append(f"**Saturation ratio chi²(100) / chi²(50) = {sat_ratio:.3f}** "
                 "(< 1.3 = saturating; > 2.0 = still growing)")
    lines.append("")
    lines.append("Per-step growth ratios:")
    lines.append("")
    lines.append("| range | per-step r |")
    lines.append("|---|---:|")
    pairs = [(50, 60), (60, 80), (80, 100)]
    for (k_lo, k_hi) in pairs:
        c_lo = next(r["cpdf"] for r in rows_data if r["k"] == k_lo)
        c_hi = next(r["cpdf"] for r in rows_data if r["k"] == k_hi)
        if c_lo > 0:
            r_step = (c_hi / c_lo) ** (1.0 / (k_hi - k_lo))
        else:
            r_step = float("nan")
        lines.append(f"| {k_lo}→{k_hi} | {r_step:.4f} |")
    lines.append("")
    lines.append("**Top 5 cells at k = 60, 80, 100:**")
    for k in [60, 80, 100]:
        lines.append("")
        lines.append(f"### k = {k}")
        lines.append("")
        lines.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell |")
        lines.append("|---:|---:|---:|---:|---:|---:|---:|")
        for rank, c in enumerate(top_per_depth[k][:5], 1):
            lines.append(f"| {rank} | {c['r2']} | {c['r3']} | "
                         f"{c['P']:.4e} | {c['Q']:.4e} | "
                         f"{c['D']:+.3e} | {c['z_cell']:+.2f} |")

    lines.append("")

    # CHECK 3
    lines.append("## CHECK 3 — Diagonal mechanism (cells (r,r) at k=30)")
    lines.append("")
    lines.append("| r | n_total | n_eq_r | n_le_100 | n_le_1e4 | n_le_1e8 | "
                 "v_min | v_median | v_max |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for vr in val_rows:
        lines.append(f"| {vr['r']} | {vr['n_total']} | {vr['n_eq_r']} | "
                     f"{vr['n_le_100']} | {vr['n_le_1e4']} | "
                     f"{vr['n_le_1e8']} | {vr['v_min']} | {vr['v_med']:,} | "
                     f"{vr['v_max']:,} |")
    lines.append("")
    lines.append("**Interpretation:** if `n_eq_r ≈ n_total` for any cell, that "
                 "cell's mass is ~entirely from depth-30 trajectories whose "
                 "VALUE equals r exactly (small-attractor concentration). "
                 "If `v_median ≫ r`, the over-mass is from large-value "
                 "trajectories landing on residue r (genuine joint-measure "
                 "structure).")
    lines.append("")

    MD_PATH.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
