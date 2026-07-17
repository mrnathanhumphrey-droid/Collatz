"""
result_bohr_probe_strat.py — bracket-stratified chi² test at k=50.

For each integer-value bracket A..E, build a separate joint histogram on
(Z/32)* x (Z/81)* and run an independent chi² test of independence.

Question: once you strip out bracket-A trajectories (small-attractor funnel,
v ≤ 100), does the v > 100 subpopulation still show non-CRT joint structure?

If H_BE chi²/df ≈ 1 (z near 0), the diagonal Bohr signal is purely descent-
funnel kinematics — no residue-level obstruction. If H_BE has substantial z,
genuine joint structure exists at large v and is the candidate obstruction.
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
TARGET_K = 50
A_, B_ = 5, 4
M2, M3 = 32, 81
ROWS, COLS = 16, 54
DF = (ROWS - 1) * (COLS - 1)  # 795
LOW, HIGH = 500_000_000, 500_000_000_000
INT64_GUARD = 1 << 60
DIAGONAL_RS = {5, 7, 11, 13, 17, 19, 23, 25, 29, 31}
BRACKET_NAMES = ["A:[1,100]", "B:(100,1e4]", "C:(1e4,1e6]",
                 "D:(1e6,1e9]", "E:>1e9"]

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_strat.csv"
TOP_PATH = OUTDIR / "result_bohr_probe_strat_top.csv"
MD_PATH = OUTDIR / "result_bohr_probe_strat.md"


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


def chi2_indep(H):
    N_eff = H.sum()
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
    return chi2, cpdf, z, (P, Q, N_eff)


def top_cells(P, Q, N_eff, k_count=10):
    D = P - Q
    r2_list = odd_residues_mod_2a(A_)
    r3_list = coprime_residues_mod_3b(B_)
    flat = np.argsort(np.abs(D).ravel())[::-1][:k_count]
    out = []
    for fi in flat:
        i = fi // COLS
        j = fi % COLS
        z_c = float(D[i, j] * np.sqrt(N_eff / max(Q[i, j], 1e-30)))
        out.append({
            "r2": int(r2_list[i]), "r3": int(r3_list[j]),
            "P": float(P[i, j]), "Q": float(Q[i, j]),
            "D": float(D[i, j]), "z_cell": z_c,
        })
    return out


def bracket_of(v):
    if v <= 100:
        return 0
    if v <= 10_000:
        return 1
    if v <= 1_000_000:
        return 2
    if v <= 1_000_000_000:
        return 3
    return 4


def main():
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] strat probe; N={N:,}; "
          f"k={TARGET_K}", flush=True)

    inv3 = build_inv3()
    rng = np.random.default_rng(seed=SEED)
    x = rng.integers(low=LOW, high=HIGH, size=N, dtype=np.int64)
    n_arr = (2 * x + 1).astype(np.int64)
    alive = np.ones(N, dtype=bool)
    overflow = np.zeros(N, dtype=bool)
    py_vals = {}
    n_promoted = 0

    for depth in range(1, TARGET_K + 1):
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
                for li in np.where(big_mask)[0]:
                    gi = sub[li]
                    overflow[gi] = True
                    py_vals[gi] = int(m[li])
                    n_promoted += 1
            n_arr[sub] = m
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

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

    n_alive = int(alive.sum())
    print(f"[t={time.time()-t0:6.1f}s] depth={TARGET_K} alive={n_alive:,} "
          f"({100*n_alive/N:.2f}%)", flush=True)

    # Bracket-stratified histograms
    H = {b: np.zeros((ROWS, COLS), dtype=np.int64) for b in range(5)}
    int64_alive_idx = np.where(alive & ~overflow)[0]
    if int64_alive_idx.size > 0:
        vals = n_arr[int64_alive_idx]
        # vectorized bracket via np.searchsorted
        bounds = np.array([100, 10_000, 1_000_000, 1_000_000_000],
                          dtype=np.int64)
        b_arr = np.searchsorted(bounds, vals, side="left")
        # b_arr: 0 for v <= 100, 1 for 100<v<=1e4, ..., 4 for v > 1e9
        # NOTE searchsorted(side="left") gives index where v would be inserted
        # to keep sorted; for v == bound, that's the index of the bound.
        # For our buckets, A: v<=100 → b=0; we want searchsorted([100, 1e4, ...], v) but with v=100 → 0 (correct).
        # Actually searchsorted([100], 100, side="left") = 0, good.
        # searchsorted([100], 101) = 1, good.
        r2 = (vals & (M2 - 1)).astype(np.int64)
        r3 = (vals % M3).astype(np.int64)
        valid = (r2 % 2 == 1) & (r3 % 3 != 0)
        rows_idx = (r2 >> 1)
        cols_idx = inv3[r3]
        for b in range(5):
            mask = (b_arr == b) & valid
            if mask.any():
                flat = rows_idx[mask] * COLS + cols_idx[mask]
                bc = np.bincount(flat, minlength=ROWS * COLS)
                H[b] += bc.reshape(ROWS, COLS)

    # py overflow trajectories
    for gi in [g for g in py_vals if alive[g]]:
        v = py_vals[gi]
        b = bracket_of(v)
        r2 = v & (M2 - 1)
        r3 = v % M3
        if r2 % 2 == 1 and r3 % 3 != 0:
            H[b][r2 >> 1, inv3[r3]] += 1

    # Aggregate H_BE = B + C + D + E
    H_BE = sum(H[b] for b in [1, 2, 3, 4])
    H_full = sum(H[b] for b in range(5))

    # Per-bracket chi²
    print(f"[t={time.time()-t0:6.1f}s] per-bracket chi²", flush=True)
    bracket_results = []
    for b in range(5):
        chi2, cpdf, z, PQN = chi2_indep(H[b])
        N_b = int(H[b].sum())
        if PQN is not None:
            P, Q, _ = PQN
            top = top_cells(P, Q, N_b, k_count=10)
            diag_n = sum(1 for c in top
                         if c["r2"] == c["r3"]
                         and c["r2"] in DIAGONAL_RS)
        else:
            top = []
            diag_n = 0
        bracket_results.append({
            "label": BRACKET_NAMES[b],
            "N": N_b,
            "chi2": chi2, "cpdf": cpdf, "z": z,
            "top": top, "diag_in_top10": diag_n,
        })
        print(f"  {BRACKET_NAMES[b]:<14}  N={N_b:>9,}  cpdf={cpdf:>10.3f}  "
              f"z={z:+11.2f}  diag={diag_n}/10", flush=True)

    # Aggregate BE
    chi2_BE, cpdf_BE, z_BE, PQN_BE = chi2_indep(H_BE)
    N_BE = int(H_BE.sum())
    P_BE, Q_BE, _ = PQN_BE
    top_BE = top_cells(P_BE, Q_BE, N_BE, k_count=10)
    diag_BE = sum(1 for c in top_BE
                  if c["r2"] == c["r3"] and c["r2"] in DIAGONAL_RS)
    print(f"  BE = B+C+D+E    N={N_BE:>9,}  cpdf={cpdf_BE:>10.3f}  "
          f"z={z_BE:+11.2f}  diag={diag_BE}/10", flush=True)

    # Aggregate full (sanity check vs prior runs)
    chi2_full, cpdf_full, z_full, PQN_full = chi2_indep(H_full)
    N_full = int(H_full.sum())
    print(f"  full           N={N_full:>9,}  cpdf={cpdf_full:>10.3f}  "
          f"z={z_full:+11.2f}", flush=True)

    # Outputs
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bracket", "N", "chi2", "chi2_per_df", "z_score",
                    "diag_in_top10"])
        for r in bracket_results:
            w.writerow([r["label"], r["N"], f"{r['chi2']:.4f}",
                        f"{r['cpdf']:.6f}", f"{r['z']:.4f}",
                        r["diag_in_top10"]])
        w.writerow(["BE_aggregate", N_BE, f"{chi2_BE:.4f}",
                    f"{cpdf_BE:.6f}", f"{z_BE:.4f}", diag_BE])
        w.writerow(["full", N_full, f"{chi2_full:.4f}",
                    f"{cpdf_full:.6f}", f"{z_full:.4f}", "—"])

    with open(TOP_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bracket", "rank", "r_2", "r_3", "P", "Q",
                    "P_minus_Q", "z_cell"])
        for r in bracket_results:
            for rank, c in enumerate(r["top"][:5], 1):
                w.writerow([r["label"], rank, c["r2"], c["r3"],
                            f"{c['P']:.6e}", f"{c['Q']:.6e}",
                            f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])
        for rank, c in enumerate(top_BE[:10], 1):
            w.writerow(["BE_aggregate", rank, c["r2"], c["r3"],
                        f"{c['P']:.6e}", f"{c['Q']:.6e}",
                        f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])

    # Markdown
    L = []
    L.append("# Bohr Probe — Bracket-Stratified chi² at k=50")
    L.append("")
    L.append(f"N = {N:,}; k = {TARGET_K}; (a, b) = ({A_}, {B_}); d.f. = {DF}; "
             f"N_alive = {n_alive:,} ({100*n_alive/N:.2f}%)")
    L.append("")
    L.append("Each bracket gets its OWN independent chi² test of joint "
             "independence on (Z/32)* × (Z/81)*. Marginals are computed within "
             "each bracket — so the test is whether the residue distribution "
             "factorizes *conditional on* the integer-value range.")
    L.append("")

    L.append("## Per-bracket chi²/df")
    L.append("")
    L.append("| bracket | N | chi²/df | z | diag in top 10 |")
    L.append("|---|---:|---:|---:|---:|")
    for r in bracket_results:
        L.append(f"| {r['label']} | {r['N']:,} | {r['cpdf']:.3f} | "
                 f"{r['z']:+.2f} | {r['diag_in_top10']}/10 |")
    L.append(f"| **BE = B+C+D+E** | **{N_BE:,}** | **{cpdf_BE:.3f}** | "
             f"**{z_BE:+.2f}** | **{diag_BE}/10** |")
    L.append(f"| full (sanity) | {N_full:,} | {cpdf_full:.3f} | "
             f"{z_full:+.2f} | — |")
    L.append("")

    L.append("## Top 5 cells per bracket")
    for r in bracket_results:
        L.append("")
        L.append(f"### {r['label']}  (N={r['N']:,})")
        L.append("")
        L.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell |")
        L.append("|---:|---:|---:|---:|---:|---:|---:|")
        for rank, c in enumerate(r["top"][:5], 1):
            L.append(f"| {rank} | {c['r2']} | {c['r3']} | {c['P']:.4e} | "
                     f"{c['Q']:.4e} | {c['D']:+.3e} | {c['z_cell']:+.2f} |")

    L.append("")
    L.append("## BE aggregate top 10 cells")
    L.append("")
    L.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell | diag? |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|:---:|")
    for rank, c in enumerate(top_BE[:10], 1):
        mark = ("✓" if (c["r2"] == c["r3"]
                        and c["r2"] in DIAGONAL_RS) else "")
        L.append(f"| {rank} | {c['r2']} | {c['r3']} | {c['P']:.4e} | "
                 f"{c['Q']:.4e} | {c['D']:+.3e} | {c['z_cell']:+.2f} | "
                 f"{mark} |")
    L.append("")

    L.append("## Decision")
    L.append("")
    if abs(z_BE) < 5:
        verdict = ("**The v > 100 subpopulation is CRT-product** "
                   f"(z_BE = {z_BE:+.2f}, within noise floor). The diagonal "
                   "Bohr signal at full chi² is **entirely** an artifact of "
                   "the descent funnel. There is no residue-level joint "
                   "obstruction.")
    elif abs(z_BE) < 50:
        verdict = (f"**Substantial residue-level signal in v > 100** "
                   f"(z_BE = {z_BE:+.2f}). Stripping the funnel removes most "
                   "but not all of the joint structure; the residual is the "
                   "candidate obstruction.")
    else:
        verdict = (f"**Strong residue-level signal in v > 100** "
                   f"(z_BE = {z_BE:+.2f}). Joint structure exists "
                   "independently of small-attractor effects; the obstruction "
                   "lives in the large-value regime.")
    L.append(verdict)
    L.append("")

    MD_PATH.write_text("\n".join(L), encoding="utf-8")
    print(f"\n[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
