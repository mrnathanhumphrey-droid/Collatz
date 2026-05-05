"""
result_bohr_probe_asymptote.py — fine-grained chi²/df at k = 40..70 (steps 5)
on the high-N restricted ensemble. Fits the saturating model
   chi²/df(k) = A − B · exp(−C·k)
to extract the bounded asymptote A. Plus diagonal-extension probe at k=50.

d.f. note: standard chi² of independence on 16×54 is (16−1)(54−1) = 795.
The user-supplied 8·54−1 = 431 appears to be a typo (a=4 not a=5). This run
uses 795 throughout for comparability with prior probes.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

sys.stdout.reconfigure(encoding="utf-8")

SEED = 20260504
N = 1_000_000
DEPTHS = [40, 45, 50, 55, 60, 65, 70]
A, B = 5, 4
M2, M3 = 32, 81
ROWS, COLS = 16, 54
DF = (ROWS - 1) * (COLS - 1)  # 795
LOW, HIGH = 500_000_000, 500_000_000_000
INT64_GUARD = 1 << 60
KNOWN_DIAGONALS = {5, 11, 13, 17, 19, 29}
EXT_CANDIDATES = [7, 23, 25, 31, 37]
EXT_K = 50

OUTDIR = Path(r"C:\Collatz")
CSV_PATH = OUTDIR / "result_bohr_probe_asymptote.csv"
TOP_PATH = OUTDIR / "result_bohr_probe_asymptote_top.csv"
EXT_PATH = OUTDIR / "result_bohr_probe_asymptote_ext.csv"
MD_PATH = OUTDIR / "result_bohr_probe_asymptote.md"


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


def chi2_stuff(H, N_eff):
    P = H.astype(np.float64) / N_eff
    Mr = P.sum(axis=1, keepdims=True)
    Mc = P.sum(axis=0, keepdims=True)
    Q = Mr * Mc
    with np.errstate(divide="ignore", invalid="ignore"):
        contrib = np.where(Q > 0, N_eff * (P - Q) ** 2 / Q, 0.0)
    chi2 = float(contrib.sum())
    return chi2, P, Q


def top_cells_with_z(P, Q, N_eff, k_count=10):
    D = P - Q
    r2_list = odd_residues_mod_2a(A)
    r3_list = coprime_residues_mod_3b(B)
    flat = np.argsort(np.abs(D).ravel())[::-1][:k_count]
    out = []
    for fi in flat:
        i = fi // COLS
        j = fi % COLS
        r2 = int(r2_list[i])
        r3 = int(r3_list[j])
        z_cell = float(D[i, j] * np.sqrt(N_eff / max(Q[i, j], 1e-30)))
        out.append({"r2": r2, "r3": r3, "P": float(P[i, j]),
                    "Q": float(Q[i, j]), "D": float(D[i, j]),
                    "z_cell": z_cell})
    return out


def main():
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] asymptote probe; N={N:,}; depths={DEPTHS}",
          flush=True)

    inv3 = build_inv3()
    rng = np.random.default_rng(seed=SEED)
    x = rng.integers(low=LOW, high=HIGH, size=N, dtype=np.int64)
    n_arr = (2 * x + 1).astype(np.int64)
    alive = np.ones(N, dtype=bool)
    overflow = np.zeros(N, dtype=bool)
    py_vals = {}

    H = {k: np.zeros((ROWS, COLS), dtype=np.int64) for k in DEPTHS}
    n_alive = {k: 0 for k in DEPTHS}
    n_promoted = 0
    depth_set = set(DEPTHS)
    depth_max = max(DEPTHS)

    for depth in range(1, depth_max + 1):
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

        if depth in depth_set:
            int64_alive_idx = np.where(alive & ~overflow)[0]
            if int64_alive_idx.size > 0:
                H[depth] += histogram(n_arr[int64_alive_idx], inv3)
            for gi in [g for g in py_vals if alive[g]]:
                v = py_vals[gi]
                r2 = v & (M2 - 1)
                r3 = v % M3
                if r2 % 2 == 1 and r3 % 3 != 0:
                    H[depth][r2 >> 1, inv3[r3]] += 1
            n_alive[depth] = int(alive.sum())
            print(f"  [t={time.time()-t0:6.1f}s] depth={depth:>3}  "
                  f"alive={n_alive[depth]:>9,} ({100*n_alive[depth]/N:6.2f}%)  "
                  f"overflow={n_promoted}", flush=True)

    print(f"[t={time.time()-t0:6.1f}s] computing chi2 + tops + fit", flush=True)

    rows_data = []
    top_per_depth = {}
    P_per_depth = {}
    Q_per_depth = {}
    for k in DEPTHS:
        chi2, P, Q = chi2_stuff(H[k], n_alive[k])
        cpdf = chi2 / DF
        z = (chi2 - DF) / np.sqrt(2.0 * DF)
        top = top_cells_with_z(P, Q, n_alive[k], k_count=10)
        diag_count = sum(1 for c in top
                         if c["r2"] == c["r3"] and c["r2"] in KNOWN_DIAGONALS)
        rows_data.append({
            "k": k, "n_alive": n_alive[k], "chi2": chi2, "cpdf": cpdf,
            "z": z, "diag": diag_count,
        })
        top_per_depth[k] = top
        P_per_depth[k] = P
        Q_per_depth[k] = Q
        print(f"  k={k:>2}  alive={n_alive[k]:>8,}  cpdf={cpdf:>9.3f}  "
              f"z={z:+10.2f}  diag={diag_count}/6", flush=True)

    # Asymptote fit
    ks = np.array([r["k"] for r in rows_data], dtype=float)
    cpdfs = np.array([r["cpdf"] for r in rows_data], dtype=float)

    A_fit = B_fit = C_fit = float("nan")
    fit_method = "none"
    if HAS_SCIPY:
        def model(k, A_, B_, C_):
            return A_ - B_ * np.exp(-C_ * k)
        for A0_factor in [1.5, 2.0, 3.0, 5.0]:
            A0 = float(cpdfs.max()) * A0_factor
            B0 = A0 - cpdfs[0]
            C0 = 0.05
            try:
                popt, _ = curve_fit(model, ks, cpdfs, p0=[A0, B0, C0],
                                    maxfev=50000)
                A_fit, B_fit, C_fit = popt
                fit_method = f"scipy curve_fit (A0_factor={A0_factor})"
                break
            except Exception as e:
                print(f"[fit] curve_fit failed at A0_factor={A0_factor}: {e}",
                      flush=True)
    if not (A_fit == A_fit):  # NaN check
        # Grid search fallback on A
        best = None
        for A_try in np.linspace(cpdfs.max() * 1.05, cpdfs.max() * 10, 500):
            if A_try <= cpdfs.max():
                continue
            y = np.log(A_try - cpdfs)
            slope, intercept = np.polyfit(ks, y, 1)
            B_try = float(np.exp(intercept))
            C_try = float(-slope)
            if C_try <= 0:
                continue
            pred = A_try - B_try * np.exp(-C_try * ks)
            ss_res = float(np.sum((cpdfs - pred) ** 2))
            if best is None or ss_res < best[0]:
                best = (ss_res, A_try, B_try, C_try)
        if best is not None:
            _, A_fit, B_fit, C_fit = best
            fit_method = "grid-search fallback"

    print(f"[fit:{fit_method}] A = {A_fit:.2f}, B = {B_fit:.2f}, "
          f"C = {C_fit:.5f}", flush=True)

    # EXT diagonal extension at k=50
    print(f"[t={time.time()-t0:6.1f}s] EXT check at k={EXT_K}", flush=True)
    P_ext = P_per_depth[EXT_K]
    Q_ext = Q_per_depth[EXT_K]
    N_ext = n_alive[EXT_K]
    ext_results = []
    for r in EXT_CANDIDATES:
        r2 = r % M2
        r3 = r % M3
        if r2 % 2 == 0 or r3 % 3 == 0:
            print(f"  r={r}: invalid (not in coprime classes)", flush=True)
            continue
        i = r2 >> 1
        j = inv3[r3]
        P_ij = float(P_ext[i, j])
        Q_ij = float(Q_ext[i, j])
        D_ij = P_ij - Q_ij
        z_c = D_ij * np.sqrt(N_ext / max(Q_ij, 1e-30))
        ext_results.append({
            "r": r, "r2": int(r2), "r3": int(r3),
            "is_diag": (r2 == r3),
            "P": P_ij, "Q": Q_ij, "D": D_ij, "z_cell": float(z_c),
        })
        print(f"  r={r:>2} → cell ({r2},{r3})  z_cell={z_c:+.2f}", flush=True)

    # Outputs
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "n_alive", "chi2", "chi2_per_df", "z_score",
                    "diag_in_top10"])
        for r in rows_data:
            w.writerow([r["k"], r["n_alive"], f"{r['chi2']:.4f}",
                        f"{r['cpdf']:.6f}", f"{r['z']:.4f}", r["diag"]])

    with open(TOP_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "rank", "r_2", "r_3", "P", "Q", "P_minus_Q",
                    "z_cell"])
        for k in DEPTHS:
            for rank, c in enumerate(top_per_depth[k][:5], 1):
                w.writerow([k, rank, c["r2"], c["r3"],
                            f"{c['P']:.6e}", f"{c['Q']:.6e}",
                            f"{c['D']:+.6e}", f"{c['z_cell']:.4f}"])

    with open(EXT_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["r", "r_2", "r_3", "is_diagonal", "P", "Q",
                    "P_minus_Q", "z_cell"])
        for er in ext_results:
            w.writerow([er["r"], er["r2"], er["r3"], er["is_diag"],
                        f"{er['P']:.6e}", f"{er['Q']:.6e}",
                        f"{er['D']:+.6e}", f"{er['z_cell']:.4f}"])

    # Markdown
    L = []
    L.append("# Bohr Probe — Asymptote Characterization (k = 40..70)")
    L.append("")
    L.append(f"N = {N:,}; (a, b) = ({A}, {B}); d.f. = {DF}; "
             f"starts in [{2*LOW+1:,}, {2*HIGH-1:,}]")
    L.append("")
    L.append(f"d.f. note: standard chi² of independence on {ROWS}×{COLS} is "
             f"(rows−1)(cols−1) = **{DF}**. User's `8·54−1 = 431` looks like a "
             "typo (a=4 vs a=5). Using 795 to stay comparable with prior runs.")
    L.append("")

    L.append("## Per-k metrics")
    L.append("")
    L.append("| k | N_alive | survival | chi²/df | z | diag in top 10 |")
    L.append("|---:|---:|---:|---:|---:|---:|")
    for r in rows_data:
        L.append(f"| {r['k']} | {r['n_alive']:,} | "
                 f"{100*r['n_alive']/N:.2f}% | {r['cpdf']:.3f} | "
                 f"{r['z']:+.2f} | {r['diag']}/6 |")
    L.append("")

    L.append("## Asymptote fit chi²/df(k) = A − B·exp(−C·k)")
    L.append("")
    L.append(f"Fit method: `{fit_method}`")
    L.append("")
    L.append(f"- **A (asymptote)** = **{A_fit:.2f}**")
    L.append(f"- B = {B_fit:.2f}")
    L.append(f"- C = {C_fit:.5f}")
    L.append("")
    L.append("| k | actual chi²/df | model chi²/df | residual |")
    L.append("|---:|---:|---:|---:|")
    for r in rows_data:
        pred = A_fit - B_fit * np.exp(-C_fit * r["k"])
        L.append(f"| {r['k']} | {r['cpdf']:.3f} | {pred:.3f} | "
                 f"{r['cpdf']-pred:+.3f} |")
    L.append("")

    L.append("## Top 5 cells per depth")
    for k in DEPTHS:
        L.append("")
        L.append(f"### k = {k}")
        L.append("")
        L.append("| rank | r_2 | r_3 | P | Q | P−Q | z_cell | diag? |")
        L.append("|---:|---:|---:|---:|---:|---:|---:|:---:|")
        for rk, c in enumerate(top_per_depth[k][:5], 1):
            mark = "✓" if (c["r2"] == c["r3"] and
                           c["r2"] in KNOWN_DIAGONALS) else ""
            L.append(f"| {rk} | {c['r2']} | {c['r3']} | {c['P']:.4e} | "
                     f"{c['Q']:.4e} | {c['D']:+.3e} | {c['z_cell']:+.2f} | "
                     f"{mark} |")
    L.append("")

    L.append("## EXT — Diagonal extension at k=50")
    L.append("")
    L.append("Candidates r ∈ {7, 23, 25, 31, 37}. Each is the cell "
             "(r mod 32, r mod 81); r=37 lands at (5, 37) — not on the "
             "r2=r3 diagonal.")
    L.append("")
    L.append("| r | (r2, r3) | on r2=r3 diag? | P | Q | P−Q | z_cell |")
    L.append("|---:|---|:---:|---:|---:|---:|---:|")
    for er in ext_results:
        di = "✓" if er["is_diag"] else ""
        L.append(f"| {er['r']} | ({er['r2']}, {er['r3']}) | {di} | "
                 f"{er['P']:.4e} | {er['Q']:.4e} | {er['D']:+.3e} | "
                 f"{er['z_cell']:+.2f} |")
    L.append("")
    over_5 = sum(1 for er in ext_results if abs(er["z_cell"]) > 5)
    over_3 = sum(1 for er in ext_results if abs(er["z_cell"]) > 3)
    L.append(f"**Cells with |z_cell| > 5:** {over_5} of {len(ext_results)}")
    L.append(f"**Cells with |z_cell| > 3:** {over_3} of {len(ext_results)}")
    L.append("")

    MD_PATH.write_text("\n".join(L), encoding="utf-8")
    print(f"\n[done] runtime = {time.time()-t0:.1f}s; outputs in {OUTDIR}",
          flush=True)


if __name__ == "__main__":
    main()
