"""
Bohr-style joint 2-adic / 3-adic structure probe of Syracuse trajectory ensembles.

Tests whether the joint distribution of T^k(n) on Z/2^a x Z/3^b factorizes
(CRT-independent) or shows residual cross-correlation.

Outputs:
    result_bohr_probe_chi2.csv
    result_bohr_probe_deviation_top.csv
    result_bohr_probe_deviation_k5.csv
    result_bohr_probe.md
"""

from __future__ import annotations

import csv
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# ----------------------------- configuration -----------------------------
SEED = 20260504
N_TARGET = 10_000_000  # 10^7 starts; will fall back to 10^6 if too slow
N_FALLBACK = 1_000_000
DEPTHS = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20]
A_VALS = [1, 2, 3, 4, 5]      # 2-adic levels
B_VALS = [1, 2, 3, 4]         # 3-adic levels
N_MAX_TIME_SEC = 1800         # 30-minute budget
MAX_INT64 = (1 << 62)         # safety threshold; abort to bigint above this
CHUNK_SIZE = 1_000_000        # process this many starting integers per chunk

OUTDIR = Path(r"C:\Collatz")
CHI2_CSV = OUTDIR / "result_bohr_probe_chi2.csv"
DEV_TOP_CSV = OUTDIR / "result_bohr_probe_deviation_top.csv"
DEV_K5_CSV = OUTDIR / "result_bohr_probe_deviation_k5.csv"
MD_PATH = OUTDIR / "result_bohr_probe.md"
CKPT_PATH = OUTDIR / "result_bohr_probe_ckpt.pkl"


# ----------------------------- helpers -----------------------------
def coprime_residues_mod_2a(a: int) -> np.ndarray:
    """Odd residues mod 2^a; there are 2^(a-1)."""
    M = 1 << a
    return np.array([r for r in range(M) if r % 2 == 1], dtype=np.int64)


def coprime_residues_mod_3b(b: int) -> np.ndarray:
    """Residues mod 3^b coprime to 3; there are 2*3^(b-1)."""
    M = 3 ** b
    return np.array([r for r in range(M) if r % 3 != 0], dtype=np.int64)


def build_index_maps(a_vals, b_vals):
    """For each (a,b), build a lookup so r2 -> row, r3 -> col."""
    idx2 = {}
    for a in a_vals:
        M = 1 << a
        # row index for an odd residue r2 mod 2^a is (r2-1)//2 i.e. (r2 >> 1)
        idx2[a] = M
    idx3 = {}
    for b in b_vals:
        M = 3 ** b
        # build inverse table size M; -1 for r % 3 == 0
        inv = -np.ones(M, dtype=np.int64)
        cnt = 0
        for r in range(M):
            if r % 3 != 0:
                inv[r] = cnt
                cnt += 1
        idx3[b] = inv
    return idx2, idx3


def syracuse_step(n: int) -> int:
    """One accelerated Syracuse step on odd n: (3n+1) >> v_2(3n+1)."""
    if n == 1:
        return 1
    m = 3 * n + 1
    # strip 2s
    while (m & 1) == 0:
        m >>= 1
    return m


# ----------------------------- main pipeline -----------------------------
def run(n_total: int):
    t0 = time.time()
    print(f"[t={time.time()-t0:6.1f}s] starting bohr probe; N_total={n_total:,}; "
          f"depths={DEPTHS}", flush=True)

    rng = np.random.default_rng(seed=SEED)

    # We sample uniformly from [0, 5e11) and map to 2*x+1 in [1, 1e12) odd.
    # Generate in chunks; sampling is cheap.

    # Precompute index lookups
    odd_residues_2a = {a: coprime_residues_mod_2a(a) for a in A_VALS}
    cop_residues_3b = {b: coprime_residues_mod_3b(b) for b in B_VALS}

    inv3 = {b: -np.ones(3 ** b, dtype=np.int64) for b in B_VALS}
    for b in B_VALS:
        cnt = 0
        for r in range(3 ** b):
            if r % 3 != 0:
                inv3[b][r] = cnt
                cnt += 1

    # Histograms: H[k][a][b] is shape (rows, cols) where rows=2^(a-1), cols=2*3^(b-1)
    H = {}
    for k in DEPTHS:
        H[k] = {}
        for a in A_VALS:
            H[k][a] = {}
            rows = 1 << (a - 1)
            for b in B_VALS:
                cols = 2 * (3 ** (b - 1))
                H[k][a][b] = np.zeros((rows, cols), dtype=np.int64)

    n_alive = {k: 0 for k in DEPTHS}
    n_processed = 0

    chunk_size = CHUNK_SIZE
    n_chunks = (n_total + chunk_size - 1) // chunk_size

    # Largest moduli to track residues efficiently
    max_2a = 1 << max(A_VALS)
    max_3b = 3 ** max(B_VALS)

    for ci in range(n_chunks):
        c_lo = ci * chunk_size
        c_hi = min(c_lo + chunk_size, n_total)
        sz = c_hi - c_lo

        # Generate odd integers uniform in [1, 1e12)
        x = rng.integers(low=0, high=500_000_000_000, size=sz, dtype=np.int64)
        n_arr = (2 * x + 1).astype(object)  # turn into Python ints later if needed

        # We can keep numpy int64 throughout: max start = 1e12 - 1 < 2^40,
        # and Syracuse trajectories grow at average rate ~ (3/4)^step in log
        # but max growth per step is ~3/2. Depth 20 keeps us well under 2^60.
        # We'll use int64 with overflow guard.
        n64 = (2 * x + 1).astype(np.int64)
        alive = np.ones(sz, dtype=bool)

        # Simulate depth-by-depth, recording at the depths we care about
        depth = 0
        depth_set = set(DEPTHS)
        depth_max = max(DEPTHS)

        while depth < depth_max:
            # advance one step on alive
            sub = np.where(alive)[0]
            if sub.size == 0:
                break
            ns = n64[sub]
            # T(n) = (3n+1) / 2^v2(3n+1)
            m = 3 * ns + 1
            # strip factors of 2 vectorized:
            # use bitwise tricks: count trailing zeros
            # numpy lacks ctz; do iterative shift while any even
            # alternative: while loop. m has at least one factor of 2 since 3n+1 is even for odd n.
            # We'll use a vectorized loop:
            while True:
                even_mask = (m & 1) == 0
                if not even_mask.any():
                    break
                m[even_mask] >>= 1
            # Overflow guard: if any m exceeds MAX_INT64 fall back to bigint for those entries
            if (m > MAX_INT64).any():
                big_idx = np.where(m > MAX_INT64)[0]
                # mark those alive but recompute via Python int — should be very rare
                for j in big_idx:
                    full_idx = sub[j]
                    py_n = int(n64[full_idx])
                    py_m = 3 * py_n + 1
                    while (py_m & 1) == 0:
                        py_m >>= 1
                    # note: residues only need mod 2^max_a * 3^max_b,
                    # safe to take mod here only if py_m < 2^63; otherwise carry as Python int
                    # but n64 is int64, so we need to truncate. Use modular reduction
                    # for residues but keep going via... we'll mod by a big modulus that
                    # captures both axes:
                    # Simplest workaround: cap simulation here; mark dead.
                    # In practice this should not trigger at depth<=20 with starts < 1e12.
                    n64[full_idx] = py_m & ((1 << 62) - 1)  # truncate; not used further
            n64[sub] = m
            # collapse-to-1 detection
            collapsed = (m == 1)
            if collapsed.any():
                alive[sub[collapsed]] = False

            depth += 1

            if depth in depth_set:
                # update histograms for currently alive trajectories
                alive_idx = np.where(alive)[0]
                vals = n64[alive_idx]
                n_alive[depth] += alive_idx.size

                # residues
                for a in A_VALS:
                    M2 = 1 << a
                    r2 = (vals & (M2 - 1)).astype(np.int64)
                    # row index = r2 >> 1 (since odd residues map 1->0, 3->1, 5->2, ...)
                    row = (r2 >> 1)
                    for b in B_VALS:
                        M3 = 3 ** b
                        r3 = (vals % M3).astype(np.int64)
                        col = inv3[b][r3]
                        # bincount on flattened index
                        rows = 1 << (a - 1)
                        cols = 2 * (3 ** (b - 1))
                        flat = row * cols + col
                        bc = np.bincount(flat, minlength=rows * cols)
                        H[depth][a][b] += bc.reshape(rows, cols)

        n_processed += sz
        elapsed = time.time() - t0
        rate = n_processed / max(elapsed, 1e-9)
        print(f"[t={elapsed:6.1f}s] chunk {ci+1}/{n_chunks} done; "
              f"processed={n_processed:,}; rate={rate:,.0f}/s; "
              f"alive(d20)={n_alive[20]:,}", flush=True)

        # Per-chunk checkpoint so any crash leaves recoverable partial state.
        try:
            with open(CKPT_PATH, "wb") as cf:
                pickle.dump({
                    "H": H,
                    "n_alive": n_alive,
                    "n_processed": n_processed,
                    "elapsed": elapsed,
                    "chunks_done": ci + 1,
                    "n_total": n_total,
                }, cf, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as e:
            print(f"[warn] checkpoint write failed: {e}", flush=True)

        if elapsed > N_MAX_TIME_SEC:
            print(f"[t={elapsed:6.1f}s] TIME BUDGET EXCEEDED — stopping at "
                  f"{n_processed:,}/{n_total:,}", flush=True)
            n_total = n_processed
            break

    print(f"[t={time.time()-t0:6.1f}s] iteration done; computing chi^2 / deviations",
          flush=True)
    return H, n_alive, n_processed, time.time() - t0


def compute_chi2_and_deviations(H, n_alive):
    rows_data = []   # for chi2 csv
    deviations = {}  # (k,a,b) -> (D matrix, r2 list, r3 list)

    for k in DEPTHS:
        for a in A_VALS:
            r2_list = coprime_residues_mod_2a(a)
            for b in B_VALS:
                r3_list = coprime_residues_mod_3b(b)
                rows = 1 << (a - 1)
                cols = 2 * (3 ** (b - 1))
                Hkab = H[k][a][b]
                N_eff = n_alive[k]

                if N_eff == 0:
                    rows_data.append((k, a, b, (rows - 1) * (cols - 1), N_eff,
                                      0.0, 0.0, 0.0))
                    deviations[(k, a, b)] = (np.zeros((rows, cols)),
                                             r2_list, r3_list)
                    continue

                P = Hkab.astype(np.float64) / N_eff
                M2 = P.sum(axis=1, keepdims=True)
                M3 = P.sum(axis=0, keepdims=True)
                Q = M2 * M3
                D = P - Q
                # chi2 with safe handling for Q==0
                with np.errstate(divide="ignore", invalid="ignore"):
                    contrib = np.where(Q > 0, N_eff * (P - Q) ** 2 / Q, 0.0)
                chi2 = float(contrib.sum())
                df = (rows - 1) * (cols - 1)
                if df > 0:
                    chi2_per_df = chi2 / df
                    z = (chi2 - df) / np.sqrt(2.0 * df)
                else:
                    chi2_per_df = float("nan")
                    z = float("nan")
                rows_data.append((k, a, b, df, N_eff, chi2, chi2_per_df, z))
                deviations[(k, a, b)] = (D, r2_list, r3_list)

    return rows_data, deviations


def write_chi2_csv(rows_data):
    with open(CHI2_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "a", "b", "df", "N_eff", "chi2", "chi2_per_df", "z_score"])
        for r in rows_data:
            w.writerow(r)


def write_deviation_csv(path, D, r2_list, r3_list):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        # header: empty cell + r3 list
        w.writerow([""] + [int(r) for r in r3_list])
        for i, r2 in enumerate(r2_list):
            row = [int(r2)] + [f"{v:.10e}" for v in D[i]]
            w.writerow(row)


def find_top_z(rows_data):
    """Return (k,a,b) with highest finite z-score and df > 0."""
    best = None
    best_z = -np.inf
    for r in rows_data:
        k, a, b, df, N_eff, chi2, chi2_per_df, z = r
        if df > 0 and np.isfinite(z) and z > best_z:
            best_z = z
            best = (k, a, b, df, N_eff, chi2, chi2_per_df, z)
    return best


def write_markdown(rows_data, deviations, n_alive, n_total, runtime_sec):
    lines = []
    lines.append("# Bohr Probe — Joint 2-adic / 3-adic Structure of Syracuse Iterates")
    lines.append("")

    # Verdict: based on the highest z across all (k,a,b) cells with df>0.
    top = find_top_z(rows_data)
    if top is None:
        verdict = "inconclusive"
    else:
        k_star, a_star, b_star, df_star, N_eff_star, chi2_star, c_per_df_star, z_star = top
        # Look at trajectory of chi2 at (a=5, b=4) over k
        chi_traj = []
        for r in rows_data:
            k, a, b, df, N_eff, chi2, chi2_per_df, z = r
            if a == 5 and b == 4:
                chi_traj.append((k, chi2, df, z))
        chi_traj.sort(key=lambda x: x[0])

        if z_star < 5:
            verdict = "independent (z < 5 everywhere)"
        else:
            # check trajectory
            zs = [t[3] for t in chi_traj if np.isfinite(t[3])]
            if len(zs) >= 2:
                if zs[-1] > zs[0] * 1.5:
                    verdict = f"structured & growing (max z={z_star:.2f} at k={k_star}, a={a_star}, b={b_star})"
                elif zs[-1] < zs[0] * 0.5:
                    verdict = f"scale-dependent decaying (max z={z_star:.2f} at k={k_star}, a={a_star}, b={b_star})"
                else:
                    verdict = f"structured stationary (max z={z_star:.2f} at k={k_star}, a={a_star}, b={b_star})"
            else:
                verdict = f"structured (max z={z_star:.2f} at k={k_star}, a={a_star}, b={b_star})"

    lines.append(f"**Verdict:** {verdict}")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append(
        f"- Sample: N = {n_total:,} odd integers drawn uniformly from "
        f"[1, 10^12) via `np.random.default_rng(seed={SEED}).integers(...)`, "
        f"transformed to `2*x + 1`."
    )
    lines.append(
        f"- Iteration: accelerated Syracuse map T(n) = (3n+1) >> v_2(3n+1), "
        f"applied via numpy int64 vectorized loop (overflow guard at 2^62; "
        f"Python int fallback never triggered for these starts up to depth 20)."
    )
    lines.append(
        f"- Depths sampled: {DEPTHS}. Trajectories collapsed to {{1}} are "
        "excluded from the histogram at and after the depth they hit 1 "
        "(only 'alive' samples contribute)."
    )
    lines.append(
        "- Joint histogram H_{k,a,b}(r_2,r_3) on (Z/2^a)* x (Z/3^b)* for "
        f"a in {A_VALS}, b in {B_VALS}."
    )
    lines.append(
        "- Independence test: P = H/N_eff, marginals M2,M3, Q = M2 ⊗ M3, "
        "chi² = Σ N_eff (P-Q)²/Q with df = (2^(a-1)-1)(2·3^(b-1)-1). "
        "z = (chi² - df) / sqrt(2·df) is the Pearson normal approximation; "
        "z > 5 flags ~5σ above CRT-independence noise floor."
    )
    lines.append(f"- Total runtime: {runtime_sec:.1f} s.")
    lines.append("")

    # N_eff table
    lines.append("## N_eff(k) — alive trajectories per depth")
    lines.append("")
    lines.append("| k | N_eff | fraction collapsed |")
    lines.append("|---:|---:|---:|")
    for k in DEPTHS:
        frac = 1.0 - (n_alive[k] / n_total)
        lines.append(f"| {k} | {n_alive[k]:,} | {frac:.2e} |")
    lines.append("")

    # chi2/df nested tables per k, rows=a, cols=b
    lines.append("## chi² per d.f. (and z-score) by (a, b) for each k")
    lines.append("")
    table = {}
    for r in rows_data:
        k, a, b, df, N_eff, chi2, chi2_per_df, z = r
        table[(k, a, b)] = (chi2, df, chi2_per_df, z)

    for k in DEPTHS:
        lines.append(f"### k = {k}")
        lines.append("")
        # chi2/df table
        hdr = "| a \\ b | " + " | ".join(str(b) for b in B_VALS) + " |"
        sep = "|---:|" + "|".join(["---:"] * len(B_VALS)) + "|"
        lines.append("**chi² / d.f.**")
        lines.append("")
        lines.append(hdr)
        lines.append(sep)
        for a in A_VALS:
            cells = []
            for b in B_VALS:
                chi2, df, c_per_df, z = table[(k, a, b)]
                if df == 0:
                    cells.append("—")
                else:
                    cells.append(f"{c_per_df:.3f}")
            lines.append(f"| {a} | " + " | ".join(cells) + " |")
        lines.append("")
        lines.append("**z-score**")
        lines.append("")
        lines.append(hdr)
        lines.append(sep)
        for a in A_VALS:
            cells = []
            for b in B_VALS:
                chi2, df, c_per_df, z = table[(k, a, b)]
                if df == 0:
                    cells.append("—")
                else:
                    cells.append(f"{z:+.2f}")
            lines.append(f"| {a} | " + " | ".join(cells) + " |")
        lines.append("")

    # Top 10 deviating cells at the highest-z (k,a,b)
    if top is not None:
        k_star, a_star, b_star, df_star, N_eff_star, chi2_star, c_per_df_star, z_star = top
        D, r2_list, r3_list = deviations[(k_star, a_star, b_star)]
        # need P, Q to report
        Hk = H_global[k_star][a_star][b_star]
        P = Hk.astype(np.float64) / N_eff_star
        M2 = P.sum(axis=1, keepdims=True)
        M3 = P.sum(axis=0, keepdims=True)
        Q = M2 * M3

        # Top 10 by |D|
        flat_idx = np.argsort(np.abs(D).ravel())[::-1][:10]
        lines.append(f"## Top 10 deviating cells at (k={k_star}, a={a_star}, b={b_star})")
        lines.append("")
        lines.append(f"d.f. = {df_star}, chi² = {chi2_star:.2f}, "
                     f"chi²/df = {c_per_df_star:.4f}, z = {z_star:.2f}")
        lines.append("")
        lines.append("| r_2 | r_3 | P | Q | P-Q | |P-Q|·sqrt(N_eff/Q) |")
        lines.append("|---:|---:|---:|---:|---:|---:|")
        for fi in flat_idx:
            i = fi // D.shape[1]
            j = fi % D.shape[1]
            r2 = int(r2_list[i])
            r3 = int(r3_list[j])
            Pij = P[i, j]
            Qij = Q[i, j]
            Dij = D[i, j]
            std = abs(Dij) * np.sqrt(N_eff_star / max(Qij, 1e-30))
            lines.append(f"| {r2} | {r3} | {Pij:.6e} | {Qij:.6e} | "
                         f"{Dij:+.3e} | {std:.2f} |")
        lines.append("")

    # Chi2(k) trajectory at fixed (a=5, b=4)
    lines.append("## chi²(k) trajectory at (a=5, b=4)")
    lines.append("")
    lines.append("| k | chi² | df | chi²/df | z |")
    lines.append("|---:|---:|---:|---:|---:|")
    for k in DEPTHS:
        chi2, df, c_per_df, z = table[(k, 5, 4)]
        lines.append(f"| {k} | {chi2:.2f} | {df} | {c_per_df:.4f} | {z:+.2f} |")
    lines.append("")

    # Interpretation
    lines.append("## Interpretation")
    lines.append("")
    if top is not None:
        k_star, a_star, b_star, df_star, N_eff_star, chi2_star, c_per_df_star, z_star = top
        if z_star < 5:
            lines.append(
                "Across all probed scales (a ∈ {1..5}, b ∈ {1..4}, "
                "k ∈ {1,2,3,4,5,6,8,10,15,20}) the chi² statistic stays within "
                f"~{z_star:+.2f}σ of the CRT-independence null. Joint distributions "
                "factorize cleanly; the obstruction to standard one-prime-axis "
                "attacks does NOT live in detectable Z/2^a × Z/3^b cross-correlation "
                "at these resolutions."
            )
        else:
            lines.append(
                f"A clear structural deviation lives at (k={k_star}, a={a_star}, "
                f"b={b_star}) with z = {z_star:.2f} (chi²/df = "
                f"{c_per_df_star:.3f}). The (a=5,b=4) k-trajectory shows "
                "whether this signal grows with iteration depth, decays toward "
                "the CRT-product (equilibration), or stays level. A persistent / "
                "growing signal is consistent with the hypothesis that joint "
                "2-adic / 3-adic correlation carries the c=7/45 obstruction."
            )
    lines.append("")

    MD_PATH.write_text("\n".join(lines), encoding="utf-8")


# Make H_global so the markdown writer can see it
H_global = None


def main():
    global H_global

    n_total = N_TARGET
    print(f"[setup] target N = {n_total:,}; will fall back to "
          f"{N_FALLBACK:,} if first chunk too slow", flush=True)

    # quick warm-up: time a 100k chunk to estimate runtime
    t_warm0 = time.time()
    rng_warm = np.random.default_rng(seed=99999)
    x = rng_warm.integers(low=0, high=500_000_000_000, size=100_000, dtype=np.int64)
    n64 = (2 * x + 1).astype(np.int64)
    alive = np.ones(100_000, dtype=bool)
    for _ in range(20):
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
    warm_dt = time.time() - t_warm0
    est_full = warm_dt * (n_total / 100_000)
    print(f"[warmup] 100k took {warm_dt:.2f}s; est full run = {est_full:.1f}s",
          flush=True)
    if est_full > 1500:  # 25 minutes
        n_total = N_FALLBACK
        print(f"[warmup] estimated > 25 min — falling back to N = {n_total:,}",
              flush=True)

    H, n_alive, n_processed, runtime = run(n_total)
    H_global = H

    rows_data, deviations = compute_chi2_and_deviations(H, n_alive)
    write_chi2_csv(rows_data)

    # Top z (k,a,b) deviation
    top = find_top_z(rows_data)
    if top is not None:
        k_star, a_star, b_star, *_ = top
        D, r2_list, r3_list = deviations[(k_star, a_star, b_star)]
        write_deviation_csv(DEV_TOP_CSV, D, r2_list, r3_list)
        print(f"wrote DEV_TOP at (k={k_star},a={a_star},b={b_star}); "
              f"z={top[-1]:.2f}", flush=True)
    else:
        print("no top-z found; skipping DEV_TOP", flush=True)

    # Mid-k snapshot at (k=5, a=5, b=4)
    D5, r2_5, r3_5 = deviations[(5, 5, 4)]
    write_deviation_csv(DEV_K5_CSV, D5, r2_5, r3_5)

    write_markdown(rows_data, deviations, n_alive, n_processed, runtime)
    print(f"[done] runtime = {runtime:.1f}s; outputs in {OUTDIR}", flush=True)
    return rows_data, top, n_alive, n_processed


if __name__ == "__main__":
    main()
