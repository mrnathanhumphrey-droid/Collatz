"""
probe_B_composition.py
======================
Iterate the inter-level forcing operator R_k composition across multiple
levels. The single-level R_k has σ_1 ≈ 0.6706 (k-stable across k=5,6,7) which
doesn't match any known empirical rate (rho_slow ≈ 0.83, rate-1/2, etc.).

This probe tests whether the composition R_(m) = R_{k+m-1} ∘ ι_{k+m-2} ∘ ...
∘ ι_k ∘ R_k: V_k → W_{k+m} produces a different effective per-level rate as
m grows. The geometric mean σ_1(R_(m))^(1/m) is the candidate Lyapunov rate
of the composition.

Embedding ι_k: W_{k+1} ↪ V_{k+1} is the orthonormal embedding (P_{W_{k+1}}^T).

Choose k_start = 4 to allow up to m=4 compositions while keeping K_8 the
largest matrix (n=4374, fits in memory).

Levels in composition: k_start = 4, with R applied at k=4, 5, 6, 7. So:
- R_(1) = R_4: V_4 → W_5
- R_(2) = R_5 ∘ ι_4 ∘ R_4: V_4 → W_6
- R_(3) = R_6 ∘ ι_5 ∘ R_5 ∘ ι_4 ∘ R_4: V_4 → W_7
- R_(4) = R_7 ∘ ι_6 ∘ R_6 ∘ ι_5 ∘ R_5 ∘ ι_4 ∘ R_4: V_4 → W_8

Output:
- probe_R_operator/probe_B_composition.csv
- probe_R_operator/probe_B_findings.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_R_operator"
Q = 3
K_START = 4
M_MAX = 4  # composes R_4, R_5, R_6, R_7


def order_of_two(N):
    assert N % 2 == 1
    m = 1
    v = 2 % N
    while v != 1:
        v = (v * 2) % N
        m += 1
    return m


def build_K_dense(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    pi = inv2
    for v in range(M):
        powers_inv2[v] = pi
        pi = (pi * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    weights = np.zeros(M, dtype=np.float64)
    for vv in range(min(M, 1074)):
        weights[vv] = 2.0 ** -(vv + 1)
    Z = weights.sum()
    weights /= Z
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
    return K


def build_T_lift(q, k_lower):
    k = k_lower
    N_k = q ** k
    N_kp1 = q ** (k + 1)
    coprime_k = [r for r in range(N_k) if r % q != 0]
    coprime_kp1 = [r for r in range(N_kp1) if r % q != 0]
    n_k = len(coprime_k)
    n_kp1 = len(coprime_kp1)
    state_idx_k = {r: i for i, r in enumerate(coprime_k)}
    state_idx_kp1 = {r: i for i, r in enumerate(coprime_kp1)}
    T = np.zeros((n_k, n_kp1), dtype=np.float64)
    inv_q = 1.0 / q
    for r_p in coprime_kp1:
        i_p = state_idx_kp1[r_p]
        r = r_p % N_k
        i = state_idx_k[r]
        T[i, i_p] = inv_q
    return T


def build_W_basis(T):
    Q_, _ = np.linalg.qr(T.T, mode="complete")
    n_k = T.shape[0]
    P_W = Q_[:, n_k:].T
    return P_W


def build_R_block(K_kp1, T, P_W):
    """R_k as col-vec linear map V_k → W_{k+1}, shape (dim_W, n_k)."""
    return P_W @ K_kp1.T @ T.T


def main():
    print("=" * 78)
    print(f"Probe B: composition of R operators starting at k={K_START}, "
          f"up to m={M_MAX} levels")
    print("=" * 78)

    # Build R_block, P_W (for embedding) at each level k = K_START, ..., K_START + M_MAX - 1
    R_blocks = {}    # k → R_k (col-vec form, shape (2n_k, n_k))
    P_W_basis = {}   # k → P_W rows orth basis of W_{k+1} in V_{k+1} (shape (2n_k, n_{k+1}))
    dims = {}        # k → (n_k, dim_W_{k+1})
    timings = {}

    for k in range(K_START, K_START + M_MAX):
        t0 = time.time()
        K_kp1 = build_K_dense(Q, k + 1)
        T_k = build_T_lift(Q, k)
        P_W = build_W_basis(T_k)
        R_k = build_R_block(K_kp1, T_k, P_W)
        timings[k] = time.time() - t0
        R_blocks[k] = R_k
        P_W_basis[k] = P_W
        dims[k] = (T_k.shape[0], P_W.shape[0])
        print(f"  k={k}: V_{k} dim={dims[k][0]}, W_{k+1} dim={dims[k][1]}, "
              f"R_{k} shape={R_k.shape}, t={timings[k]:.1f}s")

    # Compose: R_(m) = R_{K_START+m-1} @ P_{W_{K_START+m-2}}^T @ ... @ R_{K_START}
    # In col-vec form, with v ∈ V_{K_START}:
    #   m=1: out = R_{K_START} v ∈ W_{K_START+1}, shape (2 n_{K_START}, n_{K_START})
    #   m=2: out = R_{K_START+1} (P_W_{K_START}^T (R_{K_START} v)) ∈ W_{K_START+2}
    #   ...
    # As matrix product:
    #   R_(m) = R_{K_START+m-1} @ P_{W_{K_START+m-2}}^T @ R_{K_START+m-2} @ P_{W_{K_START+m-3}}^T @ ... @ R_{K_START}
    # No P_W transition for the last R since the output stays in W (not embedded).
    # And no P_W transition before the first R since input is already in V.

    print(f"\nBuilding composed forcing operators R_(1..{M_MAX}):")
    R_comp = {}  # m → composed matrix
    composed_sigma1 = {}
    composed_top10 = {}
    R_comp_full = R_blocks[K_START].copy()  # R_(1)
    R_comp[1] = R_comp_full
    for m in range(2, M_MAX + 1):
        # Embed previous output (in W) back into V via P_W^T, then apply next R
        # R_(m) = R_{K_START+m-1} @ P_W_{K_START+m-2}^T @ R_(m-1)
        prev_k = K_START + m - 2
        embed = P_W_basis[prev_k].T  # shape (n_{prev_k+1}, dim_W_{prev_k+1}) = (n_{K_START+m-1}, 2 n_{K_START+m-2})
        R_next = R_blocks[K_START + m - 1]
        # R_comp_full has shape (dim_W_{K_START+m-2}, n_{K_START}) = (2 n_{K_START+m-2}, n_{K_START})
        # embed shape: (n_{K_START+m-1}, 2 n_{K_START+m-2})
        # R_next shape: (2 n_{K_START+m-1}, n_{K_START+m-1})
        # Result shape: (2 n_{K_START+m-1}, n_{K_START})
        t0 = time.time()
        R_comp_full = R_next @ embed @ R_comp_full
        t_compose = time.time() - t0
        R_comp[m] = R_comp_full
        print(f"  m={m}: composed shape {R_comp_full.shape}, t={t_compose:.1f}s")

    # SVD of each composition; report sigma_1 + top values
    print(f"\nSVD of composed operators:")
    print(f"{'m':>3} {'shape':>14} {'sigma_1':>14} {'sigma_1^(1/m)':>14} "
          f"{'sigma_2':>14} {'rank':>6}")
    rows = []
    for m in sorted(R_comp.keys()):
        R = R_comp[m]
        t0 = time.time()
        sigma = np.linalg.svd(R, compute_uv=False)
        t_svd = time.time() - t0
        sigma_1 = float(sigma[0])
        per_level = sigma_1 ** (1.0 / m)
        rank_e = int(np.sum(sigma > 1e-10))
        composed_sigma1[m] = sigma_1
        composed_top10[m] = sigma[:10].tolist()
        print(f"  {m:>3} {str(R.shape):>14} {sigma_1:>14.6f} "
              f"{per_level:>14.6f} "
              f"{sigma[1]:>14.6f} {rank_e:>6} "
              f"(SVD t={t_svd:.1f}s)")
        rows.append({
            "m": m, "shape_rows": R.shape[0], "shape_cols": R.shape[1],
            "sigma_1": sigma_1, "sigma_1_per_level": per_level,
            "sigma_2": float(sigma[1]) if len(sigma) > 1 else float("nan"),
            "sigma_3": float(sigma[2]) if len(sigma) > 2 else float("nan"),
            "rank": rank_e,
        })

    # Save CSV
    out_csv = os.path.join(OUT_DIR, "probe_B_composition.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "m", "shape_rows", "shape_cols",
            "sigma_1", "sigma_1_per_level", "sigma_2", "sigma_3", "rank"
        ])
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10e}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"\n[csv: {out_csv}]")

    # Reference rates for matching
    rho_slow = 0.826934
    references = {
        "rho_slow (order-3 real root)": rho_slow,
        "rate-1/2": 0.5,
        "0.6706 (single-level σ_1)": 0.6706,
        "order-3 complex magnitude": 0.192080,
    }
    def closest(x, candidates):
        return min(candidates.items(), key=lambda kv: abs(x - kv[1]))

    # Markdown
    md_path = os.path.join(OUT_DIR, "probe_B_findings.md")
    lines = []
    lines.append("# Probe B — composition of R_k forcing operators")
    lines.append("")
    lines.append(f"Composing R_(m) = R_{K_START+M_MAX-1} ∘ ι_{K_START+M_MAX-2} ∘ "
                 f"... ∘ ι_{K_START} ∘ R_{K_START}: V_{K_START} → W_{K_START+M_MAX}, "
                 f"for m = 1..{M_MAX}.")
    lines.append("")
    lines.append("ι_k: W_{k+1} ↪ V_{k+1} is the orthonormal embedding "
                 "(P_{W_{k+1}}^T). Each composition step adds one renormalization "
                 "level. The geometric mean σ_1(R_(m))^(1/m) is the candidate "
                 "Lyapunov rate per level.")
    lines.append("")
    lines.append("## Per-m results")
    lines.append("")
    lines.append("| m | output level | shape | σ_1 | σ_1^(1/m) | σ_2 | rank |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in rows:
        out_level = K_START + r["m"]
        lines.append(f"| {r['m']} | W_{out_level} | "
                     f"{r['shape_rows']}×{r['shape_cols']} | "
                     f"{r['sigma_1']:.6f} | {r['sigma_1_per_level']:.6f} | "
                     f"{r['sigma_2']:.6f} | {r['rank']} |")
    lines.append("")

    # Lyapunov rate analysis
    lines.append("## Lyapunov rate of the composition")
    lines.append("")
    lines.append("If the composition has a stable per-level multiplicative rate "
                 "λ, then σ_1(R_(m))^(1/m) → λ as m grows.")
    lines.append("")
    lines.append("| m | σ_1^(1/m) | closest reference | distance |")
    lines.append("|---|---|---|---|")
    for r in rows:
        per_level = r["sigma_1_per_level"]
        label, refval = closest(per_level, references)
        lines.append(f"| {r['m']} | {per_level:.6f} | {label} ({refval:.4f}) | "
                     f"{abs(per_level - refval):.4f} |")
    lines.append("")

    # Cumulative scaling
    lines.append("## Decay rate σ_1(R_(m+1)) / σ_1(R_(m))")
    lines.append("")
    lines.append("Per-step σ_1 ratio (multiplicative gain when adding one "
                 "level to the composition).")
    lines.append("")
    lines.append("| m → m+1 | σ_1(R_(m)) | σ_1(R_(m+1)) | ratio |")
    lines.append("|---|---|---|---|")
    for i in range(len(rows) - 1):
        a = rows[i]["sigma_1"]
        b = rows[i + 1]["sigma_1"]
        lines.append(f"| {rows[i]['m']} → {rows[i+1]['m']} | "
                     f"{a:.6f} | {b:.6f} | {b/a:.6f} |")
    lines.append("")

    # Verdict
    lines.append("## Verdict")
    lines.append("")
    if M_MAX >= 2:
        sigma1_per_level_seq = [r["sigma_1_per_level"] for r in rows]
        last = sigma1_per_level_seq[-1]
        d_to_rho = abs(last - rho_slow)
        d_to_065 = abs(last - 0.6706)
        d_to_half = abs(last - 0.5)
        lines.append(f"σ_1^(1/m) sequence: " +
                     ", ".join(f"{x:.4f}" for x in sigma1_per_level_seq))
        lines.append("")
        if d_to_rho < 0.05:
            lines.append(f"**σ_1^(1/m) approaches ρ_slow ≈ {rho_slow:.4f}** "
                         f"(distance at m={M_MAX}: {d_to_rho:.4f}). The "
                         "composition's per-level Lyapunov rate matches the "
                         "order-3 recurrence's slow root. Strong evidence "
                         "the multi-level R chain carries the rate-determining "
                         "structure for ε_k.")
        elif d_to_065 < 0.02:
            lines.append(f"**σ_1^(1/m) stays near 0.6706** "
                         f"(distance at m={M_MAX}: {d_to_065:.4f}). The "
                         "composition is essentially a power of the single-level "
                         "R, suggesting the R_k operators commute "
                         "approximately or share dominant directions across "
                         "levels. Composition does not reveal a new rate.")
        elif d_to_half < 0.02:
            lines.append(f"**σ_1^(1/m) drops toward 0.5** (distance at m={M_MAX}: "
                         f"{d_to_half:.4f}). Per-level rate matches the "
                         "legacy rate-1/2 envelope, suggesting the composition "
                         "decays faster than single-level σ_1.")
        else:
            lines.append(f"**No clean match.** σ_1^(1/m) at m={M_MAX} = {last:.4f}, "
                         f"closest reference {closest(last, references)[0]} at "
                         f"distance {abs(last - closest(last, references)[1]):.4f}. "
                         "Reporting and pausing for analysis as per the brief's "
                         "third walk-back gate.")
    lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append("- `probe_B_composition.csv` — per-m sigma values")
    lines.append("- `probe_B_findings.md` — this writeup")

    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"[md: {md_path}]")
    print("Done.")


if __name__ == "__main__":
    main()
