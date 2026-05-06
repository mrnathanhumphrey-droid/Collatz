"""
mode_amplitudes_probe.py
========================
Decompose pi_k onto top-20 K_k eigenvectors at k=5,6,7. Track mode amplitudes
across k via lift-bridge. K_k built per the preflight machinery (q=3, Syracuse
3x+1 chain on coprime mod 3^k).

For each k:
  - Build K_k (dense float64; n in {162, 486, 1458} all tractable)
  - Compute pi_k via power iteration if not cached on disk
  - Compute top-20 right eigenpairs via scipy.linalg.eig + take top by |lam|
    (dense full eig is faster than sparse Arnoldi at these sizes; brief
    permits dense fallback)
  - Project a_i = <v_i, pi_k> (Hermitian inner product, brief's prescription)
  - Report reconstruction error ||pi_k - sum_i a_i v_i||_inf

Mode tracking across k via lift L_k: v[r'] = v[r' mod 3^k] for r' coprime
in Z/3^(k+1). Mode i at level k matches mode j at level k+1 by argmax_j
|<L_k v_i^(k) / ||L_k v_i||, v_j^(k+1)>|.

Outputs in C:\\Collatz\\probe_mode_amplitudes\\:
  pi_k{k}.npy  per k (5,6,7)
  mode_amplitudes_k{k}.csv  per k
  mode_amplitudes_comparison.csv  (cross-k, lift-tracked)
  mode_amplitudes_findings.md
"""
from __future__ import annotations

import csv
import os
import sys
import time

import numpy as np
import scipy.linalg as la

sys.stdout.reconfigure(encoding="utf-8")

OUTDIR = r"C:\Collatz\probe_mode_amplitudes"
os.makedirs(OUTDIR, exist_ok=True)

# Empirical eps_k from prior probes
EPS = {
    1: 0.2,
    2: 1/105,
    3: -5191/1019445,
    4: -2.4522582483e-3,
    5: -1.1517469151e-3,
    6: -4.9790566522e-4,
    7: -1.1752368304e-3,
}


# -------- K_k build (from result_epsilon_7.py) --------

def order_of_two(N):
    m = 1; v = 2 % N
    while v != 1:
        v = (v * 2) % N; m += 1
    return m


def build_K_float(q, k):
    N = q ** k
    M = order_of_two(N)
    inv2 = pow(2, -1, N)
    powers_inv2 = np.empty(M, dtype=np.int64)
    p = inv2
    for v in range(M):
        powers_inv2[v] = p
        p = (p * inv2) % N
    coprime = np.array([r for r in range(N) if r % q != 0], dtype=np.int64)
    n = len(coprime)
    state_idx = -np.ones(N, dtype=np.int64)
    for i, r in enumerate(coprime):
        state_idx[r] = i
    K = np.zeros((n, n), dtype=np.float64)
    Z_v = 1.0 - 2.0 ** (-M)
    weights = np.array([(2.0 ** (-v)) / Z_v for v in range(1, M + 1)],
                       dtype=np.float64)
    for i_r, r in enumerate(coprime):
        base = (q * int(r) + 1) % N
        for v in range(1, M + 1):
            tgt = (base * int(powers_inv2[v - 1])) % N
            j = int(state_idx[tgt])
            K[i_r, j] += weights[v - 1]
    return K, coprime


def stationary_power_iter(K, max_iter=10000, tol=1e-12):
    n = K.shape[0]
    pi = np.full(n, 1.0 / n, dtype=np.float64)
    for it in range(max_iter):
        pi_new = pi @ K
        pi_new /= pi_new.sum()
        residual = float(np.max(np.abs(pi_new - pi)))
        pi = pi_new
        if residual < tol:
            return pi, it + 1, residual
    return pi, max_iter, residual


def top20_eigpairs(K):
    """Return top-20 eigenvalues and right eigenvectors of K, by |lambda| desc."""
    lams, vecs = la.eig(K)
    idx = np.argsort(-np.abs(lams))[:20]
    return lams[idx], vecs[:, idx]


def project_amplitudes(vecs, pi):
    """a_i = <v_i, pi>; return complex array of shape (20,)."""
    return vecs.conj().T @ pi.astype(np.complex128)


def lift_vector(v_lower, k):
    """L_k v at level k+1: v[r'] = v[r' mod 3^k] for r' coprime in Z/3^(k+1)."""
    Nk = 3 ** k
    Nk1 = 3 ** (k + 1)
    coprime_low = [r for r in range(Nk) if r % 3 != 0]
    coprime_up = [r for r in range(Nk1) if r % 3 != 0]
    idx_low = {r: i for i, r in enumerate(coprime_low)}
    out = np.zeros(len(coprime_up), dtype=v_lower.dtype)
    for i_up, r_up in enumerate(coprime_up):
        out[i_up] = v_lower[idx_low[r_up % Nk]]
    return out


def main():
    print("=" * 78)
    print("Mode amplitudes probe: pi_k decomposed onto top-20 K_k eigenvectors")
    print("=" * 78)

    data = {}
    for k in [5, 6, 7]:
        print(f"\n--- k = {k} ---")
        # Build K
        t0 = time.time()
        K, coprime = build_K_float(3, k)
        n = K.shape[0]
        print(f"  build K_{k} ({n}x{n}): {time.time()-t0:.2f}s")

        # Stationary
        pi_path = os.path.join(OUTDIR, f"pi_k{k}.npy")
        if os.path.exists(pi_path):
            pi = np.load(pi_path)
            print(f"  pi_{k}: loaded cached from {pi_path}")
            pi_iters = "cached"
            pi_residual = float("nan")
        else:
            t0 = time.time()
            pi, pi_iters, pi_residual = stationary_power_iter(K)
            np.save(pi_path, pi)
            print(f"  pi_{k}: power iter {pi_iters} steps, residual "
                  f"{pi_residual:.2e}, {time.time()-t0:.2f}s, saved {pi_path}")

        # Top-20 eigenpairs
        t0 = time.time()
        lams, vecs = top20_eigpairs(K)
        print(f"  top-20 eigpairs (dense full eig): {time.time()-t0:.2f}s")
        print(f"    lam_1 = {lams[0].real:+.6f} + {lams[0].imag:+.6f}i  "
              f"|.| = {abs(lams[0]):.10f}")
        print(f"    lam_2 = {lams[1].real:+.6e} + {lams[1].imag:+.6e}i  "
              f"|.| = {abs(lams[1]):.6e}")

        # Projection
        amps = project_amplitudes(vecs, pi)
        recon = (vecs @ amps).real
        recon_err = float(np.max(np.abs(pi - recon)))
        recon_rel = recon_err / float(np.max(np.abs(pi)))
        print(f"  projection a_i = <v_i, pi>:")
        print(f"    a_1 = {amps[0].real:+.6e} + {amps[0].imag:+.6e}i  "
              f"|.| = {abs(amps[0]):.6e}")
        print(f"    |a_2| = {abs(amps[1]):.6e}")
        print(f"    reconstruction err ||pi - sum a_i v_i||_inf = "
              f"{recon_err:.6e}  (rel to ||pi||_inf: {recon_rel:.4f})")

        data[k] = {
            "n": n, "K": K, "coprime": coprime,
            "pi": pi, "pi_iters": pi_iters, "pi_residual": pi_residual,
            "lams": lams, "vecs": vecs,
            "amps": amps, "recon_err": recon_err, "recon_rel": recon_rel,
        }

        # Per-k CSV
        var_total = float(np.sum(np.abs(amps) ** 2))
        out_csv = os.path.join(OUTDIR, f"mode_amplitudes_k{k}.csv")
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["rank", "lambda_real", "lambda_imag", "magnitude",
                        "argument_rad", "amplitude_real", "amplitude_imag",
                        "abs_amplitude", "abs_amplitude_sq",
                        "variance_fraction"])
            for i in range(20):
                lam = lams[i]; a = amps[i]
                w.writerow([i + 1, lam.real, lam.imag,
                            abs(lam), float(np.angle(lam)),
                            a.real, a.imag, abs(a), abs(a) ** 2,
                            (abs(a) ** 2) / var_total])
        print(f"    saved {out_csv}")

    # ---- Mode tracking across k via lift bridge ----
    print("\n--- Mode tracking via lift L_k ---")
    overlaps = {}
    for k_low in [5, 6]:
        k_up = k_low + 1
        n_low = data[k_low]["vecs"].shape[1]
        n_up = data[k_up]["vecs"].shape[1]
        M = np.zeros((n_low, n_up), dtype=np.complex128)
        for i in range(n_low):
            v_lift = lift_vector(data[k_low]["vecs"][:, i], k_low)
            v_lift /= np.linalg.norm(v_lift)
            for j in range(n_up):
                M[i, j] = v_lift.conj() @ data[k_up]["vecs"][:, j]
        overlaps[(k_low, k_up)] = M
        # Best match per mode at k_low
        print(f"  k={k_low} -> k={k_up} match table (best j per i, |overlap|):")
        print(f"    {'i':>3}  {'j_best':>6}  {'|overlap|':>10}  {'flag':>10}")
        for i in range(20):
            row = np.abs(M[i])
            j_best = int(np.argmax(row))
            ov = float(row[j_best])
            flag = "ok" if ov >= 0.5 else "AMBIGUOUS"
            print(f"    {i+1:>3}  {j_best+1:>6}  {ov:>10.4f}  {flag:>10}")

    # Comparison CSV
    out_cmp = os.path.join(OUTDIR, "mode_amplitudes_comparison.csv")
    with open(out_cmp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank_at_k5", "lambda_at_k5", "|a_5|^2",
                    "j_match_at_k6", "overlap_56", "lambda_at_k6", "|a_6|^2",
                    "j_match_at_k7", "overlap_67", "lambda_at_k7", "|a_7|^2",
                    "trend_5_to_7"])
        for i in range(20):
            a5_sq = abs(data[5]["amps"][i]) ** 2
            row56 = np.abs(overlaps[(5, 6)][i])
            j6 = int(np.argmax(row56))
            ov56 = float(row56[j6])
            lam6 = data[6]["lams"][j6]
            a6_sq = abs(data[6]["amps"][j6]) ** 2
            row67 = np.abs(overlaps[(6, 7)][j6])
            j7 = int(np.argmax(row67))
            ov67 = float(row67[j7])
            lam7 = data[7]["lams"][j7]
            a7_sq = abs(data[7]["amps"][j7]) ** 2
            # trend
            if a5_sq > 0 and a6_sq > 0 and a7_sq > 0:
                if a7_sq > a6_sq > a5_sq:
                    trend = "growing"
                elif a7_sq < a6_sq < a5_sq:
                    trend = "decaying"
                else:
                    trend = "non-monotone"
            else:
                trend = "zero"
            w.writerow([i + 1, abs(data[5]["lams"][i]), a5_sq,
                        j6 + 1, ov56, abs(lam6), a6_sq,
                        j7 + 1, ov67, abs(lam7), a7_sq,
                        trend])
    print(f"\n  saved {out_cmp}")

    # ---- Pre-registered questions ----
    print("\n--- Pre-registered questions ---")

    # Q1: Mode-crossing
    print("\nQ1: Mode-crossing of dominant non-trivial mode (rank 2/3 at k=5)")
    a5_nt = np.abs(data[5]["amps"])[1:]  # skip Perron
    dom_nt_5 = int(np.argmax(a5_nt))  # rank in non-trivial set
    print(f"  k=5 dominant non-trivial: rank {dom_nt_5 + 2} (|a|² = "
          f"{a5_nt[dom_nt_5]**2:.4e})")
    # Track this mode through lifts
    j6_dom = int(np.argmax(np.abs(overlaps[(5, 6)][dom_nt_5 + 1])))
    print(f"  -> at k=6 lifts to mode {j6_dom + 1} "
          f"(overlap {abs(overlaps[(5,6)][dom_nt_5+1, j6_dom]):.4f})")
    # What's the dominant non-trivial at k=6 via direct projection?
    a6_nt = np.abs(data[6]["amps"])[1:]
    dom_nt_6 = int(np.argmax(a6_nt))
    print(f"  k=6 dominant non-trivial (direct): rank {dom_nt_6 + 2} "
          f"(|a|² = {a6_nt[dom_nt_6]**2:.4e})")
    a7_nt = np.abs(data[7]["amps"])[1:]
    dom_nt_7 = int(np.argmax(a7_nt))
    print(f"  k=7 dominant non-trivial (direct): rank {dom_nt_7 + 2} "
          f"(|a|² = {a7_nt[dom_nt_7]**2:.4e})")

    crossing = "NO" if (dom_nt_5 == dom_nt_6 == dom_nt_7) else "YES"
    print(f"  -> mode-crossing across k=5,6,7: {crossing}")

    # Q2: amplitude growth/decay with k
    print("\nQ2: Non-trivial mode amplitude trends")
    # via the comparison CSV's tracked rows
    growing = 0; decaying = 0; non_mono = 0
    with open(out_cmp) as f:
        for row in csv.DictReader(f):
            if int(row["rank_at_k5"]) == 1:
                continue
            t = row["trend_5_to_7"]
            if t == "growing": growing += 1
            elif t == "decaying": decaying += 1
            else: non_mono += 1
    print(f"  of 19 non-trivial modes: growing={growing}, decaying={decaying}, "
          f"non-monotone={non_mono}")

    # Q3: connection to eps_k
    print("\nQ3: Does any single mode's |a_i|² track |eps_k|?")
    eps_vec = np.array([abs(EPS[5]), abs(EPS[6]), abs(EPS[7])])
    print(f"  reference: |eps_5|, |eps_6|, |eps_7| = "
          f"{eps_vec[0]:.4e}, {eps_vec[1]:.4e}, {eps_vec[2]:.4e}")
    print(f"  ratio eps_6/eps_5 = {eps_vec[1]/eps_vec[0]:.4f}, "
          f"eps_7/eps_6 = {eps_vec[2]/eps_vec[1]:.4f}")
    print()
    # For each lift-tracked mode, compute Pearson correlation of (|a_5|^2,
    # |a_6|^2, |a_7|^2) with (|eps_5|^2, |eps_6|^2, |eps_7|^2)? With only
    # 3 points correlation is weak; instead just print which modes have a
    # qualitatively-similar non-monotone pattern (smaller at k=6 than k=5
    # or k=7).
    matches = []
    with open(out_cmp) as f:
        for row in csv.DictReader(f):
            i = int(row["rank_at_k5"])
            if i == 1:
                continue
            a5_sq = float(row["|a_5|^2"])
            a6_sq = float(row["|a_6|^2"])
            a7_sq = float(row["|a_7|^2"])
            # eps pattern: |eps_6| < |eps_5|, |eps_7| > |eps_6|, |eps_7| ≈ |eps_5|
            if a5_sq > 0 and a6_sq > 0 and a7_sq > 0:
                if a6_sq < a5_sq and a7_sq > a6_sq:
                    rel_eps5 = abs(np.log(a5_sq + 1e-30) - np.log(eps_vec[0]**2))
                    matches.append((i, a5_sq, a6_sq, a7_sq, rel_eps5))
    if matches:
        print(f"  Modes with qualitatively similar non-monotone shape (a_6 < a_5, a_7 > a_6):")
        for i, a5, a6, a7, _ in matches[:10]:
            print(f"    rank {i}: |a|² = {a5:.4e}, {a6:.4e}, {a7:.4e}")
    else:
        print("  no modes found with non-monotone shape matching eps_k")
    print()

    # Findings markdown
    md = []
    md.append("# Result: mode amplitude decomposition of pi_k onto K_k eigenvectors")
    md.append("")
    md.append("**Date:** 2026-05-05.")
    md.append("")
    md.append("## Setup")
    md.append("")
    md.append("| k | n_states | pi power-iter steps | top-20 eigvalue range |")
    md.append("|---|---|---|---|")
    for k in [5, 6, 7]:
        d = data[k]
        lam_max = abs(d["lams"][0]); lam_min = abs(d["lams"][19])
        md.append(f"| {k} | {d['n']} | {d['pi_iters']} | "
                  f"{lam_max:.4e} ↓ {lam_min:.4e} |")
    md.append("")
    md.append("## Reconstruction quality (top-20 spans pi_k how well)")
    md.append("")
    md.append("| k | ||pi - Σ a_i v_i||_inf | rel to ||pi||_inf |")
    md.append("|---|---|---|")
    for k in [5, 6, 7]:
        md.append(f"| {k} | {data[k]['recon_err']:.4e} | "
                  f"{data[k]['recon_rel']:.4f} |")
    md.append("")
    md.append("## Per-k top-5 amplitudes |a_i|² and variance fractions")
    md.append("")
    for k in [5, 6, 7]:
        d = data[k]
        var_total = float(np.sum(np.abs(d["amps"]) ** 2))
        md.append(f"### k = {k}")
        md.append("")
        md.append("| rank | |λ| | arg(λ) | |a_i| | |a_i|² | variance frac |")
        md.append("|---|---|---|---|---|---|")
        for i in range(5):
            lam = d["lams"][i]; a = d["amps"][i]
            md.append(f"| {i+1} | {abs(lam):.4e} | {float(np.angle(lam)):+.4f} | "
                      f"{abs(a):.4e} | {abs(a)**2:.4e} | "
                      f"{(abs(a)**2)/var_total:.4f} |")
        md.append("")

    md.append("## Pre-registered questions")
    md.append("")
    md.append(f"### Q1: Mode-crossing")
    md.append("")
    md.append(f"- k=5 dominant non-trivial: rank {dom_nt_5 + 2}, "
              f"|a|² = {a5_nt[dom_nt_5]**2:.4e}")
    md.append(f"- k=6 dominant non-trivial: rank {dom_nt_6 + 2}, "
              f"|a|² = {a6_nt[dom_nt_6]**2:.4e}")
    md.append(f"- k=7 dominant non-trivial: rank {dom_nt_7 + 2}, "
              f"|a|² = {a7_nt[dom_nt_7]**2:.4e}")
    md.append(f"- Mode-crossing observed: **{crossing}**")
    md.append("")
    md.append(f"### Q2: Amplitude growth/decay with k")
    md.append("")
    md.append(f"Of 19 non-trivial lift-tracked modes:")
    md.append(f"- growing (a_5 < a_6 < a_7): **{growing}**")
    md.append(f"- decaying (a_5 > a_6 > a_7): **{decaying}**")
    md.append(f"- non-monotone: **{non_mono}**")
    md.append("")
    md.append(f"### Q3: Connection to ε_k")
    md.append("")
    md.append(f"Reference: |ε_5|, |ε_6|, |ε_7| = "
              f"{eps_vec[0]:.4e}, {eps_vec[1]:.4e}, {eps_vec[2]:.4e}.")
    md.append(f"|ε_k| has a local minimum at k=6 and bounces back at k=7 "
              f"(ratio 0.43, 2.36).")
    md.append("")
    if matches:
        md.append(f"Modes with same qualitative non-monotone shape (a_6 < a_5, "
                  f"a_7 > a_6): {len(matches)} found:")
        for i, a5, a6, a7, _ in matches[:10]:
            md.append(f"- rank {i} at k=5: |a|² = {a5:.4e}, {a6:.4e}, "
                      f"{a7:.4e}")
    else:
        md.append("No modes found with non-monotone shape matching ε_k. The "
                  "eps_k non-monotone pattern is NOT carried by a single mode "
                  "in the top-20.")
    md.append("")
    md.append("## Files")
    md.append("")
    md.append("- `mode_amplitudes_probe.py` — script")
    md.append("- `pi_k{5,6,7}.npy` — cached stationary distributions")
    md.append("- `mode_amplitudes_k{5,6,7}.csv` — per-k amplitude tables")
    md.append("- `mode_amplitudes_comparison.csv` — cross-k tracked via lift")
    md.append("- `mode_amplitudes_findings.md` — this writeup")

    out_md = os.path.join(OUTDIR, "mode_amplitudes_findings.md")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    print(f"saved {out_md}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
