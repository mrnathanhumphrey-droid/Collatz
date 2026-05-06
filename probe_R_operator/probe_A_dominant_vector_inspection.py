"""
probe_A_dominant_vector_inspection.py
======================================
Inspect the dominant left singular vector u_1 in V_{k+1} natural basis
and the dominant right singular vector v_1 in V_k natural basis, for
the R_k forcing operator at k = 5, 6, 7.

For each k, compute:
- u_1 in W_{k+1} coordinates (from npz)
- u_1 in V_{k+1} natural basis: u_V = P_W^T @ u_W
- v_1 in V_k natural basis (already in that basis from npz)
- Per-residue concentration: |u_V|^2 grouped by r mod 3, mod 9, mod 27...
- Sign pattern + sparsity
- Comparison with pi_{k+1} (the stationary)

Outputs:
- probe_R_operator/probe_A_u1_support_k{k}.csv  (per-residue energy)
- probe_R_operator/probe_A_v1_support_k{k}.csv
- probe_R_operator/probe_A_findings.md
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
LEVELS = [5, 6, 7]
Q = 3
TOP_DIRECTIONS = 1  # u_1 and v_1 only for now (top singular vector pair)


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
    Z_v = 1.0 - 2.0 ** (-min(M, 1074))
    weights = np.zeros(M, dtype=np.float64)
    for vv in range(min(M, 1074)):
        weights[vv] = (2.0 ** -(vv + 1)) / Z_v
    K = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (q * r + 1) % N
        targets = (base * powers_inv2) % N
        js = state_idx[targets]
        K[i_r] = np.bincount(js, weights=weights, minlength=n)
    return K, coprime, state_idx


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
    return T, coprime_k, coprime_kp1


def build_W_basis(T):
    Q_, _ = np.linalg.qr(T.T, mode="complete")
    n_k = T.shape[0]
    P_W = Q_[:, n_k:].T
    return P_W


def stationary_vec(K, tol=1e-15, max_iter=20000):
    """Stationary measure: pi @ K = pi."""
    n = K.shape[0]
    pi = np.full(n, 1.0 / n)
    for it in range(max_iter):
        pi_new = pi @ K
        s = pi_new.sum()
        if s != 0:
            pi_new /= s
        r = float(np.linalg.norm(pi_new - pi, ord=1))
        pi = pi_new
        if r < tol:
            return pi
    return pi


def per_residue_energy(vec, coprime, modulus):
    """Sum of vec[i]^2 grouped by coprime[i] mod modulus.
    Returns dict residue → energy fraction."""
    energy = float(np.sum(vec ** 2))
    if energy == 0:
        return {}
    out = {}
    for r in range(modulus):
        if r % Q == 0:
            continue
        idx = np.where(coprime % modulus == r)[0]
        if len(idx) == 0:
            continue
        out[r] = float(np.sum(vec[idx] ** 2)) / energy
    return out


def per_residue_signed(vec, coprime, modulus):
    """Sum of vec[i] (signed) grouped by coprime[i] mod modulus."""
    out = {}
    for r in range(modulus):
        if r % Q == 0:
            continue
        idx = np.where(coprime % modulus == r)[0]
        if len(idx) == 0:
            continue
        out[r] = float(np.sum(vec[idx]))
    return out


def sparsity_summary(vec, tol=1e-10):
    n = len(vec)
    energy = float(np.sum(vec ** 2))
    if energy == 0:
        return dict(n=n, n_nonzero=0, l_inf=0.0, l1_over_l2=0.0,
                    energy_top1=0.0, energy_top10=0.0,
                    energy_top_pct1=0.0, energy_top_pct10=0.0)
    abs_v = np.abs(vec)
    nz = int(np.sum(abs_v > tol))
    sorted_sq = np.sort(vec ** 2)[::-1]
    cum = np.cumsum(sorted_sq) / energy
    return dict(
        n=n,
        n_nonzero=nz,
        l_inf=float(np.max(abs_v)),
        l1_over_l2=float(np.sum(abs_v) / math.sqrt(energy)),  # ~ sqrt(eff. support size)
        energy_top1=float(cum[0]),
        energy_top10=float(cum[min(9, n - 1)]),
        energy_top_pct1=float(cum[max(0, n // 100 - 1)]),
        energy_top_pct10=float(cum[max(0, n // 10 - 1)]),
    )


def main():
    print("=" * 78)
    print("Probe A: dominant singular vector inspection")
    print("=" * 78)

    findings = {}

    for k in LEVELS:
        print(f"\n--- k = {k} ---")
        npz_path = os.path.join(OUT_DIR, f"R_k{k}_dominant_vectors.npz")
        if not os.path.exists(npz_path):
            print(f"  [missing] {npz_path}; skipping")
            continue
        data = np.load(npz_path)
        u_W = data["u_top"][:, 0]   # u_1 in W coordinates, length dim_W
        v_V = data["v_top"][0, :]   # v_1 in V coordinates, length n_k
        sigma_1 = float(data["sigma_top"][0])
        print(f"  loaded sigma_1 = {sigma_1:.10f}, "
              f"|u_W| = {len(u_W)}, |v_V| = {len(v_V)}")

        # Need P_W to map u_W back into V_{k+1} natural basis
        T, coprime_k, coprime_kp1 = build_T_lift(Q, k)
        P_W = build_W_basis(T)
        n_k = T.shape[0]
        n_kp1 = T.shape[1]
        coprime_k_arr = np.array(coprime_k)
        coprime_kp1_arr = np.array(coprime_kp1)

        # u_1 in V_{k+1} natural basis: u_W is row-vector, embed via @ P_W
        # which gives length n_{k+1}.  In column-vec form: P_W.T @ u_W.
        u_V = P_W.T @ u_W
        # Sanity check: |u_V|^2 should equal |u_W|^2 (orthonormal embedding)
        u_norm_diff = abs(np.linalg.norm(u_V) - np.linalg.norm(u_W))
        print(f"  embedding |u_V| - |u_W| = {u_norm_diff:.2e}")

        # Stationary at level k+1 for comparison
        K_kp1, _, _ = build_K_dense(Q, k + 1)
        pi_kp1 = stationary_vec(K_kp1)
        K_k, _, _ = build_K_dense(Q, k)
        pi_k = stationary_vec(K_k)

        # Per-residue energy at multiple moduli
        moduli_uV = [3, 9, 27, 81]
        moduli_vV = [3, 9, 27]

        print(f"  u_1 in V_{k+1} natural basis (length {n_kp1}):")
        spar_u = sparsity_summary(u_V)
        print(f"    sparsity: {spar_u['n_nonzero']}/{spar_u['n']} nonzero, "
              f"||u||_inf={spar_u['l_inf']:.4e}, "
              f"top-1% energy={spar_u['energy_top_pct1']:.4f}, "
              f"top-10% energy={spar_u['energy_top_pct10']:.4f}")
        for m in moduli_uV:
            if 3 ** k >= m:  # modulus must be ≤ N_{k+1} = 3^{k+1}
                e = per_residue_energy(u_V, coprime_kp1_arr, m)
                top3 = sorted(e.items(), key=lambda kv: -kv[1])[:6]
                print(f"    energy by r mod {m}: top 6 = "
                      + ", ".join(f"r={r}:{p:.4f}" for r, p in top3))

        print(f"  v_1 in V_{k} natural basis (length {n_k}):")
        spar_v = sparsity_summary(v_V)
        print(f"    sparsity: {spar_v['n_nonzero']}/{spar_v['n']} nonzero, "
              f"||v||_inf={spar_v['l_inf']:.4e}, "
              f"top-10% energy={spar_v['energy_top_pct10']:.4f}")
        for m in moduli_vV:
            if 3 ** (k - 1) >= m:
                e = per_residue_energy(v_V, coprime_k_arr, m)
                top = sorted(e.items(), key=lambda kv: -kv[1])
                print(f"    energy by r mod {m}: "
                      + ", ".join(f"r={r}:{p:.4f}" for r, p in top))

        # Inner products with pi_{k+1} and pi_k (alignment with stationary)
        ip_u_pi = float(u_V @ pi_kp1) / (np.linalg.norm(u_V) * np.linalg.norm(pi_kp1))
        ip_v_pi = float(v_V @ pi_k) / (np.linalg.norm(v_V) * np.linalg.norm(pi_k))
        print(f"  cos(u_1, pi_{k+1}) = {ip_u_pi:+.6f}")
        print(f"  cos(v_1, pi_{k}) = {ip_v_pi:+.6f}")

        # Save per-residue CSVs
        out_u = os.path.join(OUT_DIR, f"probe_A_u1_support_k{k}.csv")
        with open(out_u, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["residue", "modulus", "energy_fraction", "signed_sum"])
            for m in moduli_uV:
                if 3 ** k >= m:
                    e = per_residue_energy(u_V, coprime_kp1_arr, m)
                    s = per_residue_signed(u_V, coprime_kp1_arr, m)
                    for r in sorted(e.keys()):
                        w.writerow([r, m, f"{e[r]:.10e}", f"{s[r]:+.10e}"])

        out_v = os.path.join(OUT_DIR, f"probe_A_v1_support_k{k}.csv")
        with open(out_v, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["residue", "modulus", "energy_fraction", "signed_sum"])
            for m in moduli_vV:
                if 3 ** (k - 1) >= m:
                    e = per_residue_energy(v_V, coprime_k_arr, m)
                    s = per_residue_signed(v_V, coprime_k_arr, m)
                    for r in sorted(e.keys()):
                        w.writerow([r, m, f"{e[r]:.10e}", f"{s[r]:+.10e}"])

        findings[k] = {
            "sigma_1": sigma_1,
            "n_k": n_k, "n_kp1": n_kp1,
            "spar_u": spar_u, "spar_v": spar_v,
            "ip_u_pi": ip_u_pi, "ip_v_pi": ip_v_pi,
            "u_V": u_V, "v_V": v_V,
            "coprime_k": coprime_k_arr, "coprime_kp1": coprime_kp1_arr,
        }

    # Cross-k consistency: if u_1 has the same residue-class structure across
    # k, that's evidence the dominant singular direction is structurally fixed.
    print("\n=== Cross-k consistency of u_1's mod-3 structure ===")
    for k in sorted(findings.keys()):
        f = findings[k]
        e3 = per_residue_energy(f["u_V"], f["coprime_kp1"], 3)
        e9 = per_residue_energy(f["u_V"], f["coprime_kp1"], 9)
        print(f"  k={k}: u_1 mod 3 = " +
              ", ".join(f"r={r}:{p:.4f}" for r, p in sorted(e3.items())) +
              f"  | mod 9 top: " +
              ", ".join(f"r={r}:{p:.3f}" for r, p in
                       sorted(e9.items(), key=lambda kv: -kv[1])[:4]))

    # ---- Markdown writeup ----
    md_path = os.path.join(OUT_DIR, "probe_A_findings.md")
    lines = []
    lines.append("# Probe A — dominant singular vector support inspection")
    lines.append("")
    lines.append("Inspecting u_1 (left singular vector in W_{k+1} embedded back "
                 "into V_{k+1} natural basis) and v_1 (right singular vector "
                 "in V_k natural basis) for R_k = P_{W_{k+1}} ∘ K_{k+1} ∘ L_k "
                 "at k = 5, 6, 7.")
    lines.append("")
    lines.append("## Per-k summary")
    lines.append("")
    lines.append("| k | sigma_1 | u_1 nz/n_{k+1} | v_1 nz/n_k | top-10% E(u_1) "
                 "| top-10% E(v_1) | cos(u_1, pi_{k+1}) | cos(v_1, pi_k) |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for k in sorted(findings.keys()):
        f = findings[k]
        u_s, v_s = f["spar_u"], f["spar_v"]
        lines.append(f"| {k} | {f['sigma_1']:.6f} | "
                     f"{u_s['n_nonzero']}/{u_s['n']} | "
                     f"{v_s['n_nonzero']}/{v_s['n']} | "
                     f"{u_s['energy_top_pct10']:.4f} | "
                     f"{v_s['energy_top_pct10']:.4f} | "
                     f"{f['ip_u_pi']:+.4f} | {f['ip_v_pi']:+.4f} |")
    lines.append("")

    lines.append("## u_1 mod-3 residue structure across k")
    lines.append("")
    lines.append("Energy fraction of u_1 in V_{k+1} on residues r mod 3 ∈ {1, 2}:")
    lines.append("")
    lines.append("| k | r=1 mod 3 | r=2 mod 3 |")
    lines.append("|---|---|---|")
    for k in sorted(findings.keys()):
        f = findings[k]
        e3 = per_residue_energy(f["u_V"], f["coprime_kp1"], 3)
        lines.append(f"| {k} | {e3.get(1, 0):.6f} | {e3.get(2, 0):.6f} |")
    lines.append("")

    lines.append("## u_1 mod-9 residue structure")
    lines.append("")
    lines.append("Energy fraction of u_1 grouped by r mod 9 (only coprime "
                 "residues r ∈ {1,2,4,5,7,8}):")
    lines.append("")
    lines.append("| k | r=1 | r=2 | r=4 | r=5 | r=7 | r=8 |")
    lines.append("|---|---|---|---|---|---|---|")
    for k in sorted(findings.keys()):
        f = findings[k]
        e9 = per_residue_energy(f["u_V"], f["coprime_kp1"], 9)
        lines.append(f"| {k} | "
                     f"{e9.get(1, 0):.6f} | {e9.get(2, 0):.6f} | "
                     f"{e9.get(4, 0):.6f} | {e9.get(5, 0):.6f} | "
                     f"{e9.get(7, 0):.6f} | {e9.get(8, 0):.6f} |")
    lines.append("")

    lines.append("## v_1 mod-3 residue structure (input direction in V_k)")
    lines.append("")
    lines.append("| k | r=1 mod 3 | r=2 mod 3 |")
    lines.append("|---|---|---|")
    for k in sorted(findings.keys()):
        f = findings[k]
        e3v = per_residue_energy(f["v_V"], f["coprime_k"], 3)
        lines.append(f"| {k} | {e3v.get(1, 0):.6f} | {e3v.get(2, 0):.6f} |")
    lines.append("")

    lines.append("## Stationary alignment")
    lines.append("")
    lines.append("cos(u_1, pi_{k+1}) and cos(v_1, pi_k): if either is near ±1, "
                 "the dominant SVD direction is the stationary itself; if near "
                 "0, it's orthogonal (a non-trivial deviation mode).")
    lines.append("")
    for k in sorted(findings.keys()):
        f = findings[k]
        lines.append(f"- k={k}: cos(u_1, pi_{k+1}) = {f['ip_u_pi']:+.6f}, "
                     f"cos(v_1, pi_k) = {f['ip_v_pi']:+.6f}")
    lines.append("")

    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"\n[md: {md_path}]")
    print("Done.")


if __name__ == "__main__":
    main()
