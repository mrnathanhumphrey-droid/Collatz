"""
qsd_late_t_analysis.py — Diagnostic: at late t (130-190), does D(r,t)
stabilize? If so, what is the asymptotic shape, and does it match any QSD?

Combines extended-horizon empirical data with QSD eigenvectors.
"""
import sys
import numpy as np
from fractions import Fraction
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


def build_chang_kernel():
    odd_residues = list(range(1, 64, 2))
    idx = {r: i for i, r in enumerate(odd_residues)}
    P = [[Fraction(0)] * 32 for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(128):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm & 1 == 0:
                mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], 128)
    return P, odd_residues, idx


def stationary_exact(P_frac):
    n = 32
    A = [[P_frac[j][i] - (Fraction(1) if i == j else Fraction(0))
          for j in range(n)] for i in range(n)]
    A[n - 1] = [Fraction(1)] * n
    b = [Fraction(0)] * n
    b[n - 1] = Fraction(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if A[row][col] != 0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            b[col], b[pivot] = b[pivot], b[col]
        piv = A[col][col]
        for j in range(col, n):
            A[col][j] /= piv
        b[col] /= piv
        for row in range(n):
            if row != col and A[row][col] != 0:
                factor = A[row][col]
                for j in range(col, n):
                    A[row][j] -= factor * A[col][j]
                b[row] -= factor * b[col]
    return b


def build_qsd(P_frac, odd_residues, idx, absorbing_set):
    surviving = [r for r in odd_residues if r not in absorbing_set]
    n_sub = len(surviving)
    sub_idx = {r: i for i, r in enumerate(surviving)}
    P_sub = np.zeros((n_sub, n_sub))
    for r_from in surviving:
        for r_to in surviving:
            P_sub[sub_idx[r_from], sub_idx[r_to]] = float(
                P_frac[idx[r_from]][idx[r_to]]
            )
    eigvals, eigvecs = np.linalg.eig(P_sub.T)
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    lambda_PF = eigvals[0].real
    v = eigvecs[:, 0].real
    if v.sum() < 0:
        v = -v
    v = v / v.sum()
    v_dict = {r: v[sub_idx[r]] for r in surviving}
    return lambda_PF, v_dict, eigvals, surviving


def project_to_mod32(dist_dict_64):
    out = {}
    for r32 in range(1, 32, 2):
        out[r32] = dist_dict_64.get(r32, 0.0) + dist_dict_64.get(r32 + 32, 0.0)
    return out


def main():
    log("=" * 78)
    log("QSD late-t analysis: empirical asymptote vs QSD eigenvectors")
    log("=" * 78)

    # Load extended-horizon empirical data
    rho_data = {}  # t -> {r: rho}
    n_alive = {}
    odd_r32 = list(range(1, 32, 2))
    with open(OUT / "qsd_extended.csv") as f:
        header = f.readline().strip().split(",")
        for line in f:
            parts = line.strip().split(",")
            t = int(parts[0])
            n_alive[t] = int(parts[1])
            rho_data[t] = {r: float(parts[2 + i]) for i, r in enumerate(odd_r32)}

    # Build kernel
    P_frac, odd_residues, idx = build_chang_kernel()
    pi = stationary_exact(P_frac)
    pi_64 = {r: float(pi[idx[r]]) for r in odd_residues}
    pi_32 = project_to_mod32(pi_64)

    # Compute D(r, t) = rho(r, t) / pi(r)
    D_data = {}
    for t in rho_data:
        D_data[t] = {r: rho_data[t][r] / pi_32[r] for r in odd_r32}

    # Average late-t D (t=130 onwards) to reduce noise
    late_t = [t for t in sorted(D_data.keys()) if t >= 130]
    log(f"\nLate-t snapshots used for averaging: {late_t}")
    log(f"  n_alive at these snapshots: " +
        ", ".join(f"t={t}: {n_alive[t]}" for t in late_t))

    # Weight by n_alive (sampling-noise weighting)
    D_avg = {}
    D_var = {}
    for r in odd_r32:
        # weighted mean by n_alive (each rho has variance ~ rho/n)
        weights = [n_alive[t] for t in late_t]
        vals = [D_data[t][r] for t in late_t]
        D_avg[r] = np.average(vals, weights=weights)
        # Sample variance across snapshots
        D_var[r] = np.var(vals, ddof=1)

    log(f"\nAveraged D(r) over late-t snapshots (weighted by n_alive):")
    log(f"  {'r':>3}  {'D_avg':>7}  {'std':>7}  {'individual snapshots':>40}")
    for r in odd_r32:
        snap_str = " ".join(f"{D_data[t][r]:.2f}" for t in late_t)
        log(f"  {r:>3}  {D_avg[r]:>7.4f}  {np.sqrt(D_var[r]):>7.4f}  {snap_str}")

    # Compute QSDs for relevant absorption sets
    log("\n" + "=" * 78)
    log("Compare D_avg(r) to QSD predictions")
    log("=" * 78)

    absorption_sets = {
        "(a) {21}": [21],
        "(b) {21, 53}": [21, 53],
        "(c) {1}": [1],
        "(d) {1, 33}": [1, 33],
        "(e) {5, 37}": [5, 37],
        "(b+e) {5,21,37,53}": [5, 21, 37, 53],
        "(b+c) {1,21,53}": [1, 21, 53],
        "(b+e+c) {1,5,21,37,53}": [1, 5, 21, 37, 53],
        "all m_j cyl {1,5,21,33,37,53}": [1, 5, 21, 33, 37, 53],
    }

    log(f"\n{'absorption':<32}  {'lambda_PF':>10}  {'sum |D_avg-DQSD|':>17}  "
        f"{'spectral gap':>13}")
    qsd_results = {}
    for label, A in absorption_sets.items():
        lambda_PF, v_dict_64, eigvals, surviving = build_qsd(
            P_frac, odd_residues, idx, A
        )
        v_32 = project_to_mod32(v_dict_64)
        # D_QSD at mod-32 level: only for residues whose BOTH lifts survive
        D_qsd = {}
        for r32 in odd_r32:
            if v_32[r32] > 0 and pi_32[r32] > 0:
                D_qsd[r32] = v_32[r32] / pi_32[r32]
        # Compare to D_avg
        common = [r for r in odd_r32 if r in D_qsd]
        total_dev = sum(abs(D_avg[r] - D_qsd[r]) for r in common)
        spec_gap = abs(eigvals[1]) / lambda_PF
        qsd_results[label] = {
            "lambda_PF": lambda_PF,
            "v_32": v_32,
            "D_qsd": D_qsd,
            "total_dev": total_dev,
            "spec_gap": spec_gap,
        }
        log(f"{label:<32}  {lambda_PF:>10.6f}  {total_dev:>17.4f}  "
            f"{spec_gap:>13.4f}")

    # Find best
    best = min(qsd_results.keys(), key=lambda k: qsd_results[k]["total_dev"])
    log(f"\nBest match: {best} with total dev {qsd_results[best]['total_dev']:.4f}")

    # Per-residue comparison for best
    log(f"\nPer-residue comparison for {best}:")
    log(f"  {'r':>3}  {'D_avg':>7}  {'D_QSD':>7}  {'diff':>8}  "
        f"{'D_avg/D_QSD':>11}")
    for r in odd_r32:
        d_avg = D_avg[r]
        d_qsd = qsd_results[best]["D_qsd"].get(r)
        if d_qsd is not None:
            diff = d_avg - d_qsd
            ratio = d_avg / d_qsd
            log(f"  {r:>3}  {d_avg:>7.4f}  {d_qsd:>7.4f}  {diff:>+8.4f}  "
                f"{ratio:>11.4f}")
        else:
            log(f"  {r:>3}  {d_avg:>7.4f}  {'(absorb)':>7}  {'—':>8}  {'—':>11}")

    # Empirical asymptotic survival rate
    log("\n" + "=" * 78)
    log("Empirical per-step survival rate at late t")
    log("=" * 78)
    log("\nFrom qsd_extended.csv n_alive:")
    log(f"  {'t-range':>10}  {'survival rate (per step)':>26}")
    sorted_t = sorted(rho_data.keys())
    for i in range(1, len(sorted_t)):
        t0, t1 = sorted_t[i-1], sorted_t[i]
        if n_alive[t0] > 0:
            r = (n_alive[t1] / n_alive[t0]) ** (1 / (t1 - t0))
            log(f"  {t0:>3}-{t1:>3}    {r:>20.6f}")

    # Save CSVs
    with open(OUT / "qsd_late_t_avg.csv", "w") as f:
        f.write("r,pi,D_avg,D_std," +
                ",".join(f"D_QSD_{label.replace(',','_').replace(' ','')}"
                         for label in qsd_results.keys()) + "\n")
        for r in odd_r32:
            row = [str(r), f"{pi_32[r]:.6f}",
                   f"{D_avg[r]:.6f}", f"{np.sqrt(D_var[r]):.6f}"]
            for label in qsd_results.keys():
                d_qsd = qsd_results[label]["D_qsd"].get(r)
                row.append(f"{d_qsd:.6f}" if d_qsd is not None else "")
            f.write(",".join(row) + "\n")

    # Save log
    (OUT / "qsd_late_t_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n[wrote] qsd_late_t_avg.csv, qsd_late_t_log.txt")


if __name__ == "__main__":
    main()
