"""
qsd_computation.py — Result 49: Formal QSD on Chang's mod-64 kernel
under {m_j} absorption.

Tests three absorption sets:
  (a) {21} mod 64               — single m_j chain entry
  (b) {21, 53} mod 64           — full r=21 mod 32 cylinder (both lifts)
  (c) {1} mod 64                — trivial cycle (matches Door 3 empirical)

For each:
  - Build P_sub (substochastic restriction)
  - Eigendecompose -> Perron eigenvalue lambda_PF + QSD eigenvector v
  - Compute spectral gap |lambda_2 / lambda_PF|
  - Compute D_QSD(r) = v(r) / pi(r)
  - Compare to empirical D(r, t) at t = 50, 70, 90, 110

Survival-rate diagnostic: lambda_PF^110 should ~= 0.12 if the absorption
set matches what defines the empirical survivor population.

Functional fits to empirical D(r, t):
  - exp:  D = D_inf + A exp(-t/tau)
  - power: D = D_inf + A t^(-alpha)
  - log:  D = a + b log(t)
  - lin:  D = a + b t

Verdict gates:
  alpha: empirical D(r,t) -> v(r)/pi(r) (one of the three QSDs)
  beta:  D(r,t) converges but not to QSD prediction
  gamma: D(r,t) does not converge (logarithmic / linear drift)
"""
import sys
import numpy as np
from fractions import Fraction
from pathlib import Path
from scipy.optimize import curve_fit

sys.stdout.reconfigure(encoding="utf-8")

OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


# ----------------------------------------------------------------------
# Build Chang's cylinder-averaged kernel (Definition C.5) exactly
# ----------------------------------------------------------------------
def build_chang_kernel():
    odd_residues = list(range(1, 64, 2))  # 32 odd residues mod 64
    idx = {r: i for i, r in enumerate(odd_residues)}
    N_LIFTS = 128  # 2^7 lifts -> mod 64 + 7 = mod 8192

    P = [[Fraction(0)] * 32 for _ in range(32)]
    for i, r in enumerate(odd_residues):
        counts = [0] * 32
        for k in range(N_LIFTS):
            n = r + 64 * k
            mm = 3 * n + 1
            while mm & 1 == 0:
                mm >>= 1
            counts[idx[mm % 64]] += 1
        for j in range(32):
            P[i][j] = Fraction(counts[j], N_LIFTS)
    return P, odd_residues, idx


def stationary_exact(P_frac):
    """Solve pi P = pi exactly. P is row-stochastic."""
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


# ----------------------------------------------------------------------
# QSD: build P_sub by deleting absorbing rows/cols
# ----------------------------------------------------------------------
def build_qsd(P_frac, odd_residues, idx, absorbing_set):
    """
    Construct P_sub (substochastic) by removing absorbing residues.
    Compute Perron eigenvalue and right (Yaglom-limit) eigenvector.
    Returns (lambda_PF, v_dict, lambda_2_over_PF, surviving_residues, P_sub_float).
    """
    surviving = [r for r in odd_residues if r not in absorbing_set]
    n_sub = len(surviving)
    sub_idx = {r: i for i, r in enumerate(surviving)}

    # Float matrix for eigendecomp (32x32 over rationals -> float OK at this size)
    P_sub = np.zeros((n_sub, n_sub))
    for r_from in surviving:
        for r_to in surviving:
            P_sub[sub_idx[r_from], sub_idx[r_to]] = float(
                P_frac[idx[r_from]][idx[r_to]]
            )

    # Eigendecomposition
    # Yaglom limit / QSD = LEFT Perron eigenvector of P_sub (interpretable as
    # row vector v with v P_sub = lambda v).  Equivalently, right Perron of P_sub.T.
    eigvals, eigvecs = np.linalg.eig(P_sub.T)
    # Sort by |eigvalue| descending
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    lambda_PF = eigvals[0].real  # should be real positive
    lambda_2 = eigvals[1]
    v = eigvecs[:, 0].real

    # Normalize: should be all same sign; flip if needed
    if v.sum() < 0:
        v = -v
    v = v / v.sum()  # probability normalized

    # Convert to dict
    v_dict = {r: v[sub_idx[r]] for r in surviving}

    # Spectral gap
    spectral_gap_ratio = abs(lambda_2) / lambda_PF

    return lambda_PF, v_dict, spectral_gap_ratio, surviving, P_sub


# ----------------------------------------------------------------------
# Project mod-64 distribution to mod-32 (sums {r, r+32})
# ----------------------------------------------------------------------
def project_to_mod32(dist_dict_64, odd_residues_64):
    """Sum dist[r] + dist[r+32] for each r in 1..31 odd."""
    out = {}
    for r32 in range(1, 32, 2):
        out[r32] = dist_dict_64.get(r32, 0.0) + dist_dict_64.get(r32 + 32, 0.0)
    return out


# ----------------------------------------------------------------------
# Empirical D(r, t) from Door 3 CSV
# ----------------------------------------------------------------------
def load_empirical_D():
    """Returns dict t -> {r: D(r,t)} and dict t -> n_alive."""
    rows = OUT / "experiments_output" / "chang_qsd_test.csv"
    D_emp = {}
    n_alive = {}
    with open(rows) as f:
        header = f.readline().strip().split(",")
        # Format: t, n_alive, D_r1, D_r3, ..., D_r31
        residues_h = [int(h.replace("D_r", "")) for h in header[2:]]
        for line in f:
            parts = line.strip().split(",")
            t = int(parts[0])
            n_alive[t] = int(parts[1])
            D_emp[t] = {r: float(parts[2 + i]) for i, r in enumerate(residues_h)}
    return D_emp, n_alive


# ----------------------------------------------------------------------
# Functional fits to D(r, t)
# ----------------------------------------------------------------------
def fit_drift_models(t_arr, D_arr):
    """
    Fit four models. Return dict with R^2, params, asymptote (or None).
    """
    out = {}
    t = np.asarray(t_arr, dtype=float)
    D = np.asarray(D_arr, dtype=float)

    # exponential approach: D = D_inf + A * exp(-t / tau)
    try:
        def fexp(t, D_inf, A, tau):
            return D_inf + A * np.exp(-t / tau)
        # Initial guesses
        p0 = [D[-1], D[0] - D[-1], 50.0]
        popt, _ = curve_fit(fexp, t, D, p0=p0, maxfev=10000)
        D_pred = fexp(t, *popt)
        ss_res = np.sum((D - D_pred) ** 2)
        ss_tot = np.sum((D - D.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        out["exp"] = {"params": popt, "asymptote": popt[0], "r2": r2}
    except Exception as e:
        out["exp"] = {"params": None, "asymptote": None, "r2": -np.inf}

    # power-law approach: D = D_inf + A * t^(-alpha) (skip t=0)
    try:
        mask = t > 0
        def fpow(tt, D_inf, A, alpha):
            return D_inf + A * tt ** (-alpha)
        p0 = [D[-1], D[0] - D[-1], 0.5]
        popt, _ = curve_fit(fpow, t[mask], D[mask], p0=p0, maxfev=10000)
        D_pred = fpow(t[mask], *popt)
        ss_res = np.sum((D[mask] - D_pred) ** 2)
        ss_tot = np.sum((D[mask] - D[mask].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        out["pow"] = {"params": popt, "asymptote": popt[0], "r2": r2}
    except Exception:
        out["pow"] = {"params": None, "asymptote": None, "r2": -np.inf}

    # logarithmic drift: D = a + b * log(t) (no asymptote)
    try:
        mask = t > 0
        slope, intercept = np.polyfit(np.log(t[mask]), D[mask], 1)
        D_pred = intercept + slope * np.log(t[mask])
        ss_res = np.sum((D[mask] - D_pred) ** 2)
        ss_tot = np.sum((D[mask] - D[mask].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        out["log"] = {"params": (intercept, slope), "asymptote": None, "r2": r2}
    except Exception:
        out["log"] = {"params": None, "asymptote": None, "r2": -np.inf}

    # linear drift
    try:
        slope, intercept = np.polyfit(t, D, 1)
        D_pred = intercept + slope * t
        ss_res = np.sum((D - D_pred) ** 2)
        ss_tot = np.sum((D - D.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        out["lin"] = {"params": (intercept, slope), "asymptote": None, "r2": r2}
    except Exception:
        out["lin"] = {"params": None, "asymptote": None, "r2": -np.inf}

    # Pick best
    best = max(out.keys(), key=lambda k: out[k]["r2"])
    out["best"] = best
    return out


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    log("=" * 78)
    log("QSD computation on Chang cylinder-averaged kernel mod 64")
    log("=" * 78)

    log("\n[Step 1] Build Chang kernel P (Definition C.5) and stationary pi...")
    P_frac, odd_residues, idx = build_chang_kernel()

    # Verify row sums
    for i in range(32):
        assert sum(P_frac[i]) == 1
    log("  Row sums OK")

    pi_64 = stationary_exact(P_frac)
    pi_dict_64 = {r: float(pi_64[idx[r]]) for r in odd_residues}
    pi_mod32 = project_to_mod32(pi_dict_64, odd_residues)

    # I_2 verification
    I2 = [7, 27, 31, 59, 63]
    mass_I2 = sum(pi_64[idx[r]] for r in I2)
    log(f"  pi(I_2) = {mass_I2} = {float(mass_I2):.7f}")
    log(f"  Chang reports 10121/65280 = {float(Fraction(10121, 65280)):.7f}")
    log(f"  Exact match: {mass_I2 == Fraction(10121, 65280)}")

    log("\n  pi mod-32 (uniform = 1/16 = 0.0625):")
    for r32 in range(1, 32, 2):
        log(f"    r={r32:>2}: pi = {pi_mod32[r32]:.6f}, "
            f"pi/unif = {pi_mod32[r32] * 16:.4f}")

    # ------------------------------------------------------------------
    # Step 2: Build P_sub for three absorption choices
    # ------------------------------------------------------------------
    log("\n[Step 2] QSD eigendecomposition for three absorption sets")

    absorption_sets = {
        "(a) {21}": [21],
        "(b) {21, 53}": [21, 53],
        "(c) {1}": [1],
    }

    # Bonus extras for sensitivity
    absorption_sets["(d) {1, 33}"] = [1, 33]  # full mod-32 r=1 cylinder
    absorption_sets["(e) {5, 37}"] = [5, 37]  # m_2 = 5 cylinder
    absorption_sets["(f) {1, 5, 21}"] = [1, 5, 21]  # all three small m_j
    # Full Chang invariant core I_2 (control)
    absorption_sets["(g) I_2={7,27,31,59,63}"] = [7, 27, 31, 59, 63]

    qsd_results = {}
    for label, A in absorption_sets.items():
        lambda_PF, v_dict_64, spec_gap, surviving, _ = build_qsd(
            P_frac, odd_residues, idx, A
        )
        v_mod32 = project_to_mod32(v_dict_64, odd_residues)
        qsd_results[label] = {
            "absorbing": A,
            "lambda_PF": lambda_PF,
            "spec_gap_ratio": spec_gap,
            "v_64": v_dict_64,
            "v_32": v_mod32,
        }
        log(f"\n  {label}:")
        log(f"    lambda_PF = {lambda_PF:.6f}")
        log(f"    |lambda_2|/|lambda_PF| = {spec_gap:.4f}")
        log(f"    survival^110 = lambda_PF^110 = {lambda_PF ** 110:.4f}")

    # ------------------------------------------------------------------
    # Step 7 (early): survival-rate diagnostic
    # ------------------------------------------------------------------
    log("\n[Step 7 (early)] Survival-rate diagnostic")
    log("  Empirical: 236389/2000000 = 0.1182 survivors at t=110")
    log("  Required: lambda_PF^110 = 0.1182 -> lambda_PF = 0.1182^(1/110)")
    target_lambda = 0.1182 ** (1 / 110)
    log(f"  Target lambda_PF = {target_lambda:.6f}")
    log("  Comparing computed vs target:")
    for label, r in qsd_results.items():
        diff = r["lambda_PF"] - target_lambda
        match = "MATCH" if abs(diff) < 0.005 else ""
        log(f"    {label}: lambda_PF={r['lambda_PF']:.6f}, "
            f"diff={diff:+.4f} {match}")

    # ------------------------------------------------------------------
    # Step 3: Compare D_QSD to empirical D(r, t)
    # ------------------------------------------------------------------
    log("\n[Step 3] D_QSD(r) = v(r) / pi(r) vs empirical D(r, t)")

    D_emp, n_alive = load_empirical_D()
    log(f"  Loaded empirical D(r, t) for t in {sorted(D_emp.keys())}")
    log(f"  n_alive at each t: " +
        ", ".join(f"t={t}:{n}" for t, n in sorted(n_alive.items())))

    # Compute D_QSD for each absorption choice
    for label, r in qsd_results.items():
        log(f"\n  --- {label} ---")
        v32 = r["v_32"]
        # pi32 may have value 0 at absorbing residues; but absorbing residues
        # were in the v_dict so they map cleanly.  For projection to mod 32,
        # absorbing residues at mod-64 level might leave the mod-32 with only
        # half the mass.  We need to compute D = v(r)/pi(r) only for r where
        # both are nonzero.
        D_qsd = {}
        for r32 in range(1, 32, 2):
            if v32[r32] > 0 and pi_mod32[r32] > 0:
                # Renormalize v32 to be a probability over surviving mod-32
                # residues — but only project absorbing residues in mod-32
                # if BOTH their lifts are absorbed
                D_qsd[r32] = v32[r32] / pi_mod32[r32]
            else:
                D_qsd[r32] = None

        # Header
        line_h = f"    {'r':>3}  {'D_QSD':>8}"
        for t in [50, 70, 90, 110]:
            line_h += f"  {'D(t='+str(t)+')':>10}"
        line_h += f"  {'|D110-DQSD|':>11}"
        log(line_h)
        for r32 in range(1, 32, 2):
            d_qsd = D_qsd[r32]
            d_qsd_str = f"{d_qsd:.4f}" if d_qsd is not None else "  —  "
            line = f"    {r32:>3}  {d_qsd_str:>8}"
            for t in [50, 70, 90, 110]:
                line += f"  {D_emp[t].get(r32, float('nan')):>10.4f}"
            if d_qsd is not None:
                line += f"  {abs(D_emp[110][r32] - d_qsd):>11.4f}"
            else:
                line += f"  {'—':>11}"
            log(line)

        # Total deviation between D_QSD and D_empirical(t=110)
        residues_with_qsd = [r32 for r32 in range(1, 32, 2)
                             if D_qsd[r32] is not None]
        total_dev = sum(
            abs(D_emp[110][r32] - D_qsd[r32])
            for r32 in residues_with_qsd
        )
        log(f"    Total |D(t=110) - D_QSD|: {total_dev:.4f} "
            f"over {len(residues_with_qsd)} non-absorbed residues")

    # ------------------------------------------------------------------
    # Step 4: Functional fits to D(r, t)
    # ------------------------------------------------------------------
    log("\n[Step 4] Functional fits to empirical D(r, t)")

    t_vals = sorted(D_emp.keys())
    drift_fits = {}
    log(f"\n  {'r':>3}  {'best':>5}  {'r2_exp':>7}  {'r2_pow':>7}  "
        f"{'r2_log':>7}  {'r2_lin':>7}  {'D_inf_exp':>10}  "
        f"{'D_inf_pow':>10}")
    for r32 in range(1, 32, 2):
        D_arr = [D_emp[t][r32] for t in t_vals]
        fits = fit_drift_models(t_vals, D_arr)
        drift_fits[r32] = fits
        log(f"  {r32:>3}  {fits['best']:>5}  "
            f"{fits['exp']['r2']:>7.4f}  {fits['pow']['r2']:>7.4f}  "
            f"{fits['log']['r2']:>7.4f}  {fits['lin']['r2']:>7.4f}  "
            f"{fits['exp']['asymptote']:>10.4f}  "
            f"{fits['pow']['asymptote']:>10.4f}")

    # Aggregate best-fit family
    fit_counts = {"exp": 0, "pow": 0, "log": 0, "lin": 0}
    for r32, f in drift_fits.items():
        fit_counts[f["best"]] += 1
    log(f"\n  Best-fit family count: {fit_counts}")

    # ------------------------------------------------------------------
    # Step 6: Absorption sensitivity - find best-matching QSD
    # ------------------------------------------------------------------
    log("\n[Step 6] Absorption-set sensitivity: which QSD best matches "
        "empirical D(t=110)?")
    log(f"\n  {'absorption':<28}  {'lambda_PF':>10}  "
        f"{'lam^110':>9}  {'sum |D-DQSD|':>13}")
    best_label = None
    best_total_dev = float("inf")
    for label, r in qsd_results.items():
        D_qsd = {}
        for r32 in range(1, 32, 2):
            v32 = r["v_32"]
            if v32[r32] > 0 and pi_mod32[r32] > 0:
                D_qsd[r32] = v32[r32] / pi_mod32[r32]
        residues_with_qsd = list(D_qsd.keys())
        total_dev = sum(
            abs(D_emp[110][r32] - D_qsd[r32])
            for r32 in residues_with_qsd
        )
        log(f"  {label:<28}  {r['lambda_PF']:>10.6f}  "
            f"{r['lambda_PF'] ** 110:>9.4f}  {total_dev:>13.4f}")
        if total_dev < best_total_dev:
            best_total_dev = total_dev
            best_label = label
    log(f"\n  BEST MATCH: {best_label} (total dev = {best_total_dev:.4f})")

    # ------------------------------------------------------------------
    # Save CSVs
    # ------------------------------------------------------------------
    log("\n[Save] CSVs")

    # qsd_eigendecomp.csv
    with open(OUT / "qsd_eigendecomp.csv", "w") as f:
        f.write("absorption,r,v_64,v_32,lambda_PF,spec_gap_ratio\n")
        for label, r in qsd_results.items():
            for r64 in odd_residues:
                v64 = r["v_64"].get(r64, 0.0)
                f.write(f"{label},{r64},{v64:.10f},,"
                        f"{r['lambda_PF']:.10f},"
                        f"{r['spec_gap_ratio']:.10f}\n")
            for r32 in range(1, 32, 2):
                v32 = r["v_32"][r32]
                f.write(f"{label},{r32}_mod32,,{v32:.10f},,\n")

    # qsd_vs_empirical.csv
    with open(OUT / "qsd_vs_empirical.csv", "w") as f:
        cols = ["r", "pi_mod32", "D_emp_t50", "D_emp_t70",
                "D_emp_t90", "D_emp_t110"]
        for label in qsd_results.keys():
            cols.append(f"D_QSD_{label}")
        f.write(",".join(cols) + "\n")
        for r32 in range(1, 32, 2):
            row = [str(r32), f"{pi_mod32[r32]:.6f}",
                   f"{D_emp[50][r32]:.6f}",
                   f"{D_emp[70][r32]:.6f}",
                   f"{D_emp[90][r32]:.6f}",
                   f"{D_emp[110][r32]:.6f}"]
            for label, r in qsd_results.items():
                v32 = r["v_32"][r32]
                if v32 > 0 and pi_mod32[r32] > 0:
                    row.append(f"{v32 / pi_mod32[r32]:.6f}")
                else:
                    row.append("")
            f.write(",".join(row) + "\n")

    # qsd_drift_fits.csv
    with open(OUT / "qsd_drift_fits.csv", "w") as f:
        f.write("r,best_model,r2_exp,r2_pow,r2_log,r2_lin,"
                "D_inf_exp,A_exp,tau_exp,D_inf_pow,A_pow,alpha_pow\n")
        for r32 in range(1, 32, 2):
            fits = drift_fits[r32]
            row = [str(r32), fits["best"]]
            for k in ["exp", "pow", "log", "lin"]:
                row.append(f"{fits[k]['r2']:.6f}")
            # Exp params
            if fits["exp"]["params"] is not None:
                row.extend(f"{p:.6f}" for p in fits["exp"]["params"])
            else:
                row.extend(["", "", ""])
            # Pow params
            if fits["pow"]["params"] is not None:
                row.extend(f"{p:.6f}" for p in fits["pow"]["params"])
            else:
                row.extend(["", "", ""])
            f.write(",".join(row) + "\n")

    log("[wrote] qsd_eigendecomp.csv, qsd_vs_empirical.csv, "
        "qsd_drift_fits.csv")

    (OUT / "qsd_computation_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("[wrote] qsd_computation_log.txt")


if __name__ == "__main__":
    main()
