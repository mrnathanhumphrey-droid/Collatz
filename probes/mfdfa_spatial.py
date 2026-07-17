"""
mfdfa_spatial.py — proper multifractal analysis of the trajectory measure
as a spatial measure on Z_2 (NOT a temporal time series).

Result 59's mass-dim_q2 sweep showed multifractality lives in the SPATIAL
distribution of mass across cylinders mod 2^k, not in temporal v_t/log m_t
along orbits (which are essentially Brownian, MF-DFA gives h~0.5 for q>=0).

This script extends Result 59 to the FULL multifractal spectrum:
  Z_q(k) = Sigma_{cylinders C mod 2^k} mu(C)^q ~ 2^(-k * tau(q))
  D_q    = tau(q) / (q - 1)              Renyi dimension
  alpha  = d tau / d q                   singularity strength
  f(alpha) = q * alpha - tau(q)         singularity spectrum

For mu, use Result 58's inverse-tree subtree-size weighting (Pearson +0.86
with D_emp; established as the best single-weight family).
"""
import math
import sys
import time
from collections import deque, defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)


# ----------------------------------------------------------------
# Inverse tree (from zadic_measure_framework / Result 52)
# ----------------------------------------------------------------
def build_inverse_tree(max_value):
    tree = {1: {'parent': None, 'depth': 0}}
    q = deque([1])
    while q:
        m = q.popleft()
        d = tree[m]['depth']
        if m % 3 == 0:
            continue
        v_start = 2 if (m % 3 == 1) else 1
        for v in range(v_start, 64, 2):
            num = m * (1 << v) - 1
            if num <= 0:
                continue
            pred = num // 3
            if pred > max_value:
                break
            if pred & 1 == 0 or pred == m:
                continue
            if pred not in tree:
                tree[pred] = {'parent': m, 'depth': d + 1}
                q.append(pred)
    return tree


def compute_subtree_sizes(tree):
    """Subtree size = number of descendants in tree (including self)."""
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    # Post-order DFS
    order = []
    stack = [1]
    visited = set()
    while stack:
        m = stack[-1]
        if m in visited:
            stack.pop()
            order.append(m)
            continue
        visited.add(m)
        for c in children.get(m, []):
            stack.append(c)
    sizes = {}
    for m in order:
        sizes[m] = 1 + sum(sizes[c] for c in children.get(m, []))
    return sizes


# ----------------------------------------------------------------
# Multifractal partition function Z_q(2^k)
# ----------------------------------------------------------------
def partition_function_Zq(weights_by_m, ks, qs):
    """
    For each modulus M = 2^k, group odd m by residue, compute mass
    mu(C) normalized to sum 1. Then Z_q(k) = Sigma mu(C)^q.

    Returns dict: Z[q][k_idx]
    """
    odd_m = [m for m in weights_by_m if m & 1]
    Z = {q: np.zeros(len(ks)) for q in qs}
    n_cyls = np.zeros(len(ks), dtype=np.int64)
    for ki, k in enumerate(ks):
        mod = 1 << k
        cyls = defaultdict(float)
        for m in odd_m:
            cyls[m % mod] += weights_by_m[m]
        total = sum(cyls.values())
        if total <= 0:
            continue
        masses = np.array([c / total for c in cyls.values()])
        n_cyls[ki] = len(masses)
        for q in qs:
            if abs(q) < 1e-9:
                # Z(0) = N(s), the number of cells with nonzero mass
                Z[q][ki] = float(len(masses))
            elif abs(q - 1) < 1e-9:
                # q=1 limit: D_1 = entropy / log(1/s);
                # Z_1 ~ exp(-Sigma p log p) so D_1 = log Z_1 / log(1/s) = -log Z_1 / log s
                Z[q][ki] = float(np.exp(-np.sum(masses *
                                                np.log(masses + 1e-30))))
            else:
                Z[q][ki] = float(np.sum(masses ** q))
    return Z, n_cyls


def extract_tau_q(Z_q, ks, q):
    """tau(q) = -slope of log Z_q(k) vs k * log 2 (since 2^(-k) is the scale)."""
    log_Z = np.log(Z_q)
    log_scale = np.array(ks, dtype=np.float64) * math.log(2.0)
    valid = np.isfinite(log_Z)
    if valid.sum() < 3:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(log_scale[valid], log_Z[valid], 1)
    pred = slope * log_scale[valid] + intercept
    ss_res = np.sum((log_Z[valid] - pred) ** 2)
    ss_tot = np.sum((log_Z[valid] - log_Z[valid].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    tau = -slope
    return float(tau), float(r2), float(slope)


def legendre_transform_tau(qs, taus):
    """alpha = dtau/dq, f(alpha) = q*alpha - tau."""
    qs = np.asarray(qs)
    taus = np.asarray(taus)
    alphas = np.gradient(taus, qs)
    f_alphas = qs * alphas - taus
    return alphas, f_alphas


# ----------------------------------------------------------------
# Weighted Cantor model fits
# ----------------------------------------------------------------
def one_scale_cantor_Dq(q, p, l=0.5):
    if abs(q - 1) < 1e-9:
        h = -(p*np.log(p+1e-12) + (1-p)*np.log(1-p+1e-12))
        return h / np.log(1/l)
    return np.log(p**q + (1-p)**q) / ((1 - q) * np.log(l))


def two_scale_cantor_Dq(q, p1, p2, l1, l2):
    def f(D):
        return p1**q * l1**((q-1)*D) + p2**q * l2**((q-1)*D) - 1
    try:
        return brentq(f, -5, 5)
    except ValueError:
        return np.nan


def fit_one_scale(qs, Dq_emp, l=0.5):
    qs = np.array(qs); Dq_emp = np.array(Dq_emp)
    def loss(params):
        p = params[0]
        if not (0.01 < p < 0.99):
            return 1e9
        Dq_mod = np.array([one_scale_cantor_Dq(q, p, l) for q in qs])
        valid = np.isfinite(Dq_mod) & np.isfinite(Dq_emp)
        return float(np.sum((Dq_mod[valid] - Dq_emp[valid]) ** 2))
    best = None
    for p0 in [0.2, 0.4, 0.5, 0.6, 0.8]:
        res = minimize(loss, [p0], method='Nelder-Mead',
                       options={'xatol': 1e-7, 'fatol': 1e-9})
        if best is None or res.fun < best.fun:
            best = res
    return best.x[0], best.fun


def fit_two_scale(qs, Dq_emp):
    qs = np.array(qs); Dq_emp = np.array(Dq_emp)
    def loss(params):
        p1, l1, l2 = params
        p2 = 1 - p1
        if not (0.01 < p1 < 0.99 and 0.05 < l1 < 0.95 and 0.05 < l2 < 0.95):
            return 1e9
        Dq_mod = np.array([two_scale_cantor_Dq(q, p1, p2, l1, l2) for q in qs])
        valid = np.isfinite(Dq_mod) & np.isfinite(Dq_emp)
        if valid.sum() < 4:
            return 1e9
        return float(np.sum((Dq_mod[valid] - Dq_emp[valid]) ** 2))
    best = None
    for x0 in [(0.5, 0.4, 0.4), (0.7, 0.3, 0.5), (0.3, 0.5, 0.3),
              (0.6, 0.45, 0.35), (0.5, 0.5, 0.5)]:
        res = minimize(loss, x0, method='Nelder-Mead',
                       options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 4000})
        if best is None or res.fun < best.fun:
            best = res
    return tuple(best.x), best.fun


def main():
    log("=" * 78)
    log("Spatial multifractal analysis of trajectory measure on Z_2")
    log("Z_q(2^k) = sum_C mu(C)^q ~ 2^(-k*tau(q))")
    log("=" * 78)

    # ---------------- Build inverse tree ----------------
    max_value = 1 << 20  # ~1M, fits in memory; tree gives ~380K odd nodes
    log(f"\nStep 1: building inverse tree at max_value = 2^20 = {max_value}")
    t0 = time.time()
    tree = build_inverse_tree(max_value)
    log(f"  tree size: {len(tree):,} nodes, {time.time()-t0:.1f}s")
    sizes = compute_subtree_sizes(tree)
    n_odd = sum(1 for m in tree if m & 1)
    log(f"  odd nodes: {n_odd:,}")

    # Variant (a) weighting: subtree-size (Result 58 best)
    weights = {m: float(sizes[m]) for m in tree if m & 1}
    log(f"  using variant (a) subtree-size weights")

    # ---------------- Partition function ----------------
    ks = list(range(6, 17))  # 2^6 = 64 ... 2^16 = 65,536
    qs = [-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5]
    log(f"\nStep 2: partition function Z_q(2^k), ks={ks}, qs={qs}")
    t0 = time.time()
    Z, n_cyls = partition_function_Zq(weights, ks, qs)
    log(f"  Z_q computation: {time.time()-t0:.1f}s")
    log(f"  n_active_cyls per k: " +
        ", ".join(f"k={k}:{n}" for k, n in zip(ks, n_cyls)))

    # ---------------- Extract tau(q) ----------------
    log("\n" + "=" * 78)
    log("Step 3: tau(q) from log Z_q(k) vs k log 2")
    log("=" * 78)
    tau_q = {}
    r2_q = {}
    for q in qs:
        tau, r2, slope = extract_tau_q(Z[q], ks, q)
        tau_q[q] = tau
        r2_q[q] = r2

    log(f"\n  {'q':>5}  {'tau(q)':>8}  {'D_q':>8}  {'R^2':>6}")
    Dq_emp = []
    for q in qs:
        t = tau_q[q]
        if abs(q) < 1e-9:
            d = -t  # D_0 = -tau(0)
        elif abs(q - 1) < 1e-9:
            # D_1 = log Z_1 / log(1/s); since Z_1 ~ s^tau(1), D_1 = -tau(1)
            d = -t
        else:
            d = t / (q - 1)
        Dq_emp.append(d)
        log(f"  {q:>5.1f}  {t:>8.4f}  {d:>8.4f}  {r2_q[q]:>6.4f}")
    Dq_emp = np.array(Dq_emp)

    # ---------------- Legendre transform ----------------
    log("\n" + "=" * 78)
    log("Step 4: f(alpha) via Legendre transform")
    log("=" * 78)
    tau_arr = np.array([tau_q[q] for q in qs])
    alphas, f_alphas = legendre_transform_tau(qs, tau_arr)

    log(f"\n  {'q':>5}  {'tau(q)':>8}  {'alpha':>8}  {'f(alpha)':>9}")
    for i, q in enumerate(qs):
        log(f"  {q:>5.1f}  {tau_arr[i]:>8.4f}  {alphas[i]:>8.4f}  "
            f"{f_alphas[i]:>9.4f}")

    width = float(np.nanmax(alphas) - np.nanmin(alphas))
    f_max_idx = int(np.nanargmax(f_alphas))
    log(f"\n  alpha range: [{np.nanmin(alphas):.4f}, {np.nanmax(alphas):.4f}], "
        f"width = {width:.4f}")
    log(f"  Peak f(alpha) = {f_alphas[f_max_idx]:.4f} at "
        f"alpha = {alphas[f_max_idx]:.4f} (q = {qs[f_max_idx]:.2f})")
    log(f"  D_0 (support dim) = {f_alphas[f_max_idx]:.4f}, "
        f"alpha_0 = {alphas[f_max_idx]:.4f}")

    # ---------------- Locate Chang 0.68, R23 0.6755 ---------
    log("\n" + "=" * 78)
    log("Step 5: locate Chang/R23/R59 dim values on f(alpha)")
    log("=" * 78)
    targets = {
        "Chang H-dim 0.6800": 0.68,
        "R23 lambda_max 0.6755": 0.6755,
        "R59 dim_q2(k=12) 0.67": 0.67,
        "R59 dim_q2(k=15) 0.54": 0.54,
        "R59 dim_q2(k=7) 0.83": 0.83,
        "support dim alpha_0": float(alphas[f_max_idx]),
    }
    log(f"\n  Looking for q where alpha matches target dim:")
    log(f"  {'target':>27}  {'q*':>6}  {'alpha*':>7}  {'f(alpha*)':>9}")
    valid = np.isfinite(alphas) & np.isfinite(f_alphas)
    qs_arr = np.array(qs)
    for name, target in targets.items():
        if valid.sum() == 0:
            log(f"  {name:>27}  no valid spectrum")
            continue
        diffs = np.abs(alphas[valid] - target)
        idx = int(np.argmin(diffs))
        q_star = qs_arr[valid][idx]
        a_star = alphas[valid][idx]
        f_star = f_alphas[valid][idx]
        log(f"  {name:>27}  {q_star:>6.2f}  {a_star:>7.4f}  {f_star:>9.4f}")

    # Also: locate by D_q (Renyi-q dim) match
    log(f"\n  Looking for q where D_q matches target dim (Renyi):")
    log(f"  {'target':>27}  {'q*':>6}  {'D_q*':>7}")
    Dq_arr = np.array(Dq_emp)
    valid_d = np.isfinite(Dq_arr)
    for name, target in targets.items():
        if valid_d.sum() == 0:
            continue
        diffs = np.abs(Dq_arr[valid_d] - target)
        idx = int(np.argmin(diffs))
        q_star = qs_arr[valid_d][idx]
        d_star = Dq_arr[valid_d][idx]
        log(f"  {name:>27}  {q_star:>6.2f}  {d_star:>7.4f}")

    # ---------------- Cantor fits ----------------
    log("\n" + "=" * 78)
    log("Step 6: weighted Cantor set fits to D_q spectrum")
    log("=" * 78)
    qs_for_fit = [q for q in qs if abs(q - 1) > 1e-6]
    Dq_for_fit = np.array([Dq_emp[i] for i, q in enumerate(qs)
                           if abs(q - 1) > 1e-6])
    p_best, loss1 = fit_one_scale(qs_for_fit, Dq_for_fit, l=0.5)
    log(f"\n  One-scale (l=1/2): p = {p_best:.4f}, RSS = {loss1:.6f}")
    Dq_one = np.array([one_scale_cantor_Dq(q, p_best, 0.5) for q in qs])
    log(f"\n  {'q':>5}  {'D_q^emp':>9}  {'D_q^1sc':>9}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"  {q:>5.1f}  {Dq_emp[i]:>9.4f}  {Dq_one[i]:>9.4f}  "
            f"{Dq_one[i] - Dq_emp[i]:>+7.4f}")

    (p1, l1, l2), loss2 = fit_two_scale(qs_for_fit, Dq_for_fit)
    log(f"\n  Two-scale: p1={p1:.4f}, p2={1-p1:.4f}, l1={l1:.4f}, l2={l2:.4f}")
    log(f"    RSS = {loss2:.6f}")
    Dq_two = np.array([two_scale_cantor_Dq(q, p1, 1-p1, l1, l2) for q in qs])
    log(f"\n  {'q':>5}  {'D_q^emp':>9}  {'D_q^2sc':>9}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"  {q:>5.1f}  {Dq_emp[i]:>9.4f}  {Dq_two[i]:>9.4f}  "
            f"{Dq_two[i] - Dq_emp[i]:>+7.4f}")

    log(f"\n  RSS: 1-scale = {loss1:.4f}, 2-scale = {loss2:.4f}")
    if loss1 > 0:
        log(f"  Reduction by 2-scale: {(loss1 - loss2)/loss1*100:.1f}%")

    # ---------------- Save outputs ----------------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)
    with open(OUT / "mfdfa_spatial_tau_q.csv", "w") as f:
        f.write("q,tau_q,D_q,R2\n")
        for i, q in enumerate(qs):
            f.write(f"{q:.2f},{tau_q[q]:.6f},{Dq_emp[i]:.6f},{r2_q[q]:.6f}\n")
    log("  [wrote] mfdfa_spatial_tau_q.csv")

    with open(OUT / "mfdfa_spatial_f_alpha.csv", "w") as f:
        f.write("q,tau,alpha,f_alpha\n")
        for i, q in enumerate(qs):
            f.write(f"{q:.2f},{tau_arr[i]:.6f},{alphas[i]:.6f},{f_alphas[i]:.6f}\n")
    log("  [wrote] mfdfa_spatial_f_alpha.csv")

    with open(OUT / "mfdfa_spatial_Z_q.csv", "w") as f:
        f.write("q," + ",".join(f"k={k}" for k in ks) + "\n")
        for q in qs:
            f.write(f"{q:.2f}," + ",".join(f"{Z[q][i]:.6e}"
                                          for i in range(len(ks))) + "\n")
    log("  [wrote] mfdfa_spatial_Z_q.csv")

    with open(OUT / "mfdfa_spatial_cantor_fits.csv", "w") as f:
        f.write("model,params,RSS\n")
        f.write(f"one_scale,p={p_best:.6f}|l=0.5,{loss1:.6f}\n")
        f.write(f"two_scale,p1={p1:.6f}|p2={1-p1:.6f}|l1={l1:.6f}|l2={l2:.6f},"
                f"{loss2:.6f}\n")
    log("  [wrote] mfdfa_spatial_cantor_fits.csv")

    # ---------------- Verdict ---------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    r2_mean = float(np.nanmean([r2_q[q] for q in qs]))
    log(f"  Mean R^2 of log Z_q vs k log 2: {r2_mean:.4f}")
    log(f"  alpha-spectrum width: {width:.4f}")
    log(f"  D_0 (support dim): {f_alphas[f_max_idx]:.4f}")
    log(f"  D_2 (correlation dim, q=2): {Dq_emp[qs.index(2)]:.4f}")
    log(f"  One-scale RSS: {loss1:.4f}")
    log(f"  Two-scale RSS: {loss2:.4f}")
    log(f"")
    if r2_mean < 0.9:
        log(f"  Outcome: (gamma) — spatial scaling fits poor (R^2 < 0.9)")
    elif loss1 < 0.05 or loss2 < 0.02:
        log(f"  Outcome: (alpha) — clean spectrum + Cantor fit lands.")
    elif width > 0.1:
        log(f"  Outcome: (beta) — clean MF spectrum, no clean Cantor closed form.")
    else:
        log(f"  Outcome: monofractal (alpha-width < 0.1).")

    (OUT / "mfdfa_spatial_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] mfdfa_spatial_log.txt")


if __name__ == "__main__":
    main()
