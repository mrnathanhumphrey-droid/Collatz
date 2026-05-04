"""
mfdfa_trajectory_measure.py — Result 61 candidate.

Apply Kantelhardt's MF-DFA (multifractal detrended fluctuation analysis)
to the per-orbit Syracuse v_2 sequence and extract h(q), tau(q), f(alpha).

Following Kantelhardt 2002/2008 (arXiv:0804.0747), Section 6.3:
  Y(j)    = Sigma_{i<=j} (x_i - <x>)                        cumulative profile
  F^2(nu,s) = (1/s) Sigma (Y(i) - p_nu(i))^2                segment-wise detrended var (DFA1)
  F_q(s)  = ((1/N_s) Sigma_nu [F^2(nu,s)]^(q/2))^(1/q)       q-th-order fluctuation
  F_0(s)  = exp((1/(2 N_s)) Sigma_nu log F^2(nu,s))          q=0 limit
  log F_q(s) ~ h(q) log s                                    generalized Hurst
  tau(q)  = q*h(q) - 1                                       Renyi exponent
  alpha   = h(q) + q*h'(q)                                   singularity strength
  f(alpha)= q*(alpha - h(q)) + 1                             Hausdorff of mu-set

Time series: concatenate v_2 values (Syracuse step = (3m+1)/2^v) across many
orbits at N=2^32. Single 30M-point sequence, MF-DFA at scales s in {16..2048},
q in {-5..5}.

Then:
- Locate Chang H-dim 0.68 and Result 23 lambda_max-derived 0.6755 on f(alpha)
- Fit one-scale weighted Cantor set D_q(p) and Macek-Wojcik two-scale model
- sigma-band conditional MF-DFA across 5 bands (Result 22 quantile bands)
"""
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit, prange
from scipy.optimize import brentq, minimize

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(r"C:\Collatz")
log_lines = []
def log(s):
    print(s, flush=True)
    log_lines.append(s)

LOG_2 = np.log(2.0)
LOG_3 = np.log(3.0)
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


# ---------------------------------------------------------------------
# Walker: concatenate v_2 values across many orbits
# ---------------------------------------------------------------------
@njit(parallel=True, cache=True)
def walk_v_seqs(starts, max_T, max_v_per_orbit):
    """
    For each orbit, record:
      v_seq[i, 0..t_i-1] = v_2 at each Syracuse step
      v_count[i] = t_i (number of Syracuse steps)
      sigma_arr[i] = total Collatz steps (or -1 if orbit failed)
    Plus per-orbit log_n for sigma-band computation.
    """
    n = len(starts)
    v_seq = np.full((n, max_v_per_orbit), -1, dtype=np.int8)
    v_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    log_n_arr = np.zeros(n, dtype=np.float64)
    for i in prange(n):
        m = np.int64(starts[i])
        log_n_arr[i] = np.log(np.float64(m))
        sigma = 0
        T = 0
        failed = False
        while (m & 1) == 0 and m > 1:
            m >>= 1
            sigma += 1
        if m == 1:
            sigma_arr[i] = sigma
            v_count[i] = 0
            continue
        while m != 1 and T < max_T and T < max_v_per_orbit:
            if m > MAX_VAL // 3:
                failed = True
                break
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1
                v += 1
            v_seq[i, T] = v if v < 127 else 127
            sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            v_count[i] = T
    return v_seq, v_count, sigma_arr, log_n_arr


# ---------------------------------------------------------------------
# MF-DFA core
# ---------------------------------------------------------------------
@njit(cache=True, parallel=True)
def f2_segments(Y, s, order):
    """
    Compute F^2(nu, s) for all segments nu, with polynomial detrending of given order.
    Return array of length 2*N_s (forward + reverse).
    """
    N = len(Y)
    N_s = N // s
    out = np.empty(2 * N_s, dtype=np.float64)
    x_seg = np.arange(s).astype(np.float64)
    # Forward
    for nu in prange(N_s):
        seg = Y[nu*s : (nu+1)*s]
        # Polyfit degree=order via numpy is not numba-safe;
        # do explicit linear (order=1) fit by closed form
        # x_seg sum, mean, etc.
        x_mean = (s - 1) / 2.0
        y_mean = 0.0
        for j in range(s):
            y_mean += seg[j]
        y_mean /= s
        num = 0.0
        den = 0.0
        for j in range(s):
            xc = x_seg[j] - x_mean
            num += xc * (seg[j] - y_mean)
            den += xc * xc
        b = num / den
        a = y_mean - b * x_mean
        # Detrended sum of squares
        ss = 0.0
        for j in range(s):
            r = seg[j] - (a + b * x_seg[j])
            ss += r * r
        out[nu] = ss / s
    # Reverse
    Y_rev_offset = N - N_s * s  # may have leftover at start
    for nu in prange(N_s):
        # nu-th segment from the end
        end = N - nu * s
        start = end - s
        seg = Y[start:end]
        x_mean = (s - 1) / 2.0
        y_mean = 0.0
        for j in range(s):
            y_mean += seg[j]
        y_mean /= s
        num = 0.0
        den = 0.0
        for j in range(s):
            xc = x_seg[j] - x_mean
            num += xc * (seg[j] - y_mean)
            den += xc * xc
        b = num / den
        a = y_mean - b * x_mean
        ss = 0.0
        for j in range(s):
            r = seg[j] - (a + b * x_seg[j])
            ss += r * r
        out[N_s + nu] = ss / s
    return out


def mfdfa(x, scales, qs):
    """
    Standard MF-DFA1 (linear detrending).
    Returns:
      F_qs[q] = array of F_q(s) over scales
      h_q, h_se, tau_q : dicts indexed by q
      r2_q : R^2 of log-log fit per q
    """
    Y = np.cumsum(x - np.mean(x))
    N = len(Y)
    F_qs = {q: np.zeros(len(scales), dtype=np.float64) for q in qs}
    for i_s, s in enumerate(scales):
        F2 = f2_segments(Y, s, 1)
        F2 = F2[F2 > 0]  # filter degenerate segments
        if len(F2) == 0:
            for q in qs:
                F_qs[q][i_s] = np.nan
            continue
        for q in qs:
            if abs(q) < 1e-9:  # q=0 limit
                F_qs[q][i_s] = float(np.exp(0.5 * np.mean(np.log(F2))))
            else:
                F_qs[q][i_s] = float(np.mean(F2 ** (q / 2.0)) ** (1.0 / q))
    h_q = {}
    h_se = {}
    tau_q = {}
    r2_q = {}
    log_s = np.log(scales)
    for q in qs:
        log_F = np.log(F_qs[q])
        valid = np.isfinite(log_F)
        if valid.sum() < 3:
            h_q[q] = np.nan
            h_se[q] = np.nan
            tau_q[q] = np.nan
            r2_q[q] = np.nan
            continue
        slope, intercept = np.polyfit(log_s[valid], log_F[valid], 1)
        h_q[q] = float(slope)
        tau_q[q] = q * slope - 1
        # R^2 and SE
        pred = slope * log_s[valid] + intercept
        ss_res = np.sum((log_F[valid] - pred) ** 2)
        ss_tot = np.sum((log_F[valid] - log_F[valid].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        n_pts = valid.sum()
        sxx = np.sum((log_s[valid] - log_s[valid].mean()) ** 2)
        se = np.sqrt(ss_res / max(n_pts - 2, 1) / sxx) if sxx > 0 else np.nan
        h_se[q] = float(se)
        r2_q[q] = float(r2)
    return F_qs, h_q, h_se, tau_q, r2_q


def legendre_transform(qs, h_arr):
    """
    alpha = h + q*h', f = q*(alpha - h) + 1
    """
    qs = np.asarray(qs, dtype=np.float64)
    h_arr = np.asarray(h_arr, dtype=np.float64)
    h_prime = np.gradient(h_arr, qs)
    alphas = h_arr + qs * h_prime
    f_alphas = qs * (alphas - h_arr) + 1.0
    return alphas, f_alphas


# ---------------------------------------------------------------------
# Weighted Cantor model fits
# ---------------------------------------------------------------------
def one_scale_cantor_Dq(q, p, l=0.5):
    """
    D_q for one-scale weighted Cantor set: two children with weights p, 1-p
    each at scale l (default 1/2 — symmetric Cantor support).
    D_q = log(p^q + (1-p)^q) / ((1-q) * log(l))
    For q=1: D_1 = -(p log p + (1-p) log(1-p)) / log(1/l)
    """
    if abs(q - 1) < 1e-9:
        h = -(p * np.log(p + 1e-12) + (1 - p) * np.log(1 - p + 1e-12))
        return h / np.log(1 / l)
    return np.log(p**q + (1 - p)**q) / ((1 - q) * np.log(l))


def two_scale_cantor_Dq(q, p1, p2, l1, l2):
    """
    Two-scale weighted Cantor: p1^q * l1^((q-1)*D) + p2^q * l2^((q-1)*D) = 1
    Solve for D = D_q.
    """
    def f(D):
        return p1**q * l1**((q-1)*D) + p2**q * l2**((q-1)*D) - 1
    try:
        return brentq(f, -5, 5)
    except ValueError:
        return np.nan


def fit_one_scale(qs, h_emp, l=0.5):
    """
    Fit p in (0,1) by minimizing sum(D_q^model - D_q^emp)^2 where D_q = (q*h-1)/(q-1).
    """
    qs = np.array(qs, dtype=np.float64)
    h_emp = np.array(h_emp, dtype=np.float64)
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9, (qs * h_emp - 1) / (qs - 1), h_emp)
    def loss(params):
        p = params[0]
        if p <= 0 or p >= 1:
            return 1e9
        Dq_mod = np.array([one_scale_cantor_Dq(q, p, l) for q in qs])
        valid = np.isfinite(Dq_mod) & np.isfinite(Dq_emp)
        return float(np.sum((Dq_mod[valid] - Dq_emp[valid]) ** 2))
    best = None
    for p0 in [0.3, 0.5, 0.7]:
        res = minimize(loss, [p0], method='Nelder-Mead',
                       options={'xatol': 1e-7, 'fatol': 1e-9})
        if best is None or res.fun < best.fun:
            best = res
    return best.x[0], best.fun


def fit_two_scale(qs, h_emp):
    """
    Fit (p1, p2, l1, l2) with p1+p2=1 (probability) and l1, l2 in (0, 1/2 reasonable).
    """
    qs = np.array(qs, dtype=np.float64)
    h_emp = np.array(h_emp, dtype=np.float64)
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9, (qs * h_emp - 1) / (qs - 1), h_emp)
    def loss(params):
        p1, l1, l2 = params
        p2 = 1 - p1
        if not (0 < p1 < 1 and 0 < l1 < 1 and 0 < l2 < 1):
            return 1e9
        Dq_mod = np.array([two_scale_cantor_Dq(q, p1, p2, l1, l2) for q in qs])
        valid = np.isfinite(Dq_mod) & np.isfinite(Dq_emp)
        return float(np.sum((Dq_mod[valid] - Dq_emp[valid]) ** 2))
    best = None
    for x0 in [(0.5, 0.4, 0.4), (0.7, 0.3, 0.5), (0.3, 0.5, 0.3),
              (0.6, 0.45, 0.35)]:
        res = minimize(loss, x0, method='Nelder-Mead',
                       options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 4000})
        if best is None or res.fun < best.fun:
            best = res
    return tuple(best.x), best.fun


# ---------------------------------------------------------------------
# Locate prior dim values on f(alpha)
# ---------------------------------------------------------------------
def locate_dim_on_falpha(target_dim, alphas, f_alphas, qs):
    """
    For a target dim value, find the q on the spectrum where alpha closest
    to target. Returns q*, alpha*, f(alpha*).
    """
    valid = np.isfinite(alphas) & np.isfinite(f_alphas)
    if valid.sum() == 0:
        return None, None, None
    diffs = np.abs(alphas[valid] - target_dim)
    idx = np.argmin(diffs)
    qs_v = np.array(qs)[valid]
    return float(qs_v[idx]), float(alphas[valid][idx]), float(f_alphas[valid][idx])


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    log("=" * 78)
    log("MF-DFA on per-orbit Syracuse v_2 sequence (Kantelhardt 2002/2008)")
    log("=" * 78)

    # ------------------ Step 1: walk orbits, build v sequence ----------
    log2N = 32
    N = 1 << log2N
    n_orbits = 600_000  # ~30M v values total
    max_T = 200
    log(f"\nStep 1: walking {n_orbits:,} orbits at N=2^{log2N}, max_T={max_T}")

    rng = np.random.default_rng(20260503)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1

    t0 = time.time()
    v_seq, v_count, sigma_arr, log_n_arr = walk_v_seqs(starts, max_T, max_T)
    log(f"  walk: {time.time()-t0:.1f}s")

    ok = (sigma_arr > 0) & (v_count > 0)
    v_count = v_count[ok]
    sigma_arr = sigma_arr[ok]
    log_n_arr = log_n_arr[ok]
    v_seq = v_seq[ok]
    n_ok = len(v_count)
    log(f"  ok orbits: {n_ok:,}")

    # Concatenate
    total_v = int(v_count.sum())
    log(f"  total v values: {total_v:,}")
    flat_v = np.empty(total_v, dtype=np.float64)
    pos = 0
    orbit_offsets = np.zeros(n_ok + 1, dtype=np.int64)
    for i in range(n_ok):
        c = int(v_count[i])
        flat_v[pos:pos + c] = v_seq[i, :c]
        pos += c
        orbit_offsets[i + 1] = pos

    log(f"  flat_v stats: mean={flat_v.mean():.4f}, std={flat_v.std():.4f}")
    log(f"  (theoretical Geom(1/2): mean=2.000, std=1.414)")

    # ------------------ Step 2-4: MF-DFA on full sequence ----------
    scales = np.array([16, 32, 64, 128, 256, 512, 1024, 2048, 4096], dtype=np.int64)
    qs = np.array([-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5],
                  dtype=np.float64)
    log(f"\nStep 2-4: MF-DFA, scales={list(scales)}, qs={list(qs)}")
    t0 = time.time()
    F_qs, h_q, h_se, tau_q, r2_q = mfdfa(flat_v, scales, qs)
    log(f"  MF-DFA: {time.time()-t0:.1f}s")

    log(f"\n  {'q':>5}  {'h(q)':>8}  {'SE':>7}  {'R^2':>6}  {'tau(q)':>8}")
    for q in qs:
        log(f"  {q:>5.1f}  {h_q[q]:>8.4f}  {h_se[q]:>7.4f}  "
            f"{r2_q[q]:>6.4f}  {tau_q[q]:>8.4f}")

    h_arr = np.array([h_q[q] for q in qs])
    tau_arr = np.array([tau_q[q] for q in qs])
    delta_h = h_arr.max() - h_arr.min()
    log(f"\n  Multifractal width Delta h = h_min - h_max = {delta_h:.4f}")
    log(f"  (Delta h > 0.1 suggests multifractality; ~0 = monofractal)")

    # ------------------ Step 5: f(alpha) Legendre transform ----------
    log("\n" + "=" * 78)
    log("Step 5: f(alpha) via Legendre transform")
    log("=" * 78)
    alphas, f_alphas = legendre_transform(qs, h_arr)

    log(f"\n  {'q':>5}  {'h(q)':>8}  {'alpha':>8}  {'f(alpha)':>9}")
    for i, q in enumerate(qs):
        log(f"  {q:>5.1f}  {h_arr[i]:>8.4f}  {alphas[i]:>8.4f}  "
            f"{f_alphas[i]:>9.4f}")

    alpha_max = alphas.max()
    alpha_min = alphas.min()
    width = alpha_max - alpha_min
    f_max_idx = int(np.nanargmax(f_alphas))
    log(f"\n  alpha range: [{alpha_min:.4f}, {alpha_max:.4f}], width = {width:.4f}")
    log(f"  Peak f(alpha) = {f_alphas[f_max_idx]:.4f} at alpha = {alphas[f_max_idx]:.4f} "
        f"(q = {qs[f_max_idx]:.2f})")

    # ------------------ Step 6: locate Chang 0.68, R23 0.6755 ---------
    log("\n" + "=" * 78)
    log("Step 6: locate prior dim values on f(alpha) spectrum")
    log("=" * 78)
    targets = {
        "Chang H-dim 0.6800": 0.68,
        "R23 lambda_max 0.6755": 0.6755,
        "R59 dim_q2(k=12) 0.67": 0.67,
        "R59 dim_q2(k=15) 0.54": 0.54,
        "R59 dim_q2(k=7)  0.83": 0.83,
    }
    log(f"\n  {'target':>25}  {'q*':>6}  {'alpha*':>7}  {'f(alpha*)':>9}")
    for name, target in targets.items():
        q_star, alpha_star, f_star = locate_dim_on_falpha(
            target, alphas, f_alphas, qs)
        if q_star is not None:
            log(f"  {name:>25}  {q_star:>6.2f}  {alpha_star:>7.4f}  "
                f"{f_star:>9.4f}")
        else:
            log(f"  {name:>25}  outside spectrum range")

    # ------------------ Step 7: weighted Cantor model fits -----------
    log("\n" + "=" * 78)
    log("Step 7: weighted Cantor set fits (Macek-Wojcik 2026 candidates)")
    log("=" * 78)

    log("\n  One-scale weighted Cantor (l = 1/2, fit p):")
    p_best, loss1 = fit_one_scale(qs, h_arr, l=0.5)
    log(f"    p = {p_best:.4f}, RSS = {loss1:.6f}")
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9, (qs * h_arr - 1) / (qs - 1), h_arr)
    Dq_one = np.array([one_scale_cantor_Dq(q, p_best, 0.5) for q in qs])
    log(f"\n    {'q':>5}  {'D_q^emp':>8}  {'D_q^1scale':>10}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"    {q:>5.1f}  {Dq_emp[i]:>8.4f}  {Dq_one[i]:>10.4f}  "
            f"{Dq_one[i] - Dq_emp[i]:>+7.4f}")

    log("\n  Two-scale weighted Cantor (fit p1, l1, l2):")
    (p1, l1, l2), loss2 = fit_two_scale(qs, h_arr)
    log(f"    p1 = {p1:.4f}, p2 = {1-p1:.4f}, l1 = {l1:.4f}, l2 = {l2:.4f}")
    log(f"    RSS = {loss2:.6f}")
    Dq_two = np.array([two_scale_cantor_Dq(q, p1, 1-p1, l1, l2) for q in qs])
    log(f"\n    {'q':>5}  {'D_q^emp':>8}  {'D_q^2scale':>10}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"    {q:>5.1f}  {Dq_emp[i]:>8.4f}  {Dq_two[i]:>10.4f}  "
            f"{Dq_two[i] - Dq_emp[i]:>+7.4f}")

    log(f"\n  RSS comparison: 1-scale = {loss1:.6f}, 2-scale = {loss2:.6f}")
    log(f"  RSS reduction by 2-scale: {(loss1 - loss2) / loss1 * 100:.1f}%")

    # ------------------ Step 8: sigma-band conditional MF-DFA ---------
    log("\n" + "=" * 78)
    log("Step 8: sigma-band conditional MF-DFA")
    log("=" * 78)
    sigma_resid = sigma_arr - K_H * log_n_arr
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    h_band = {}
    for b in range(5):
        mask = band == b
        if mask.sum() == 0:
            continue
        # Concatenate v sequences for orbits in this band
        sel = np.where(mask)[0]
        n_band = len(sel)
        flat_b_lens = v_count[sel]
        total_b = int(flat_b_lens.sum())
        flat_b = np.empty(total_b, dtype=np.float64)
        p = 0
        for i in sel:
            c = int(v_count[i])
            flat_b[p:p + c] = v_seq[i, :c]
            p += c
        if total_b < 8 * scales[-1]:
            log(f"  band {band_names[b]}: insufficient data ({total_b:,} < "
                f"{8 * scales[-1]:,}); skipping")
            continue
        scales_b = scales[scales <= total_b // 8]
        log(f"  band {band_names[b]}: n_orbits={n_band:,}, "
            f"total_v={total_b:,}, scales={list(scales_b)}")
        F_qs_b, h_q_b, h_se_b, tau_q_b, r2_q_b = mfdfa(flat_b, scales_b, qs)
        h_band[b] = h_q_b

    # Compare h(q=2) across bands as quick proxy
    log(f"\n  Band-conditional h(q):")
    log(f"  {'band':>8}  " + "  ".join(f"q={q:>5.1f}" for q in
                                       [-3, -1, 0, 1, 2, 3, 5]))
    for b in range(5):
        if b not in h_band:
            continue
        line = f"  {band_names[b]:>8}"
        for q in [-3, -1, 0, 1, 2, 3, 5]:
            line += f"  {h_band[b].get(q, np.nan):>7.4f}"
        log(line)

    # ------------------ Save outputs ---------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)

    with open(OUT / "mfdfa_h_q.csv", "w") as f:
        f.write("q,h_q,h_se,r2,tau_q\n")
        for q in qs:
            f.write(f"{q:.2f},{h_q[q]:.6f},{h_se[q]:.6f},{r2_q[q]:.6f},"
                    f"{tau_q[q]:.6f}\n")
    log("  [wrote] mfdfa_h_q.csv")

    with open(OUT / "mfdfa_tau_q.csv", "w") as f:
        f.write("q,tau_q\n")
        for q in qs:
            f.write(f"{q:.2f},{tau_q[q]:.6f}\n")
    log("  [wrote] mfdfa_tau_q.csv")

    with open(OUT / "mfdfa_f_alpha.csv", "w") as f:
        f.write("q,h_q,alpha,f_alpha\n")
        for i, q in enumerate(qs):
            f.write(f"{q:.2f},{h_arr[i]:.6f},{alphas[i]:.6f},{f_alphas[i]:.6f}\n")
    log("  [wrote] mfdfa_f_alpha.csv")

    with open(OUT / "mfdfa_cantor_fits.csv", "w") as f:
        f.write("model,params,RSS,RSS_per_q\n")
        f.write(f"one_scale,p={p_best:.6f}|l=0.5,{loss1:.6f},{loss1/len(qs):.6f}\n")
        f.write(f"two_scale,p1={p1:.6f}|p2={1-p1:.6f}|l1={l1:.6f}|l2={l2:.6f},"
                f"{loss2:.6f},{loss2/len(qs):.6f}\n")
    log("  [wrote] mfdfa_cantor_fits.csv")

    with open(OUT / "mfdfa_band_conditional.csv", "w") as f:
        f.write("band,q,h_q\n")
        for b in range(5):
            if b not in h_band:
                continue
            for q in qs:
                v = h_band[b].get(q, np.nan)
                f.write(f"{band_names[b]},{q:.2f},{v:.6f}\n")
    log("  [wrote] mfdfa_band_conditional.csv")

    # ------------------ Verdict ------------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    r2_mean = float(np.nanmean([r2_q[q] for q in qs]))
    log(f"  Mean R^2 of log-log fits: {r2_mean:.4f}")
    log(f"  Multifractal width Delta h: {delta_h:.4f}")
    log(f"  alpha-spectrum width: {width:.4f}")
    log(f"  One-scale Cantor fit RSS: {loss1:.6f}")
    log(f"  Two-scale Cantor fit RSS: {loss2:.6f}")
    log(f"")
    if r2_mean < 0.95:
        log(f"  Outcome: (gamma) — MF-DFA fits poor (R^2 < 0.95)")
    elif loss1 < 0.05 or loss2 < 0.02:
        log(f"  Outcome: (alpha) — clean spectrum + Cantor fit lands.")
    else:
        log(f"  Outcome: (beta) — clean spectrum, no clean weighted Cantor fit.")

    (OUT / "mfdfa_log.txt").write_text("\n".join(log_lines), encoding="utf-8")
    log("\n  [wrote] mfdfa_log.txt")


if __name__ == "__main__":
    main()
