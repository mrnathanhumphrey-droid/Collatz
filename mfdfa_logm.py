"""
mfdfa_logm.py — Result 61 follow-up.

Re-run MF-DFA with log_2(m_t) as the time series instead of v_t.

Why: integer-valued v_t produces F^2 ~ 0 segments under linear detrending at
small scales, which blows up negative-q moments. log_2(m_t) is naturally
continuous-valued (m goes through many distinct integer values along an orbit)
and behaves Brownian-like with descent drift -log_2(4/3) ~ -0.415 per step.
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
K_H = 3.0 / np.log(4.0/3.0)
MAX_VAL = np.int64(2**62)


@njit(parallel=True, cache=True)
def walk_logm_seqs(starts, max_T):
    n = len(starts)
    log_m_seq = np.full((n, max_T), np.nan, dtype=np.float64)
    seq_count = np.zeros(n, dtype=np.int32)
    sigma_arr = np.full(n, -1, dtype=np.int32)
    log_n_arr = np.zeros(n, dtype=np.float64)
    LOG_2_64 = np.float64(np.log(2.0))
    for i in prange(n):
        m = np.int64(starts[i])
        log_n_arr[i] = np.log(np.float64(m)) / LOG_2_64
        sigma = 0
        T = 0
        failed = False
        while (m & 1) == 0 and m > 1:
            m >>= 1
            sigma += 1
        if m == 1:
            sigma_arr[i] = sigma
            seq_count[i] = 0
            continue
        # Record initial log_2 m
        log_m_seq[i, 0] = np.log(np.float64(m)) / LOG_2_64
        T = 1
        while m != 1 and T < max_T:
            if m > MAX_VAL // 3:
                failed = True
                break
            x = 3 * m + 1
            v = 0
            while (x & 1) == 0:
                x >>= 1
                v += 1
            log_m_seq[i, T] = np.log(np.float64(x)) / LOG_2_64
            sigma += 1 + v
            T += 1
            m = x
        if not failed and m == 1:
            sigma_arr[i] = sigma
            seq_count[i] = T
    return log_m_seq, seq_count, sigma_arr, log_n_arr


@njit(cache=True, parallel=True)
def f2_segments_dfa1(Y, s):
    """
    Compute F^2(nu, s) for non-overlapping segments with DFA1 (linear detrending).
    Forward + reverse pass.
    """
    N = len(Y)
    N_s = N // s
    out = np.empty(2 * N_s, dtype=np.float64)
    x_seg = np.arange(s).astype(np.float64)
    x_mean = (s - 1) / 2.0
    x_var = 0.0
    for j in range(s):
        xc = x_seg[j] - x_mean
        x_var += xc * xc
    for nu in prange(N_s):
        seg = Y[nu*s : (nu+1)*s]
        y_mean = 0.0
        for j in range(s):
            y_mean += seg[j]
        y_mean /= s
        cov = 0.0
        for j in range(s):
            cov += (x_seg[j] - x_mean) * (seg[j] - y_mean)
        b = cov / x_var
        a = y_mean - b * x_mean
        ss = 0.0
        for j in range(s):
            r = seg[j] - (a + b * x_seg[j])
            ss += r * r
        out[nu] = ss / s
    for nu in prange(N_s):
        end = N - nu * s
        start = end - s
        seg = Y[start:end]
        y_mean = 0.0
        for j in range(s):
            y_mean += seg[j]
        y_mean /= s
        cov = 0.0
        for j in range(s):
            cov += (x_seg[j] - x_mean) * (seg[j] - y_mean)
        b = cov / x_var
        a = y_mean - b * x_mean
        ss = 0.0
        for j in range(s):
            r = seg[j] - (a + b * x_seg[j])
            ss += r * r
        out[N_s + nu] = ss / s
    return out


def mfdfa(x, scales, qs):
    Y = np.cumsum(x - np.mean(x))
    F_qs = {q: np.zeros(len(scales), dtype=np.float64) for q in qs}
    for i_s, s in enumerate(scales):
        F2 = f2_segments_dfa1(Y, int(s))
        F2 = F2[F2 > 1e-30]
        if len(F2) == 0:
            for q in qs:
                F_qs[q][i_s] = np.nan
            continue
        for q in qs:
            if abs(q) < 1e-9:
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
        pred = slope * log_s[valid] + intercept
        ss_res = np.sum((log_F[valid] - pred) ** 2)
        ss_tot = np.sum((log_F[valid] - log_F[valid].mean()) ** 2)
        r2_q[q] = float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan
        n_pts = valid.sum()
        sxx = np.sum((log_s[valid] - log_s[valid].mean()) ** 2)
        h_se[q] = float(np.sqrt(ss_res / max(n_pts - 2, 1) / sxx)) if sxx > 0 else np.nan
    return F_qs, h_q, h_se, tau_q, r2_q


def legendre_transform(qs, h_arr):
    qs = np.asarray(qs, dtype=np.float64)
    h_arr = np.asarray(h_arr, dtype=np.float64)
    h_prime = np.gradient(h_arr, qs)
    alphas = h_arr + qs * h_prime
    f_alphas = qs * (alphas - h_arr) + 1.0
    return alphas, f_alphas


def one_scale_cantor_Dq(q, p, l=0.5):
    if abs(q - 1) < 1e-9:
        h = -(p * np.log(p + 1e-12) + (1 - p) * np.log(1 - p + 1e-12))
        return h / np.log(1 / l)
    return np.log(p**q + (1 - p)**q) / ((1 - q) * np.log(l))


def two_scale_cantor_Dq(q, p1, p2, l1, l2):
    def f(D):
        return p1**q * l1**((q-1)*D) + p2**q * l2**((q-1)*D) - 1
    try:
        return brentq(f, -5, 5)
    except ValueError:
        return np.nan


def fit_one_scale(qs, h_emp, l=0.5):
    qs = np.array(qs, dtype=np.float64)
    h_emp = np.array(h_emp, dtype=np.float64)
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9,
                      (qs * h_emp - 1) / np.where(np.abs(qs - 1) > 1e-9, qs - 1, 1.0),
                      h_emp)
    def loss(params):
        p = params[0]
        if p <= 0 or p >= 1:
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


def fit_two_scale(qs, h_emp):
    qs = np.array(qs, dtype=np.float64)
    h_emp = np.array(h_emp, dtype=np.float64)
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9,
                      (qs * h_emp - 1) / np.where(np.abs(qs - 1) > 1e-9, qs - 1, 1.0),
                      h_emp)
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
              (0.6, 0.45, 0.35), (0.5, 0.5, 0.5), (0.65, 0.4, 0.4)]:
        res = minimize(loss, x0, method='Nelder-Mead',
                       options={'xatol': 1e-6, 'fatol': 1e-8, 'maxiter': 4000})
        if best is None or res.fun < best.fun:
            best = res
    return tuple(best.x), best.fun


def locate_dim_on_falpha(target, alphas, f_alphas, qs):
    valid = np.isfinite(alphas) & np.isfinite(f_alphas)
    if valid.sum() == 0:
        return None, None, None
    diffs = np.abs(alphas[valid] - target)
    idx = np.argmin(diffs)
    qs_v = np.array(qs)[valid]
    return float(qs_v[idx]), float(alphas[valid][idx]), float(f_alphas[valid][idx])


def main():
    log("=" * 78)
    log("MF-DFA on log_2(m_t) sequence (continuous-valued, smooth Brownian-like)")
    log("=" * 78)

    # ------------------ Step 1: walk orbits, log_2(m_t) sequences ----
    log2N = 32
    N = 1 << log2N
    n_orbits = 600_000
    max_T = 200
    log(f"\nStep 1: walking {n_orbits:,} orbits at N=2^{log2N}, max_T={max_T}")

    rng = np.random.default_rng(20260503)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1
    t0 = time.time()
    log_m_seq, seq_count, sigma_arr, log_n_arr = walk_logm_seqs(starts, max_T)
    log(f"  walk: {time.time()-t0:.1f}s")

    ok = (sigma_arr > 0) & (seq_count > 0)
    seq_count = seq_count[ok]
    sigma_arr = sigma_arr[ok]
    log_n_arr = log_n_arr[ok]
    log_m_seq = log_m_seq[ok]
    n_ok = len(seq_count)
    log(f"  ok orbits: {n_ok:,}")

    total_v = int(seq_count.sum())
    log(f"  total log m values: {total_v:,}")
    flat_m = np.empty(total_v, dtype=np.float64)
    pos = 0
    for i in range(n_ok):
        c = int(seq_count[i])
        flat_m[pos:pos + c] = log_m_seq[i, :c]
        pos += c

    log(f"  flat_m stats: mean={flat_m.mean():.4f}, std={flat_m.std():.4f}")

    # Use INCREMENTS log_2(m_{t+1}) - log_2(m_t) — these are Brownian-like
    # increments with mean drift -log_2(4/3) and zero-mean fluctuations
    # First take per-orbit increments to avoid jumps between unrelated orbits
    flat_inc = np.empty(total_v - n_ok, dtype=np.float64)
    pos = 0
    for i in range(n_ok):
        c = int(seq_count[i])
        if c < 2:
            continue
        d = log_m_seq[i, 1:c] - log_m_seq[i, :c-1]
        flat_inc[pos:pos + c - 1] = d
        pos += c - 1
    flat_inc = flat_inc[:pos]
    log(f"  flat_inc (per-orbit increments) stats: "
        f"n={len(flat_inc):,}, mean={flat_inc.mean():.4f}, "
        f"std={flat_inc.std():.4f}")
    log(f"  (theoretical drift = -log_2(4/3) = {-np.log2(4/3):.4f})")

    # ------------------ MF-DFA on increment series ----
    scales = np.array([16, 32, 64, 128, 256, 512, 1024, 2048, 4096],
                      dtype=np.int64)
    qs = np.array([-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5],
                  dtype=np.float64)

    log(f"\nStep 2-4: MF-DFA on log m increment series, scales={list(scales)}")
    t0 = time.time()
    F_qs, h_q, h_se, tau_q, r2_q = mfdfa(flat_inc, scales, qs)
    log(f"  MF-DFA: {time.time()-t0:.1f}s")

    log(f"\n  {'q':>5}  {'h(q)':>8}  {'SE':>7}  {'R^2':>6}  {'tau(q)':>8}")
    for q in qs:
        log(f"  {q:>5.1f}  {h_q[q]:>8.4f}  {h_se[q]:>7.4f}  "
            f"{r2_q[q]:>6.4f}  {tau_q[q]:>8.4f}")

    h_arr = np.array([h_q[q] for q in qs])
    delta_h = h_arr.max() - h_arr.min()
    r2_mean = float(np.nanmean([r2_q[q] for q in qs]))
    log(f"\n  Multifractal width Delta h = {delta_h:.4f}")
    log(f"  Mean R^2 = {r2_mean:.4f}")

    # ------------------ Step 5: f(alpha) ------
    alphas, f_alphas = legendre_transform(qs, h_arr)

    log("\n" + "=" * 78)
    log("Step 5: f(alpha) singularity spectrum")
    log("=" * 78)
    log(f"\n  {'q':>5}  {'h(q)':>8}  {'alpha':>8}  {'f(alpha)':>9}")
    for i, q in enumerate(qs):
        log(f"  {q:>5.1f}  {h_arr[i]:>8.4f}  {alphas[i]:>8.4f}  "
            f"{f_alphas[i]:>9.4f}")
    width = alphas.max() - alphas.min()
    f_max_idx = int(np.nanargmax(f_alphas))
    log(f"\n  alpha range: [{alphas.min():.4f}, {alphas.max():.4f}], "
        f"width = {width:.4f}")
    log(f"  Peak f(alpha) = {f_alphas[f_max_idx]:.4f} at "
        f"alpha = {alphas[f_max_idx]:.4f} (q = {qs[f_max_idx]:.2f})")

    # ------------------ Step 6: locate prior dim values --
    log("\n" + "=" * 78)
    log("Step 6: locate prior dim values on f(alpha)")
    log("=" * 78)
    targets = {
        "Chang H-dim 0.6800": 0.68,
        "R23 lambda_max 0.6755": 0.6755,
        "R59 dim_q2(k=12) 0.67": 0.67,
        "R59 dim_q2(k=15) 0.54": 0.54,
        "R59 dim_q2(k=7) 0.83": 0.83,
    }
    log(f"\n  {'target':>25}  {'q*':>6}  {'alpha*':>7}  {'f(alpha*)':>9}")
    for name, target in targets.items():
        q_star, alpha_star, f_star = locate_dim_on_falpha(
            target, alphas, f_alphas, qs)
        if q_star is not None:
            log(f"  {name:>25}  {q_star:>6.2f}  {alpha_star:>7.4f}  "
                f"{f_star:>9.4f}")

    # ------------------ Step 7: weighted Cantor fits -----
    log("\n" + "=" * 78)
    log("Step 7: weighted Cantor set fits")
    log("=" * 78)
    p_best, loss1 = fit_one_scale(qs, h_arr, l=0.5)
    log(f"\n  One-scale (l=1/2): p = {p_best:.4f}, RSS = {loss1:.6f}")
    Dq_emp = np.where(np.abs(qs - 1) > 1e-9,
                      (qs * h_arr - 1) / np.where(np.abs(qs-1)>1e-9, qs-1, 1.0),
                      h_arr)
    Dq_one = np.array([one_scale_cantor_Dq(q, p_best, 0.5) for q in qs])
    log(f"\n    {'q':>5}  {'D_q^emp':>9}  {'D_q^1sc':>9}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"    {q:>5.1f}  {Dq_emp[i]:>9.4f}  {Dq_one[i]:>9.4f}  "
            f"{Dq_one[i] - Dq_emp[i]:>+7.4f}")

    (p1, l1, l2), loss2 = fit_two_scale(qs, h_arr)
    log(f"\n  Two-scale: p1={p1:.4f}, p2={1-p1:.4f}, l1={l1:.4f}, l2={l2:.4f}")
    log(f"    RSS = {loss2:.6f}")
    Dq_two = np.array([two_scale_cantor_Dq(q, p1, 1-p1, l1, l2) for q in qs])
    log(f"\n    {'q':>5}  {'D_q^emp':>9}  {'D_q^2sc':>9}  {'diff':>7}")
    for i, q in enumerate(qs):
        log(f"    {q:>5.1f}  {Dq_emp[i]:>9.4f}  {Dq_two[i]:>9.4f}  "
            f"{Dq_two[i] - Dq_emp[i]:>+7.4f}")

    log(f"\n  RSS: 1-scale = {loss1:.4f}, 2-scale = {loss2:.4f}")
    log(f"  Reduction by 2-scale: {(loss1 - loss2)/loss1*100:.1f}%")

    # ------------------ Step 8: sigma-band conditional ---
    log("\n" + "=" * 78)
    log("Step 8: sigma-band conditional MF-DFA on log m increments")
    log("=" * 78)
    sigma_resid = sigma_arr - K_H * log_n_arr * LOG_2
    edges = np.percentile(sigma_resid, [25, 50, 75, 95])
    band = np.digitize(sigma_resid, edges)
    band_names = ['0-25', '25-50', '50-75', '75-95', '95-100']

    h_band = {}
    f_alpha_band = {}
    width_band = {}
    for b in range(5):
        mask = band == b
        if mask.sum() == 0:
            continue
        sel = np.where(mask)[0]
        n_band = len(sel)
        # Build per-orbit increments for this band
        flat_b = []
        for i in sel:
            c = int(seq_count[i])
            if c < 2:
                continue
            d = log_m_seq[i, 1:c] - log_m_seq[i, :c-1]
            flat_b.append(d)
        if not flat_b:
            continue
        flat_b = np.concatenate(flat_b)
        if len(flat_b) < 8 * scales[-1]:
            log(f"  band {band_names[b]}: insufficient ({len(flat_b):,})")
            continue
        log(f"  band {band_names[b]}: n_orbits={n_band:,}, n_inc={len(flat_b):,}")
        F_qs_b, h_q_b, h_se_b, tau_q_b, r2_q_b = mfdfa(flat_b, scales, qs)
        h_band[b] = h_q_b
        h_arr_b = np.array([h_q_b[q] for q in qs])
        a_b, f_b = legendre_transform(qs, h_arr_b)
        f_alpha_band[b] = (a_b, f_b)
        width_band[b] = float(a_b.max() - a_b.min())

    log(f"\n  Band-conditional h(q) and alpha-width:")
    log(f"  {'band':>8}  " +
        "  ".join(f"q={q:>5.1f}" for q in [-3, -1, 0, 1, 2, 3, 5]) +
        "  alpha-w")
    for b in range(5):
        if b not in h_band:
            continue
        line = f"  {band_names[b]:>8}"
        for q in [-3, -1, 0, 1, 2, 3, 5]:
            line += f"  {h_band[b].get(q, np.nan):>7.4f}"
        line += f"  {width_band[b]:>7.4f}"
        log(line)

    # ------------------ Save outputs --------
    log("\n" + "=" * 78)
    log("Saving outputs")
    log("=" * 78)

    with open(OUT / "mfdfa_logm_h_q.csv", "w") as f:
        f.write("q,h_q,h_se,r2,tau_q\n")
        for q in qs:
            f.write(f"{q:.2f},{h_q[q]:.6f},{h_se[q]:.6f},{r2_q[q]:.6f},"
                    f"{tau_q[q]:.6f}\n")
    log("  [wrote] mfdfa_logm_h_q.csv")

    with open(OUT / "mfdfa_logm_f_alpha.csv", "w") as f:
        f.write("q,h_q,alpha,f_alpha\n")
        for i, q in enumerate(qs):
            f.write(f"{q:.2f},{h_arr[i]:.6f},{alphas[i]:.6f},{f_alphas[i]:.6f}\n")
    log("  [wrote] mfdfa_logm_f_alpha.csv")

    with open(OUT / "mfdfa_logm_cantor_fits.csv", "w") as f:
        f.write("model,params,RSS\n")
        f.write(f"one_scale,p={p_best:.6f}|l=0.5,{loss1:.6f}\n")
        f.write(f"two_scale,p1={p1:.6f}|p2={1-p1:.6f}|l1={l1:.6f}|l2={l2:.6f},"
                f"{loss2:.6f}\n")
    log("  [wrote] mfdfa_logm_cantor_fits.csv")

    with open(OUT / "mfdfa_logm_band_conditional.csv", "w") as f:
        f.write("band,q,h_q\n")
        for b in range(5):
            if b not in h_band:
                continue
            for q in qs:
                f.write(f"{band_names[b]},{q:.2f},"
                        f"{h_band[b].get(q, np.nan):.6f}\n")
    log("  [wrote] mfdfa_logm_band_conditional.csv")

    # ------------------ Verdict --------
    log("\n" + "=" * 78)
    log("VERDICT")
    log("=" * 78)
    log(f"  Mean R^2: {r2_mean:.4f}")
    log(f"  Multifractal width Delta h: {delta_h:.4f}")
    log(f"  alpha-spectrum width: {width:.4f}")
    log(f"  One-scale RSS: {loss1:.4f}")
    log(f"  Two-scale RSS: {loss2:.4f}")
    log(f"")
    if r2_mean < 0.95:
        log(f"  Outcome: (gamma) — MF-DFA fits poor (R^2 < 0.95)")
    elif loss1 < 0.05 or loss2 < 0.02:
        log(f"  Outcome: (alpha) — clean spectrum + Cantor fit lands.")
    elif delta_h > 0.1:
        log(f"  Outcome: (beta) — clean MF spectrum, no clean Cantor fit.")
    else:
        log(f"  Outcome: monofractal (Delta h < 0.1) — series essentially "
            f"single-Hurst.")

    (OUT / "mfdfa_logm_log.txt").write_text(
        "\n".join(log_lines), encoding="utf-8"
    )
    log("\n  [wrote] mfdfa_logm_log.txt")


if __name__ == "__main__":
    main()
