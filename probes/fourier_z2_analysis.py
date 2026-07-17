"""
fourier_z2_analysis.py — Fourier analysis of trajectory measure on Z_2.

Computes μ̂(ξ) for ξ = j/2^k dyadic frequencies, fits decay law,
compares to D_q values from R61, tests Bernoulli convolution match,
computes multifractal Fourier σ_q, looks for resonances, compares to Chang's π.

Efficient computation: for each k, μ̂(j/2^k) is the discrete Fourier transform
of M_k[r] = Σ_{m ≡ r mod 2^k} weight(m). FFT in O(2^k log 2^k).
"""
import csv
import math
import os
import sys
import time
from collections import deque, defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


# ============================================================
# Inverse tree builder + subtree sizes
# ============================================================
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


def subtree_sizes(tree):
    children = defaultdict(list)
    for m, info in tree.items():
        if info['parent'] is not None:
            children[info['parent']].append(m)
    by_depth_desc = sorted(tree.keys(), key=lambda m: -tree[m]['depth'])
    size = {m: 1 for m in tree}
    for m in by_depth_desc:
        for c in children[m]:
            size[m] += size[c]
    return size


# ============================================================
# Fourier transform via FFT on cylinder buckets
# ============================================================
def compute_M_k_and_fft(integers, weights, k):
    """For modulus 2^k, compute M_k[r] = Σ weight where integer ≡ r mod 2^k.
    FFT M_k → mu_hat at all dyadic frequencies j/2^k, j=0..2^k-1.
    Returns mu_hat array of length 2^k (complex)."""
    mod = 1 << k
    M = np.zeros(mod, dtype=np.float64)
    rs = integers % mod
    np.add.at(M, rs, weights)
    # FFT: μ̂(j/2^k) = Σ_r M[r] · exp(2πi · j · r / 2^k) / Z
    # numpy fft uses convention exp(-2πi · j · r / N), so we need to take conjugate for our convention.
    # Actually for our μ̂(ξ) = E[exp(2πi ξ x)], with x = m, ξ = j/2^k:
    # μ̂(j/2^k) = (1/Z) Σ_m w(m) exp(2πi · j m / 2^k) = (1/Z) Σ_r M[r] exp(2πi · j r / 2^k)
    # = (1/Z) · IFFT(M)[j] · 2^k (since np.fft.ifft has 1/N factor)
    Z = M.sum()
    mu_hat = np.fft.ifft(M) * mod / Z  # this gives Σ_r M[r] exp(+2πi j r / 2^k) / Z
    return mu_hat


# ============================================================
# Bernoulli convolution analog on Z_2: ν_p = product Ber(p) on bits
# ν̂_p(j/2^k) = ∏_{n=0..k-1} [(1-p) + p · exp(2πi · j · 2^n / 2^k)]
#            = ∏_{n=0..k-1} [(1-p) + p · exp(2πi · j / 2^{k-n})]
# Easier: bit n contributes value 2^n. So ν̂_p(ξ) = ∏_n [(1-p) + p · e^{2πi ξ 2^n}].
# For ξ = j / 2^k, ξ · 2^n = j · 2^n / 2^k = j / 2^{k-n}.
# For n ≥ k, ξ·2^n is integer → factor = 1.
# For n < k, factor = (1-p) + p · exp(2πi · j / 2^{k-n}).
# ============================================================
def bc_fourier(j, k, p, max_n=64):
    """ν̂_p(j/2^k) for product measure with bits ~ Ber(p)."""
    if j == 0:
        return 1.0 + 0j
    val = 1.0 + 0j
    for n in range(k):
        # contribution: (1-p) + p · exp(2πi · j / 2^{k-n})
        denom = 1 << (k - n)
        # j may have factor of 2^a, simplify
        jn = j
        dn = denom
        # Reduce j/denom
        while jn % 2 == 0 and dn > 1:
            jn //= 2; dn //= 2
        if dn == 1:
            factor = 1.0 + 0j
        else:
            angle = 2 * math.pi * jn / dn
            factor = (1 - p) + p * (math.cos(angle) + 1j * math.sin(angle))
        val *= factor
    return val


# ============================================================
# Main analysis
# ============================================================
def main():
    print("# Fourier analysis of trajectory measure on Z_2")
    N = 1 << 22
    t0 = time.perf_counter()
    tree = build_inverse_tree(N)
    sizes = subtree_sizes(tree)
    odd_ints = np.array([m for m in tree if m & 1], dtype=np.int64)
    weights = np.array([float(sizes[m]) for m in odd_ints])
    print(f"# Tree built: {len(odd_ints)} odd nodes, time {time.perf_counter()-t0:.2f}s")

    K_MAX = 20
    print(f"\n# Computing μ̂(ξ) on dyadic grid for k=2..{K_MAX}")
    mu_hat_per_k = {}
    for k in range(2, K_MAX + 1):
        mh = compute_M_k_and_fft(odd_ints, weights, k)
        mu_hat_per_k[k] = mh
    print(f"# μ̂ computed, total time {time.perf_counter()-t0:.2f}s")

    # ==================
    # Step 3: Fit decay law σ
    # |μ̂(ξ)|² ~ |ξ|^(-σ) where |ξ| ~ 2^k for ξ = j/2^k with odd j
    # For each k, compute mean and median of |μ̂|² over odd j.
    # ==================
    print("\n# Step 3: Fourier decay")
    print(f"  {'k':>3}  {'|ξ|':>8}  {'#odd j':>8}  {'mean|μ̂|²':>14}  {'median|μ̂|²':>14}")
    decay_rows = []
    for k in range(3, K_MAX + 1):
        mh = mu_hat_per_k[k]
        # Odd j = j with lowest bit 1, i.e., j = 1, 3, 5, ... 2^k - 1
        odd_idx = np.arange(1, 1 << k, 2)
        vals = np.abs(mh[odd_idx]) ** 2
        mean_sq = vals.mean()
        median_sq = float(np.median(vals))
        xi_size = 2 ** k  # |ξ| = 1/(reduced denominator) — with ξ=j/2^k, j odd, |ξ| = 1/2^? actually we use |ξ| as 2^k for binning
        decay_rows.append((k, xi_size, len(odd_idx), mean_sq, median_sq))
        print(f"  {k:>3}  {xi_size:>8}  {len(odd_idx):>8}  {mean_sq:>14.6e}  {median_sq:>14.6e}")

    # Linear fit log(mean |μ̂|²) vs log(|ξ|) (= k log 2)
    ks_arr = np.array([r[0] for r in decay_rows], dtype=float)
    xi_logs = ks_arr * math.log(2)
    mean_sqs = np.array([r[3] for r in decay_rows])
    valid = mean_sqs > 0
    log_means = np.log(mean_sqs[valid])
    log_xis = xi_logs[valid]
    # Fit only large-k (asymptotic regime)
    fit_mask = ks_arr[valid] >= 10
    if fit_mask.sum() >= 3:
        slope, intercept = np.polyfit(log_xis[fit_mask], log_means[fit_mask], 1)
        sigma = -slope
        # 95% CI estimate
        residuals = log_means[fit_mask] - (slope * log_xis[fit_mask] + intercept)
        ssr = np.sum(residuals ** 2)
        n = fit_mask.sum()
        se = math.sqrt(ssr / (n - 2)) / math.sqrt(np.var(log_xis[fit_mask]) * n)
        ci = 1.96 * se
        print(f"\n  Fitted decay σ = {sigma:.4f} ± {ci:.4f} (95% CI)")
        print(f"  Reference D_q values (R61):")
        print(f"    D_0 = 1.00 → 2·D_0 = 2.00")
        print(f"    D_1 = 0.608 → 2·D_1 = 1.216")
        print(f"    D_2 = 0.267 → 2·D_2 = 0.534")
        print(f"    D_∞ ≈ 0.15 → 2·D_∞ = 0.30")
    else:
        print("  Insufficient data for asymptotic fit")
        sigma = float('nan'); ci = float('nan')

    # Step 4: which D_q matches σ?
    if not math.isnan(sigma):
        candidates = {'2·D_0=2.00': 2.00, '2·D_1=1.216': 1.216, '2·D_2=0.534': 0.534,
                      '2·D_∞=0.30': 0.30, 'D_1=0.608': 0.608, 'D_2=0.267': 0.267}
        diffs = {name: abs(sigma - val) for name, val in candidates.items()}
        best = min(diffs, key=diffs.get)
        print(f"\n  Closest D_q candidate: {best} (Δ = {diffs[best]:.4f})")

    # Save decay table
    out_csv = os.path.join(OUTDIR, "fourier_decay.csv")
    with open(out_csv, 'w') as f:
        f.write("k,xi_2k,n_odd_j,mean_abs_sq,median_abs_sq\n")
        for r in decay_rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.10e},{r[4]:.10e}\n")
    print(f"  [save] {out_csv}")

    # ==================
    # Step 5: Bernoulli-convolution match
    # Compare log|μ̂|² to log|ν̂_p|² for various p
    # ==================
    print("\n# Step 5: BC product-measure match (Ber(p) bits)")
    p_values = [0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60]
    # For each p, compute ν̂_p on a sample of dyadic frequencies
    # Sample subset: for k = 6, 8, 10, 12, 15, take all odd j (or random subset)
    sample_freqs = []
    for k in [6, 8, 10, 12, 15]:
        for j in range(1, 1 << k, 2):
            sample_freqs.append((j, k))
    print(f"  Sampling {len(sample_freqs)} dyadic frequencies")

    # Compute log|μ̂|² at sample freqs
    mu_sq_sample = np.array([abs(mu_hat_per_k[k][j]) ** 2 for j, k in sample_freqs])
    valid_sample = mu_sq_sample > 1e-30
    log_mu_sample = np.log(mu_sq_sample[valid_sample])

    print(f"  {'p':>6}  {'Pearson(log|μ̂|², log|ν̂_p|²)':>32}  {'mean log diff':>14}")
    bc_rows = []
    for p in p_values:
        log_bc_sample = np.array([math.log(abs(bc_fourier(j, k, p)) ** 2 + 1e-30) for j, k in sample_freqs])
        log_bc_valid = log_bc_sample[valid_sample]
        # Pearson
        mx = log_mu_sample.mean(); my = log_bc_valid.mean()
        num = np.sum((log_mu_sample - mx) * (log_bc_valid - my))
        den = math.sqrt(np.sum((log_mu_sample - mx) ** 2) * np.sum((log_bc_valid - my) ** 2))
        r = num / den if den > 0 else 0
        diff = (log_mu_sample - log_bc_valid).mean()
        bc_rows.append((p, r, diff))
        print(f"  {p:>6.2f}  {r:>+32.4f}  {diff:>14.4f}")

    # Save BC results
    out_bc = os.path.join(OUTDIR, "bernoulli_convolution_match.csv")
    with open(out_bc, 'w') as f:
        f.write("p,pearson_r,mean_log_diff\n")
        for r in bc_rows:
            f.write(f"{r[0]},{r[1]:.6f},{r[2]:.6f}\n")
    print(f"  [save] {out_bc}")

    # ==================
    # Step 7: Multifractal Fourier σ_q
    # σ_q = limit as |ξ|→∞ of (-2/q) log E[|μ̂(ξ)|^q] / log |ξ|
    # Here we compute (1/k) · (1/q) · log E[|μ̂(j/2^k)|^q] across odd j
    # then take asymptotic slope
    # ==================
    print("\n# Step 7: Multifractal Fourier σ_q")
    qs = [0.5, 1.0, 2.0, 3.0, 5.0]
    print(f"  {'q':>5}  {'σ_q':>10}  {'ref 2·D_q':>10}")
    sigma_q_rows = []
    for q in qs:
        # E[|μ̂|^q] per k, log-fit slope
        log_moments = []
        ks_used = []
        for k in range(8, K_MAX + 1):
            mh = mu_hat_per_k[k]
            odd_idx = np.arange(1, 1 << k, 2)
            absvals = np.abs(mh[odd_idx])
            absvals = absvals[absvals > 1e-50]
            if len(absvals) == 0:
                continue
            moment = (absvals ** q).mean()
            log_moments.append(math.log(moment) / q)  # so it's like log(L^q-norm)
            ks_used.append(k)
        if len(log_moments) >= 3:
            xs = np.array(ks_used, dtype=float) * math.log(2)
            ys = np.array(log_moments)
            slope, _ = np.polyfit(xs, ys, 1)
            sigma_q = -2 * slope  # σ_q in convention |μ̂|^q ~ |ξ|^(-q σ_q / 2)
            # (or some convention — adjust below)
            sigma_q_rows.append((q, sigma_q))
            print(f"  {q:>5.1f}  {sigma_q:>10.4f}")

    # ==================
    # Step 8: Resonance / self-similarity scan
    # ==================
    print("\n# Step 8: Resonance peaks and self-similarity")
    # Look for high |μ̂|² at specific frequencies
    print("\n  Top 10 high-magnitude frequencies at k=10:")
    k = 10
    mh = mu_hat_per_k[k]
    odd_idx = np.arange(1, 1 << k, 2)
    abs_sq = np.abs(mh[odd_idx]) ** 2
    top_indices = np.argsort(abs_sq)[-10:][::-1]
    print(f"    {'rank':>4}  {'j':>5}  {'j/2^k':>10}  {'|μ̂|²':>14}")
    for rank, idx_in_odd in enumerate(top_indices):
        j = odd_idx[idx_in_odd]
        print(f"    {rank+1:>4}  {j:>5}  {j/(2**k):>10.5f}  {abs_sq[idx_in_odd]:>14.6e}")

    # Self-similarity check: μ̂(2ξ) vs μ̂(ξ) — but since ξ is dyadic, doubling shifts to lower-k freq
    # μ̂(j/2^k) vs μ̂(j/2^{k-1}) where j is odd in second case (so j → j*2 in original numbering)
    # Compare ratios for various j, k
    print("\n  Self-similarity ratio μ̂(2ξ) / μ̂(ξ) at k=10:")
    if k >= 2:
        mh_lower = mu_hat_per_k[k-1]
        sample_js = [1, 3, 5, 7, 11, 13, 21, 31]
        print(f"    {'j':>3}  {'ξ=j/2^10':>10}  {'2ξ=j/2^9':>10}  {'|μ̂(ξ)|':>10}  {'|μ̂(2ξ)|':>10}  {'ratio':>10}")
        for j in sample_js:
            if j < (1 << (k-1)):
                a = abs(mu_hat_per_k[k][j])
                b = abs(mu_hat_per_k[k-1][j]) if j < len(mu_hat_per_k[k-1]) else float('nan')
                ratio = b / a if a > 0 else float('nan')
                print(f"    {j:>3}  {j/(2**k):>10.6f}  {j/(2**(k-1)):>10.6f}  {a:>10.6e}  {b:>10.6e}  {ratio:>10.4f}")

    # ==================
    # Step 9: Chang π Fourier comparison
    # ==================
    print("\n# Step 9: Chang π Fourier comparison")
    # Chang's π values are mod-64 cylinder weights. Equivalently a measure on Z/64Z, lifted to Z_2.
    # Build "Chang weights": w_chang(m) = π_chang(m mod 64) for m odd.
    chang_pi = {}
    with open(r"C:\Collatz\experiments_output\chang_pi.csv") as f:
        for row in csv.DictReader(f):
            chang_pi[int(row['r'])] = float(row['pi_float'])
    chang_weights = np.array([chang_pi.get(m % 64, 0.0) for m in odd_ints])

    # Compute Chang-π Fourier at same dyadic grid (only need to bin into M_k mod 2^k)
    chang_mu_per_k = {}
    for k in range(2, K_MAX + 1):
        chang_mu_per_k[k] = compute_M_k_and_fft(odd_ints, chang_weights, k)

    print(f"  {'k':>3}  {'mean|μ̂_traj|²':>14}  {'mean|μ̂_chang|²':>14}  {'ratio':>10}")
    chang_rows = []
    for k in range(3, K_MAX + 1):
        odd_idx = np.arange(1, 1 << k, 2)
        traj = np.abs(mu_hat_per_k[k][odd_idx]) ** 2
        chang = np.abs(chang_mu_per_k[k][odd_idx]) ** 2
        traj_m = traj.mean(); chang_m = chang.mean()
        ratio = traj_m / chang_m if chang_m > 0 else float('nan')
        chang_rows.append((k, traj_m, chang_m, ratio))
        print(f"  {k:>3}  {traj_m:>14.6e}  {chang_m:>14.6e}  {ratio:>10.4f}")

    # Find frequencies where μ̂_traj and μ̂_chang differ MOST (discriminating frequencies)
    print("\n  Most discriminating frequencies at k=10 (largest |traj - chang| / chang):")
    k = 10
    mh_t = mu_hat_per_k[k]
    mh_c = chang_mu_per_k[k]
    odd_idx = np.arange(1, 1 << k, 2)
    diff = np.abs(np.abs(mh_t[odd_idx]) - np.abs(mh_c[odd_idx]))
    top = np.argsort(diff)[-10:][::-1]
    print(f"    {'rank':>4}  {'j':>5}  {'ξ':>10}  {'|μ̂_traj|':>12}  {'|μ̂_chang|':>12}")
    for rank, idx in enumerate(top):
        j = odd_idx[idx]
        print(f"    {rank+1:>4}  {j:>5}  {j/(2**k):>10.5f}  {abs(mh_t[j]):>12.6e}  {abs(mh_c[j]):>12.6e}")

    # Save grid summary
    out_grid = os.path.join(OUTDIR, "mu_hat_summary.csv")
    with open(out_grid, 'w') as f:
        f.write("k,xi_size_2k,mean_abs_sq_traj,median_abs_sq_traj,mean_abs_sq_chang,ratio\n")
        for k_row, dr_row, ch_row in zip(decay_rows, decay_rows, chang_rows):
            f.write(f"{k_row[0]},{k_row[1]},{k_row[3]:.6e},{k_row[4]:.6e},{ch_row[2]:.6e},{ch_row[3]:.4f}\n")
    print(f"\n[save] {out_grid}")
    print(f"\nTotal time: {time.perf_counter()-t0:.2f}s")


if __name__ == "__main__":
    main()
