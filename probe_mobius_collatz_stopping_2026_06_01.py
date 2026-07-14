"""
Reddit Mobius-of-Collatz-stopping-time partial sum.

Y(m) = sum_{k=1..m} mu(S(k))  where S = Collatz stopping time, mu = Mobius.

Goal: extract envelope exponent alpha and Hurst H, compare to R77->R79
subdominant-rate evidence (conjectured 1/2, k=7,8 break suggests ~0.65-0.70).
"""
from __future__ import annotations
import sys, time, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress
from numba import njit

sys.stdout.reconfigure(encoding="utf-8")

N = 10**8

@njit(cache=True)
def stopping_times(N):
    S = np.zeros(N+1, dtype=np.int32)
    for k in range(2, N+1):
        n = np.int64(k)
        steps = 0
        while n >= k or n > N:
            if n & 1:
                n = 3*n + 1
            else:
                n = n >> 1
            steps += 1
        S[k] = S[n] + steps
    return S

def mobius_sieve(M):
    mu = np.zeros(M+1, dtype=np.int8)
    mu[1] = 1
    primes = []
    is_comp = np.zeros(M+1, dtype=bool)
    for i in range(2, M+1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            ip = i * p
            if ip > M:
                break
            is_comp[ip] = True
            if i % p == 0:
                mu[ip] = 0
                break
            else:
                mu[ip] = -mu[i]
    return mu

t0 = time.time()
print(f"[1/4] Computing stopping times for k=1..{N}...", flush=True)
S = stopping_times(N)
print(f"      done in {time.time()-t0:.1f}s", flush=True)

S_view = S[1:]
S_max = int(S_view.max())
S_mean = float(S_view.mean())
S_std  = float(S_view.std())
n_distinct = int(np.unique(S_view).size)
print(f"      S_max={S_max}, mean={S_mean:.1f}, std={S_std:.1f}, distinct={n_distinct}", flush=True)

print(f"[2/4] Sieving Mobius up to {S_max}...", flush=True)
mu = mobius_sieve(S_max)

print(f"[3/4] Building Y(m)...", flush=True)
muS = mu[S].astype(np.int32)
muS[0] = 0
Y = np.cumsum(muS).astype(np.int64)
print(f"      Y(N) = {int(Y[-1])}", flush=True)
abs_max = int(np.abs(Y).max())
abs_argmax = int(np.abs(Y).argmax())
print(f"      |Y|_max = {abs_max} at m = {abs_argmax}", flush=True)

print(f"[4/4] Extracting exponents...", flush=True)
absY = np.abs(Y).astype(np.float64)
running_max = np.maximum.accumulate(absY)

m_points = np.unique(np.logspace(3, np.log10(N), 300).astype(np.int64))
m_points = m_points[(m_points > 0) & (m_points <= N)]
rmax_at_m = running_max[m_points]
mask = rmax_at_m > 0
log_m = np.log10(m_points[mask].astype(np.float64))
log_rm = np.log10(rmax_at_m[mask])

slope_all, intercept, r_all, _, stderr_all = linregress(log_m, log_rm)
mid = len(log_m) // 2
slope_tail, intercept_tail, r_tail, _, stderr_tail = linregress(log_m[mid:], log_rm[mid:])

print(f"      Envelope alpha (full):  {slope_all:.4f} +/- {stderr_all:.4f}  r^2={r_all**2:.4f}", flush=True)
print(f"      Envelope alpha (tail):  {slope_tail:.4f} +/- {stderr_tail:.4f}  r^2={r_tail**2:.4f}", flush=True)

# DFA
print(f"      Running DFA...", flush=True)
chunk_sizes = np.unique(np.logspace(2, np.log10(N//20), 25).astype(np.int64))
chunk_sizes = chunk_sizes[chunk_sizes >= 100]
F_vals = []
cs_used = []
for n in chunk_sizes:
    n_chunks = N // n
    if n_chunks < 8:
        continue
    Yc = Y[:n_chunks*n].reshape(n_chunks, n).astype(np.float64)
    x = np.arange(n, dtype=np.float64)
    # vectorized linear detrend
    xm = x.mean()
    Sxx = ((x - xm)**2).sum()
    ym = Yc.mean(axis=1, keepdims=True)
    Sxy = ((x - xm) * (Yc - ym)).sum(axis=1)
    slope_c = Sxy / Sxx
    intercept_c = ym.squeeze() - slope_c * xm
    fit = slope_c[:, None] * x[None, :] + intercept_c[:, None]
    resid = Yc - fit
    F_chunk = np.sqrt((resid**2).mean(axis=1))
    F_vals.append(F_chunk.mean())
    cs_used.append(n)

F_vals = np.array(F_vals)
cs_used = np.array(cs_used)
mask_dfa = F_vals > 0
slope_dfa, intercept_dfa, r_dfa, _, stderr_dfa = linregress(np.log10(cs_used[mask_dfa]), np.log10(F_vals[mask_dfa]))
print(f"      DFA H:                  {slope_dfa:.4f} +/- {stderr_dfa:.4f}  r^2={r_dfa**2:.4f}", flush=True)

# Histogram
counts = np.bincount(S_view)
top_S_idx = np.argsort(counts)[::-1][:15]
print(f"\nTop-15 stopping times by frequency (k<={N}):")
print(f"  {'S':>5} {'count':>12} {'%':>7} {'mu(S)':>6}")
for s in top_S_idx:
    if counts[s] > 0:
        print(f"  {s:>5} {int(counts[s]):>12} {100*counts[s]/N:>6.2f}% {int(mu[s]):>+6d}")

# Plot
fig, ax = plt.subplots(2, 2, figsize=(14, 10))

ax[0,0].plot(np.arange(N+1), Y, linewidth=0.4, color="steelblue")
ax[0,0].set_xlabel("m")
ax[0,0].set_ylabel("Y(m)")
ax[0,0].set_title(f"Y(m) = sum_{{k=1..m}} mu(S(k)),  N={N:.0e}")
ax[0,0].grid(True, alpha=0.3)
ax[0,0].axhline(0, color="k", linewidth=0.5)

ax[0,1].loglog(m_points, rmax_at_m, "b-", linewidth=1, label="running max |Y|")
fit_y = 10**(intercept) * m_points.astype(np.float64)**slope_all
fit_y_tail = 10**(intercept_tail) * m_points.astype(np.float64)**slope_tail
ax[0,1].loglog(m_points, fit_y, "r--", linewidth=1, label=f"all: m^{slope_all:.3f}")
ax[0,1].loglog(m_points, fit_y_tail, "g--", linewidth=1, label=f"tail: m^{slope_tail:.3f}")
ref_x = m_points.astype(np.float64)
ref_y = ref_x**0.5
ref_y = ref_y * (rmax_at_m[0] / ref_y[0])
ax[0,1].loglog(ref_x, ref_y, color="gray", linestyle=":", linewidth=1, label="rand walk 1/2")
ax[0,1].set_xlabel("m")
ax[0,1].set_ylabel("max |Y|")
ax[0,1].set_title(f"Envelope exponent alpha")
ax[0,1].grid(True, which="both", alpha=0.3)
ax[0,1].legend(fontsize=9)

ax[1,0].loglog(cs_used[mask_dfa], F_vals[mask_dfa], "go-", linewidth=1, markersize=4, label="DFA F(n)")
fit_dfa = 10**intercept_dfa * cs_used[mask_dfa].astype(np.float64)**slope_dfa
ax[1,0].loglog(cs_used[mask_dfa], fit_dfa, "r--", linewidth=1, label=f"fit: n^{slope_dfa:.3f}")
ax[1,0].set_xlabel("chunk size n")
ax[1,0].set_ylabel("F(n)")
ax[1,0].set_title(f"DFA Hurst H = {slope_dfa:.4f}")
ax[1,0].grid(True, which="both", alpha=0.3)
ax[1,0].legend(fontsize=9)

S_show = min(200, len(counts))
ax[1,1].bar(range(S_show), counts[:S_show], color="navy", alpha=0.7)
ax[1,1].set_xlabel("S(k)")
ax[1,1].set_ylabel("count")
ax[1,1].set_title(f"Distribution of S(k):  distinct={n_distinct}, mean={S_mean:.0f}, std={S_std:.1f}")
ax[1,1].grid(True, alpha=0.3)

plt.tight_layout()
out_png = "C:/Collatz/mobius_collatz_Y_2026_06_01.png"
plt.savefig(out_png, dpi=120)
print(f"\nSaved plot to {out_png}", flush=True)

result = {
    "N": N,
    "S_max": S_max, "S_mean": S_mean, "S_std": S_std, "S_distinct": n_distinct,
    "Y_at_N": int(Y[-1]),
    "abs_Y_max": abs_max,
    "abs_Y_max_at": abs_argmax,
    "alpha_envelope_all": float(slope_all),
    "alpha_envelope_all_stderr": float(stderr_all),
    "alpha_envelope_tail": float(slope_tail),
    "alpha_envelope_tail_stderr": float(stderr_tail),
    "H_dfa": float(slope_dfa),
    "H_dfa_stderr": float(stderr_dfa),
    "top_S": [{"S": int(s), "count": int(counts[s]), "mu_S": int(mu[s])} for s in top_S_idx if counts[s] > 0],
}
out_json = "C:/Collatz/mobius_collatz_results_2026_06_01.json"
with open(out_json, "w") as f:
    json.dump(result, f, indent=2)
print(f"Saved results to {out_json}", flush=True)
print(f"\nTotal wall time: {time.time()-t0:.1f}s", flush=True)
