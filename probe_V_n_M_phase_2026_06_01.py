"""
V_n^M phase-twisted moment computation at q=3, n=2..5.

M_n_phase(eta, delta) = sum_{xi in (Z/3^n)*} exp(-2pi i xi delta / 3^n) muhat_n(xi) muhat_n*(xi eta)

For each level n: compute M_n_phase as a complex-valued tensor over (eta, delta).
Then project to a COMMON level (we use level n=3, dim 27 x 27 = 729 entries but only
18 x 27 = 486 with eta coprime). Compare across n.

Hypothesis: the dominant subdominant eigenvalue (predicted Cotaescu fit: rho=0.19,
theta=41°) controls how fast M_n_phase converges to M_inf_phase, as a function of
(eta, delta). Extract by per-direction decay-rate fit.
"""
from __future__ import annotations
import sys, os, time, json, cmath
from fractions import Fraction
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:/Collatz")

from T_lead_operator import build_markov_rational, stationary_rational

N_MAX = 5
COMMON_LEVEL = 3  # we'll project all M_n_phase down to level 3 indexing

# Build pi_n for n=2..N_MAX
print(f"[1/4] Compute pi_n stationary distributions at q=3, n=2..{N_MAX}", flush=True)
t0 = time.time()
pi_data = {}
for n in range(2, N_MAX + 1):
    tn = time.time()
    K, coprime = build_markov_rational(n)
    pi_q = stationary_rational(K)
    pi_data[n] = {"pi": pi_q, "coprime": coprime}
    print(f"  n={n}: |coprime|={len(coprime)}, time={time.time()-tn:.1f}s", flush=True)

# Compute muhat_n(xi) for all xi in Z/3^n
print(f"\n[2/4] Compute muhat_n(xi) for xi in Z/3^n", flush=True)
muhat = {}
for n in range(2, N_MAX + 1):
    tn = time.time()
    N = 3**n
    pi_q = pi_data[n]["pi"]
    coprime = pi_data[n]["coprime"]
    pi_dict = {coprime[i]: float(pi_q[i]) for i in range(len(coprime))}
    muh = np.zeros(N, dtype=np.complex128)
    for xi in range(N):
        for r in coprime:
            muh[xi] += pi_dict[r] * cmath.exp(-2j * cmath.pi * r * xi / N)
    muhat[n] = muh
    print(f"  n={n}: muhat computed, time={time.time()-tn:.1f}s", flush=True)

# Compute M_n_phase(eta, delta) for (eta, delta) in (Z/3^n)* x (Z/3^n)
print(f"\n[3/4] Compute M_n_phase(eta, delta) at each level", flush=True)
M_phase = {}
for n in range(2, N_MAX + 1):
    tn = time.time()
    N = 3**n
    coprime = pi_data[n]["coprime"]
    muh = muhat[n]
    # Build phase factor table phase[xi, delta] = exp(-2pi i xi delta / 3^n)
    # Output shape: M[eta_idx, delta] of dimension (len(coprime), N)
    M = np.zeros((len(coprime), N), dtype=np.complex128)
    for eta_idx, eta in enumerate(coprime):
        # Precompute muhat(xi*eta mod N) for all xi
        muh_eta = np.zeros(N, dtype=np.complex128)
        for xi in range(N):
            muh_eta[xi] = muh[(xi * eta) % N]
        for delta in range(N):
            # M(eta, delta) = sum_{xi in (Z/3^n)*} exp(-2pi i xi delta / 3^n) muhat(xi) muhat*(xi eta)
            val = 0.0 + 0.0j
            for xi in coprime:
                val += cmath.exp(-2j * cmath.pi * xi * delta / N) * muh[xi] * np.conj(muh_eta[xi])
            M[eta_idx, delta] = val
    M_phase[n] = M
    print(f"  n={n}: M_phase shape={M.shape}, time={time.time()-tn:.1f}s", flush=True)

# Project all M_n_phase down to level COMMON_LEVEL
# Indexing: for each (eta, delta) at level COMMON, sum contributions from level n where
#   (eta_n mod 3^COMMON, delta_n mod 3^COMMON) = (eta_common, delta_common)
print(f"\n[4/4] Project M_n_phase to COMMON_LEVEL={COMMON_LEVEL}", flush=True)
N_common = 3 ** COMMON_LEVEL
coprime_common = pi_data[COMMON_LEVEL]["coprime"]
common_to_idx = {r: i for i, r in enumerate(coprime_common)}

# Projected: M_proj[n, eta_common_idx, delta_common]
M_proj = {}
for n in range(2, N_MAX + 1):
    N = 3 ** n
    coprime_n = pi_data[n]["coprime"]
    Mn = M_phase[n]
    Mp = np.zeros((len(coprime_common), N_common), dtype=np.complex128)
    counts = np.zeros((len(coprime_common), N_common), dtype=np.int64)
    for eta_idx, eta in enumerate(coprime_n):
        eta_red = eta % N_common
        if eta_red == 0 or eta_red % 3 == 0:
            continue
        if eta_red not in common_to_idx:
            continue
        ec = common_to_idx[eta_red]
        for delta in range(N):
            dc = delta % N_common
            Mp[ec, dc] += Mn[eta_idx, delta]
            counts[ec, dc] += 1
    # Average per (ec, dc)
    Mp = np.divide(Mp, np.maximum(counts, 1), out=np.zeros_like(Mp), where=(counts > 0))
    M_proj[n] = Mp
    nonz = np.count_nonzero(counts)
    print(f"  n={n}: projected, {nonz} cells filled, time-cumulative")

# Treat M_proj[5] as best approximation to M_inf (asymptote)
# Per (eta_common, delta_common), look at decay (M_n - M_5) across n=2,3,4 -> n=5
print(f"\n=== Decay rate per (eta, delta) direction ===")
print(f"   Using M_proj[n={N_MAX}] as approximation to M_inf")
print()

# For each (ec, dc), find decay rate fit
results_per_direction = []
M_inf_approx = M_proj[N_MAX]
n_values = list(range(2, N_MAX))  # exclude N_MAX since that's our M_inf
for ec in range(len(coprime_common)):
    for dc in range(N_common):
        diffs = []
        for n in n_values:
            d = M_proj[n][ec, dc] - M_inf_approx[ec, dc]
            diffs.append(d)
        diffs = np.array(diffs)
        if np.max(np.abs(diffs)) < 1e-12:
            continue
        # Quick rate estimate from |diff_2| / |diff_3| ratio
        ratios = []
        for i in range(len(diffs) - 1):
            if abs(diffs[i]) > 0:
                ratios.append(abs(diffs[i+1]) / abs(diffs[i]))
        if not ratios:
            continue
        rate = np.mean(ratios)
        # Phase change between consecutive n (angle of diffs)
        if all(abs(d) > 1e-12 for d in diffs):
            angles = [np.angle(diffs[i+1]) - np.angle(diffs[i]) for i in range(len(diffs)-1)]
            angle_mean = np.mean(angles)
        else:
            angle_mean = float('nan')
        eta_v = coprime_common[ec]
        results_per_direction.append({
            "eta": int(eta_v), "delta": int(dc),
            "magnitudes": [float(abs(d)) for d in diffs],
            "rate_estimate": float(rate),
            "angle_per_step_deg": float(angle_mean * 180 / np.pi) if not np.isnan(angle_mean) else None,
        })

# Sort by largest magnitude (most informative directions)
results_per_direction.sort(key=lambda r: -r["magnitudes"][0])
print(f"Top 20 (eta, delta) directions by initial magnitude |M_2 - M_inf|:")
print(f"  {'eta':>5} {'delta':>5}  {'|M_2-M_inf|':>14} {'|M_3-M_inf|':>14} {'|M_4-M_inf|':>14}  {'rate':>8} {'angle°':>10}")
for r in results_per_direction[:20]:
    mags_str = "  ".join(f"{m:>12.4e}" for m in r["magnitudes"])
    angle_str = f"{r['angle_per_step_deg']:+8.2f}" if r["angle_per_step_deg"] is not None else "  -- "
    print(f"  {r['eta']:>5} {r['delta']:>5}  {mags_str}  {r['rate_estimate']:>8.4f} {angle_str:>10}")
print()

# Aggregate: rate distribution across all directions
rates = [r["rate_estimate"] for r in results_per_direction if r["rate_estimate"] > 0]
angles = [r["angle_per_step_deg"] for r in results_per_direction if r["angle_per_step_deg"] is not None]
print(f"Rate distribution across {len(rates)} (eta, delta) directions:")
print(f"  min={min(rates):.4f}, max={max(rates):.4f}, mean={np.mean(rates):.4f}, median={np.median(rates):.4f}")
if angles:
    print(f"Angle distribution:")
    print(f"  min={min(angles):+.2f}°, max={max(angles):+.2f}°, mean={np.mean(angles):+.2f}°, median={np.median(angles):+.2f}°")

# Histogram of (rate, angle) pairs — look for cluster near (0.5, 41°) which would
# match Cotaescu fit (rate=ρ=0.19 for amplitude, but ratio of consecutive |ε| differs)
print(f"\nDirections within target Cotaescu zone (rate ∈ [0.40, 0.55], angle ∈ [35°, 50°]):")
hits = [r for r in results_per_direction if 0.40 <= r["rate_estimate"] <= 0.55 and r["angle_per_step_deg"] is not None and 35 <= abs(r["angle_per_step_deg"]) <= 50]
for r in hits[:20]:
    print(f"  eta={r['eta']:>3}, delta={r['delta']:>3}: rate={r['rate_estimate']:.4f}, angle={r['angle_per_step_deg']:+.2f}°")
print(f"Total hits in zone: {len(hits)} of {len(results_per_direction)}")

# Save
output = {
    "N_MAX": N_MAX,
    "COMMON_LEVEL": COMMON_LEVEL,
    "results_per_direction": results_per_direction[:200],
    "rate_summary": {
        "min": float(min(rates)) if rates else None,
        "max": float(max(rates)) if rates else None,
        "mean": float(np.mean(rates)) if rates else None,
        "median": float(np.median(rates)) if rates else None,
    },
    "angle_summary": {
        "min": float(min(angles)) if angles else None,
        "max": float(max(angles)) if angles else None,
        "mean": float(np.mean(angles)) if angles else None,
        "median": float(np.median(angles)) if angles else None,
    } if angles else None,
    "cotaescu_zone_hits": len(hits),
    "total_directions": len(results_per_direction),
}
with open("C:/Collatz/V_n_M_phase_2026_06_01.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved: C:/Collatz/V_n_M_phase_2026_06_01.json")
print(f"Total wall time: {time.time()-t0:.1f}s")
