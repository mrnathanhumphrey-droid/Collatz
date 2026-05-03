"""
stopping_criteria_battery.py

Walk Collatz orbits, record sigma under multiple stopping criteria:
  - sigma_to_1    : standard, steps to reach 1
  - sigma_prime   : first time orbit value is prime
  - sigma_pow2    : first time orbit value is exact power of 2 with exponent >= L
  - sigma_logdrop_k : first time log(orbit) drops by >= k below log(start), k in {1,2,4,8}
  - sigma_res_r_p : first time orbit value ≡ r mod p, for several (r, p)

Outputs per-criterion stats, cross-correlations, per-sigma_to_1-band stratification.
N = 2^32, 200K orbits.
"""
import sys, time
from pathlib import Path
import numpy as np
import polars as pl
from numba import njit, prange

OUT = Path(r"C:\Collatz")
sys.stdout.reconfigure(encoding="utf-8")
log_lines = []
def log(s): print(s, flush=True); log_lines.append(s)


@njit(cache=True)
def mulmod(a, b, m):
    # Compute (a * b) mod m without overflow using Python-style int (numba uses int64).
    # For our orbits a, b < 2^62 and m < 2^62, so a*b can overflow int64.
    # Use Russian peasant multiplication mod m.
    result = 0
    a = a % m
    while b > 0:
        if b & 1:
            result = (result + a) % m
        a = (a + a) % m
        b >>= 1
    return result

@njit(cache=True)
def powmod(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp & 1:
            result = mulmod(result, base, mod)
        exp >>= 1
        base = mulmod(base, base, mod)
    return result

# Deterministic Miller-Rabin for n < 3.3e24 with these witnesses
@njit(cache=True)
def is_prime(n):
    if n < 2:
        return False
    SMALL = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in SMALL:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in SMALL:
        if a >= n:
            continue
        x = powmod(a, d, n)
        if x == 1 or x == n - 1:
            continue
        composite = True
        for _ in range(r - 1):
            x = mulmod(x, x, n)
            if x == n - 1:
                composite = False
                break
        if composite:
            return False
    return True


@njit(cache=True)
def is_pow2(n):
    return n > 0 and (n & (n - 1)) == 0


L_POW2 = 10  # threshold for sigma_pow2 (require value >= 2^10 = 1024)
LOG_DROPS = np.array([1.0, 2.0, 4.0, 8.0])
RESIDUE_TARGETS = np.array([
    [0, 3], [1, 8], [5, 8], [0, 5]
], dtype=np.int64)


@njit(parallel=True, cache=True)
def walk_with_criteria(starts, max_steps, max_value):
    n = len(starts)
    n_logdrops = len(LOG_DROPS)
    n_residues = len(RESIDUE_TARGETS)

    sigma_to_1 = np.full(n, -1, dtype=np.int64)
    sigma_prime = np.full(n, -1, dtype=np.int64)
    sigma_pow2 = np.full(n, -1, dtype=np.int64)
    sigma_logdrop = np.full((n, n_logdrops), -1, dtype=np.int64)
    sigma_residue = np.full((n, n_residues), -1, dtype=np.int64)

    for i in prange(n):
        start = starts[i]
        log_start = np.log(np.float64(start))
        x = start
        steps = 0
        prime_done = False
        pow2_done = False
        logdrop_done = np.zeros(n_logdrops, dtype=np.bool_)
        residue_done = np.zeros(n_residues, dtype=np.bool_)
        n_logdrop_done = 0
        n_residue_done = 0
        # Skip checking start itself for "first hit" criteria? Brief says "first reach";
        # start qualifies if start itself meets criterion (e.g., start ≡ 1 mod 8).
        # We'll allow steps >= 0; sigma=0 if start meets it.
        while True:
            # Check criteria
            if not prime_done and is_prime(x):
                sigma_prime[i] = steps
                prime_done = True
            if not pow2_done and is_pow2(x) and x >= (1 << L_POW2):
                sigma_pow2[i] = steps
                pow2_done = True
            if n_logdrop_done < n_logdrops:
                log_x = np.log(np.float64(x))
                drop = log_start - log_x
                for j in range(n_logdrops):
                    if not logdrop_done[j] and drop >= LOG_DROPS[j]:
                        sigma_logdrop[i, j] = steps
                        logdrop_done[j] = True
                        n_logdrop_done += 1
            if n_residue_done < n_residues:
                for j in range(n_residues):
                    if not residue_done[j]:
                        if x % RESIDUE_TARGETS[j, 1] == RESIDUE_TARGETS[j, 0]:
                            sigma_residue[i, j] = steps
                            residue_done[j] = True
                            n_residue_done += 1

            if x == 1:
                sigma_to_1[i] = steps
                break
            if steps >= max_steps or x > max_value:
                break

            if x & 1 == 0:
                x = x >> 1
            else:
                x = 3 * x + 1
            steps += 1
    return sigma_to_1, sigma_prime, sigma_pow2, sigma_logdrop, sigma_residue


def main():
    log("=" * 78)
    log("Stopping criteria battery — Collatz q=3, N=2^32, 200K orbits")
    log("=" * 78)

    log2N = 32
    N = 1 << log2N
    n_orbits = 200_000
    rng = np.random.default_rng(42)
    starts = 2 * rng.integers(1, (N - 1) // 2, size=n_orbits, dtype=np.int64) + 1

    log(f"\n[walk] {n_orbits:,} orbits, max_value=2^62 = {1<<62:.2e}")
    t0 = time.perf_counter()
    s1, sp, s2, sld, sres = walk_with_criteria(starts, 10_000_000, np.int64(1 << 62))
    log(f"[walk] done in {time.perf_counter()-t0:.1f}s")

    # Filter orbits that completed (sigma_to_1 >= 0)
    mask = s1 >= 0
    n_ok = int(mask.sum())
    log(f"[walk] {n_ok:,} / {n_orbits:,} orbits reached 1")

    starts_ok = starts[mask]
    s1 = s1[mask]; sp = sp[mask]; s2 = s2[mask]
    sld = sld[mask]; sres = sres[mask]

    # Per-criterion firing rate
    log(f"\n=== Step 1: Per-criterion firing rates ===")
    crit_data = {
        "sigma_to_1": s1,
        "sigma_prime": sp,
        f"sigma_pow2_L{L_POW2}": s2,
    }
    for j, k in enumerate(LOG_DROPS):
        crit_data[f"sigma_logdrop_{int(k)}"] = sld[:, j]
    for j, (r, p) in enumerate(RESIDUE_TARGETS):
        crit_data[f"sigma_res_{r}mod{p}"] = sres[:, j]

    log(f"  {'criterion':<24} {'fired':>9} {'fire_rate':>10} {'mean':>10} {'sd':>10} {'p50':>9} {'p95':>9}")
    for name, arr in crit_data.items():
        fired = arr >= 0
        n_f = int(fired.sum())
        if n_f > 0:
            vals = arr[fired].astype(np.float64)
            log(f"  {name:<24} {n_f:>9} {n_f/n_ok:>10.3%} "
                f"{vals.mean():>10.2f} {vals.std():>10.2f} "
                f"{np.percentile(vals, 50):>9.0f} {np.percentile(vals, 95):>9.0f}")
        else:
            log(f"  {name:<24} {n_f:>9} {0.0:>10.3%} {'--':>10} {'--':>10} {'--':>9} {'--':>9}")

    # Cross-correlations (only on orbits where ALL criteria fired)
    log(f"\n=== Step 4: Cross-criterion correlations (on orbits where ALL criteria fired) ===")
    all_fire = np.ones(n_ok, dtype=bool)
    for arr in crit_data.values():
        all_fire &= (arr >= 0)
    n_all = int(all_fire.sum())
    log(f"  {n_all:,} orbits where every criterion fired")

    if n_all > 100:
        names = list(crit_data.keys())
        M = np.column_stack([crit_data[n][all_fire].astype(np.float64) for n in names])
        corr = np.corrcoef(M.T)
        log(f"\n  Pearson correlation matrix:")
        log(f"  {'':<24} " + " ".join(f"{n[:11]:>11}" for n in names))
        for i, n in enumerate(names):
            log(f"  {n:<24} " + " ".join(f"{corr[i, j]:>+11.4f}" for j in range(len(names))))

    # Per-σ_to_1 band stratification
    log(f"\n=== Step 3: Per-σ_to_1-band stratification ===")
    log(f"  Bands by σ_to_1 quartile, plus q5 = top 5%")
    qs = np.percentile(s1, [25, 50, 75, 95]).astype(np.float64)
    log(f"  Band thresholds (σ_to_1): q25={qs[0]:.0f}, q50={qs[1]:.0f}, q75={qs[2]:.0f}, q95={qs[3]:.0f}")

    bands = np.zeros(n_ok, dtype=np.int8)
    bands[s1 > qs[0]] = 1
    bands[s1 > qs[1]] = 2
    bands[s1 > qs[2]] = 3
    band_q5 = s1 > qs[3]

    # For each criterion, compute mean ratio σ_crit / σ_to_1 per band
    log(f"\n  ratio = mean(σ_criterion / σ_to_1) per band; lower = criterion fires earlier")
    log(f"  {'criterion':<24}  {'q1 (low)':>9} {'q2':>9} {'q3':>9} {'q4 (hi)':>9} {'q5 (95+)':>9}")
    for name, arr in crit_data.items():
        if name == "sigma_to_1":
            continue
        ratios = []
        for b in range(4):
            in_band = (bands == b) & (arr >= 0)
            if in_band.sum() > 50:
                r = float(np.mean(arr[in_band].astype(np.float64) / s1[in_band].astype(np.float64).clip(min=1)))
                ratios.append(r)
            else:
                ratios.append(float("nan"))
        in_q5 = band_q5 & (arr >= 0)
        if in_q5.sum() > 50:
            r5 = float(np.mean(arr[in_q5].astype(np.float64) / s1[in_q5].astype(np.float64).clip(min=1)))
        else:
            r5 = float("nan")
        log(f"  {name:<24}  " + " ".join(f"{r:>9.4f}" for r in ratios) + f" {r5:>9.4f}")

    # Save
    out_rows = []
    for name, arr in crit_data.items():
        fired = arr >= 0
        if fired.sum() > 0:
            vals = arr[fired].astype(np.float64)
            out_rows.append({"criterion": name, "fire_rate": float(fired.sum()/n_ok),
                             "mean": float(vals.mean()), "sd": float(vals.std()),
                             "p50": float(np.percentile(vals, 50)),
                             "p95": float(np.percentile(vals, 95))})
    pl.DataFrame(out_rows).write_csv(OUT / "stopping_criteria_battery.csv")
    (OUT / "stopping_criteria_battery_log.txt").write_text("\n".join(log_lines), encoding="utf-8")
    log(f"\n[wrote] stopping_criteria_battery.csv, stopping_criteria_battery_log.txt")


if __name__ == "__main__":
    main()
