"""
result_77_7_v2_k8.py — extend R77.7 v2 to k=8.

Loads cached eps_1..eps_7 from result_77_7_eps_exact_through_k7_v2.json,
runs compute_pi_k_exact(8) with the CRT+modular solver, computes eps_8.

Wall-time estimate: ~14.4s/prime × 15.5× (N^3 ratio reduced by numpy
vectorization, k=6→k=7 was 15.5× not theoretical 27×) at k=7→k=8 →
~220s/prime. With ~400 primes for denominator growth: ~24hr wall.
Saves incrementally to result_77_7_eps_exact_through_k8_v2.json.
"""

import json
import os
import sys
import time
from fractions import Fraction

# Reuse from result_77_7_v2.py
sys.path.insert(0, r"C:\Collatz")
from result_77_7_v2 import (
    compute_pi_k_exact,
    eps_from_pi_sequence,
    EPS_CACHE_V2,
)

EPS_CACHE_K8 = r"C:\Collatz\experiments_output\result_77_7_eps_exact_through_k8_v2.json"


def load_cached_eps_v2():
    if not os.path.exists(EPS_CACHE_V2):
        raise RuntimeError(f"Cache not found: {EPS_CACHE_V2}. Run result_77_7_v2.py first.")
    with open(EPS_CACHE_V2, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return {int(k): Fraction(int(d["num"]), int(d["den"])) for k, d in data.items()}


def save_eps_cache_k8(eps_dict):
    out = {str(k): {"num": str(v.numerator), "den": str(v.denominator)}
           for k, v in eps_dict.items()}
    with open(EPS_CACHE_K8, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)


def main():
    print("=" * 78)
    print("R77.7 V2 K=8 extension — exact-rational eps_8 via CRT+modular solver")
    print("=" * 78)
    print(f"Output cache: {EPS_CACHE_K8}")
    print()

    # Load v1 cached eps_1..eps_7
    cached_eps = load_cached_eps_v2()
    print(f"Loaded {len(cached_eps)} cached eps_k from v2:")
    for k in sorted(cached_eps.keys()):
        print(f"  eps_{k} ~ {float(cached_eps[k]):+.10e}")
    print()

    if 7 not in cached_eps:
        raise RuntimeError("eps_7 not in v2 cache; cannot proceed to k=8.")

    # Initial-primes estimate for k=8:
    # pi_8 denominators ~3.3x bigger than pi_7's. k=7 used 150 primes for ~300-digit
    # denominators; k=8 with ~1000-digit denominators needs ~500 primes.
    # The compute_pi_k_exact loop auto-grows by 10 primes if reconstruction fails.
    n_primes_initial = 500

    print("=" * 78)
    print(f"Phase 5+: compute pi_8 with n_primes_initial = {n_primes_initial}")
    print(f"N=4374 dense Gauss; estimated per-prime ~220s, total ~24hr")
    print("=" * 78)

    t_start = time.time()
    pi_8, coprime_8, t_solves_8, np_used_8 = compute_pi_k_exact(
        8, n_primes_initial=n_primes_initial)
    t_total = time.time() - t_start

    # Compute eps_8 from eps_1..eps_7 (cached) + pi_8
    # eps_from_pi_sequence needs the full pi_dict; we only have pi_8.
    # Use the formula directly: X_8 = 3^8 * sum(pi_8^2), X_7 = 1 + 6*(7/15) + sum(eps_1..6)
    target = Fraction(7, 15)
    X_8 = Fraction(3 ** 8) * sum(p * p for p in pi_8)
    sum_eps_1_through_6 = sum(cached_eps[j] for j in range(1, 7))
    X_7 = Fraction(1) + Fraction(6) * target + sum_eps_1_through_6
    S_8 = X_8 - X_7
    eps_8 = S_8 - target

    # Update eps dict
    eps_dict = dict(cached_eps)
    eps_dict[8] = eps_8

    # Save
    save_eps_cache_k8(eps_dict)

    print()
    print(f"[k=8] TOTAL wall time {t_total:.1f}s = {t_total/60:.2f} min = {t_total/3600:.2f} hr")
    print(f"[k=8] primes used {np_used_8}")
    print(f"[k=8] mean per-prime solve {sum(t_solves_8)/len(t_solves_8):.2f}s")
    print()
    print(f"eps_8 = {eps_8.numerator} / {eps_8.denominator}")
    print(f"eps_8 (float) ~ {float(eps_8):+.15e}")
    print(f"|eps_8| * 2^8 = {abs(float(eps_8)) * 256:+.10f}")
    print()
    print("|eps_n|·2^n series:")
    for k in sorted(eps_dict.keys()):
        e = eps_dict[k]
        print(f"  n={k}: |eps_{k}|*2^{k} = {abs(float(e)) * (2 ** k):+.10f}")
    print()
    print(f"Saved eps_1..eps_8 to {EPS_CACHE_K8}")


if __name__ == "__main__":
    main()
