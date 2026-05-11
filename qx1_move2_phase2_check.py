"""
qx1_move2_phase2_check.py — Phase 2 empirical sanity for Move 2.

Two questions:

1. Does the generalized R78.3 formula |F̂_p(ξ)| = p^{(r+3)/2} hold at small
   primes p ∈ {5, 7} and small r? If yes, Route (a) of Phase 2 has a candidate
   structural derivation (Cochrane T2 + Plancherel + equidistribution, all
   p-blind).

2. Does K_p(r, c=1, m=0) saturate at √N = p^{(r-1)/2} at small primes? Even
   with R78.3 generalized, the K_p quantity is the SHORT-window character
   sum (length N = p^{r-1}) of a full-period (period p^r) signal. Its √N
   saturation is a different question from the full-period |F̂|.

This script answers (1) and (2) at p ∈ {3, 5, 7} for small r, providing data
for the Phase 2 disposition.
"""
from __future__ import annotations

import math
import sys
import csv
from pathlib import Path


def f_p(p: int, r: int, c: int = 1) -> list[complex]:
    """Return f(u) = exp(2*pi*i * c * (1+p)^u / M) for u = 0..p^r - 1
    where M = p^{r+1}."""
    M = p ** (r + 1)
    period = p ** r
    out = []
    pow_1plus_p = 1
    for u in range(period):
        out.append(complex(math.cos(2 * math.pi * c * pow_1plus_p / M),
                            math.sin(2 * math.pi * c * pow_1plus_p / M)))
        pow_1plus_p = (pow_1plus_p * (1 + p)) % M
    return out


def fhat_full_period(f_vals: list[complex], M: int, period: int) -> list[complex]:
    """Compute F_hat(ξ) = sum_{u=0..period-1} f(u) * e^{-2pi i xi u / M}
    for ξ = 0..M-1. (Note: full-period sum, NOT short-window.)"""
    # Direct evaluation — O(period * M) but small.
    Fhat = []
    for xi in range(M):
        s = 0+0j
        for u in range(period):
            phase = -2 * math.pi * xi * u / M
            s += f_vals[u] * complex(math.cos(phase), math.sin(phase))
        Fhat.append(s)
    return Fhat


def K_short_window(p: int, r: int, c: int = 1, m: int = 0) -> complex:
    """K_p(r, c, m) = sum_{u=0..N-1} exp(2pi i (c (1+p)^u - p^2 m u) / M)
    where N = p^{r-1}, M = p^{r+1}."""
    M = p ** (r + 1)
    N = p ** (r - 1)
    s = 0+0j
    pow_1plus_p = 1
    for u in range(N):
        arg = (c * pow_1plus_p - p * p * m * u) / M
        s += complex(math.cos(2 * math.pi * arg), math.sin(2 * math.pi * arg))
        pow_1plus_p = (pow_1plus_p * (1 + p)) % M
    return s


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    OUTPATH = Path("C:/Collatz/qx1_move2_phase2_check.csv")
    rows = []

    # ---- Question 1: generalized R78.3 |F̂_p(ξ)| = p^{(r+3)/2} ----
    print("# Q1: generalized R78.3 — full-period |F_hat_p(xi)| equidistribution")
    print()
    print(f"  {'p':>3} {'r':>2} {'M':>8} {'period':>8} {'|supp|':>8} "
          f"{'predicted |F_hat|':>20} {'max |F_hat| on supp':>22} "
          f"{'min |F_hat| on supp':>22} {'match':>8}")

    for p in [3, 5, 7]:
        for r in [2, 3]:
            M = p ** (r + 1)
            period = p ** r
            f_vals = f_p(p, r)
            Fhat = fhat_full_period(f_vals, M, period)
            mags = [abs(x) for x in Fhat]

            # Support: where mags are non-trivially nonzero. Predicted: {p·a}
            # with a ∈ subset of Z/p^r.
            support_indices = [xi for xi, m in enumerate(mags) if m > 1e-6]
            mags_on_supp = sorted([mags[xi] for xi in support_indices], reverse=True)

            predicted_magnitude = p ** ((r + 3) / 2)

            # Verify predicted: max should equal p^{(r+3)/2} (equidistribution)
            max_mag = max(mags) if mags else 0.0
            # Magnitudes "on supp" — restrict to predicted support {p·a : a ≡ 1 mod p}
            theoretical_supp = []
            for a in range(p ** r):
                if a % p == 1:  # {a ≡ 1 mod p}
                    theoretical_supp.append((p * a) % M)
            mags_theoretical_supp = [mags[xi] for xi in theoretical_supp]
            min_mag_on_theoretical = min(mags_theoretical_supp) if mags_theoretical_supp else 0
            max_mag_on_theoretical = max(mags_theoretical_supp) if mags_theoretical_supp else 0

            match = (
                abs(max_mag_on_theoretical - predicted_magnitude) / predicted_magnitude < 1e-3
                and abs(min_mag_on_theoretical - predicted_magnitude) / predicted_magnitude < 1e-3
            )

            print(
                f"  {p:>3} {r:>2} {M:>8} {period:>8} {len(theoretical_supp):>8} "
                f"{predicted_magnitude:>20.4f} {max_mag_on_theoretical:>22.4f} "
                f"{min_mag_on_theoretical:>22.4f} {'PASS' if match else 'FAIL':>8}"
            )

            rows.append({
                "test": "R78.3_generalized",
                "p": p, "r": r, "M": M, "period": period,
                "supp_size": len(theoretical_supp),
                "predicted_F_hat_magnitude": predicted_magnitude,
                "max_F_hat_on_supp": max_mag_on_theoretical,
                "min_F_hat_on_supp": min_mag_on_theoretical,
                "match": match,
                "K_p_short_window_magnitude": "",
                "predicted_sqrt_N": "",
                "K_over_sqrt_N": "",
            })

    print()
    print("# Q2: K_p(r, c=1, m=0) short-window saturation — empirical |K_p|/sqrt(N)")
    print()
    print(f"  {'p':>3} {'r':>2} {'N':>10} {'sqrt(N)':>12} "
          f"{'|K_p(1,0)|':>14} {'|K|/sqrt(N)':>14}")

    for p in [3, 5, 7]:
        for r in [3, 4, 5]:
            N = p ** (r - 1)
            sqrtN = math.sqrt(N)
            K = K_short_window(p, r, c=1, m=0)
            absK = abs(K)
            ratio = absK / sqrtN if sqrtN > 0 else 0
            print(f"  {p:>3} {r:>2} {N:>10} {sqrtN:>12.4f} {absK:>14.4f} {ratio:>14.4f}")

            rows.append({
                "test": "K_p_short_window",
                "p": p, "r": r, "M": p ** (r + 1), "period": p ** r,
                "supp_size": "",
                "predicted_F_hat_magnitude": "",
                "max_F_hat_on_supp": "",
                "min_F_hat_on_supp": "",
                "match": "",
                "K_p_short_window_magnitude": absK,
                "predicted_sqrt_N": sqrtN,
                "K_over_sqrt_N": ratio,
            })

    print()
    print("# Q2b: K_p at q in {11, 13, 17, 19, 23} small r (A2 supplementary data)")
    print()
    print(f"  {'p':>3} {'r':>2} {'N':>10} {'sqrt(N)':>12} "
          f"{'|K_p(1,0)|':>14} {'|K|/sqrt(N)':>14}")

    for p in [11, 13, 17, 19, 23]:
        # Compute at modest r for these larger primes (period = p^r grows fast)
        for r in [3, 4]:
            try:
                N = p ** (r - 1)
                sqrtN = math.sqrt(N)
                K = K_short_window(p, r, c=1, m=0)
                absK = abs(K)
                ratio = absK / sqrtN if sqrtN > 0 else 0
                print(f"  {p:>3} {r:>2} {N:>10} {sqrtN:>12.4f} {absK:>14.4f} {ratio:>14.4f}")
                rows.append({
                    "test": "K_p_short_window_A2",
                    "p": p, "r": r, "M": p ** (r + 1), "period": p ** r,
                    "supp_size": "",
                    "predicted_F_hat_magnitude": "",
                    "max_F_hat_on_supp": "",
                    "min_F_hat_on_supp": "",
                    "match": "",
                    "K_p_short_window_magnitude": absK,
                    "predicted_sqrt_N": sqrtN,
                    "K_over_sqrt_N": ratio,
                })
            except Exception as e:
                print(f"  {p:>3} {r:>2} ERROR: {e}")

    with open(OUTPATH, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"\n[write] {OUTPATH}")


if __name__ == "__main__":
    main()
