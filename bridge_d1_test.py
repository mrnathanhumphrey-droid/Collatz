"""
bridge_d1_test.py — Phase 3 empirical test of Candidate D1.

D1 hypothesis: at p = 3, r = n - 1 to match group scale Z/3^n,
the normalised F̂_3(ξ)/M predicts μ̂_n(ξ).

Procedure:
  1. Compute F̂_3(ξ) on Z/3^n via the verified G-FFT method (fhat_verification.py infrastructure).
  2. Compute μ̂_n(ξ) two ways:
     (a) Exact distribution of Syrac(Z/3^n) via Tao's recursion (1.7)/(1.22), feasible up to n ≤ 5.
     (b) Monte-Carlo sampling from Tao (1.26): Σ 3^{i-1} · 2^{-a_{[1,i]}} (mod 3^n) with a_i ~ iid Geom(2).
  3. Compare F̂_3(ξ)/M vs μ̂_n(ξ) at:
     - ξ ∈ supp(F̂_3) (multiples of 3 with ξ/3 ≡ 1 mod 3)
     - ξ ∉ supp(F̂_3) (everything else, including ξ with 3 ∤ ξ — the Tao 1.17 set)

Falsification criterion: if at any ξ with 3 ∤ ξ, |μ̂_n(ξ)| > 1e-3 while F̂_3(ξ)/M = 0,
then F̂_3/M cannot bound μ̂_n at those frequencies — D1 falsified.

Concurrence criterion: if F̂_3(ξ)/M ≈ μ̂_n(ξ) within MC noise at every ξ tested, that is
EMPIRICAL_PATTERN_NO_DERIVATION (A1 safeguard flag); investigate further before claiming bridge.
"""
from __future__ import annotations
import csv
import math
import time
from pathlib import Path

import numpy as np


def f_p_period(p: int, r: int, c: int = 1) -> np.ndarray:
    M = p ** (r + 1)
    period = p ** r
    out = np.empty(period, dtype=np.complex128)
    pow_val = 1
    for u in range(period):
        out[u] = np.exp(2j * np.pi * c * pow_val / M)
        pow_val = (pow_val * (1 + p)) % M
    return out


def F_full_on_ZM(p: int, r: int, c: int = 1) -> np.ndarray:
    """Return F̂_full evaluated on the entire Z/M (length M).

    Implementation: build the periodically-extended f_p on Z/M (length M = p^{r+1})
    by tiling p copies of the period-p^r block, then take length-M FFT.
    Algebraically, F̂_full(p·a) = p · G[a] and F̂_full(ξ) = 0 for ξ ∉ p·Z/M.
    We compute it directly to be honest about off-support behavior.
    """
    M = p ** (r + 1)
    period = p ** r
    f_short = f_p_period(p, r, c)
    f_full = np.tile(f_short, p)  # length M
    return np.fft.fft(f_full)


def syrac_exact_distribution(n: int) -> dict[int, float]:
    """Compute the exact distribution of Syrac(Z/3^n) via Tao's recursion (1.22).

    P(Syrac(Z/3^{n+1}) = x) = Σ_{a≥1: 2^a · x ≡ 1 mod 3} 2^{-a} · P(Syrac(Z/3^n) = (2^a · x - 1)/3 mod 3^n)

    Tao base case: Syrac(Z/3^0) = 0 mod 1 a.s. (so trivially p_0(0) = 1).
    """
    # Base case: at n = 0, point mass at 0.
    p_n: dict[int, float] = {0: 1.0}

    for level in range(n):
        M_next = 3 ** (level + 1)
        p_next: dict[int, float] = {}
        # For each x in Z/M_next, sum over a ≥ 1 with 2^a · x ≡ 1 (mod 3).
        # Truncate a at some practical cutoff (P(a > 40) < 1e-12, etc.).
        a_cutoff = 60
        for x in range(M_next):
            total = 0.0
            for a in range(1, a_cutoff):
                # Need 2^a · x ≡ 1 (mod 3), i.e. (2^a · x - 1) divisible by 3.
                num = (pow(2, a, 3 * M_next) * x - 1) % (3 * M_next)
                if num % 3 != 0:
                    continue
                y = (num // 3) % (3 ** level) if level > 0 else 0
                pr_y = p_n.get(y, 0.0)
                if pr_y > 0:
                    total += (2 ** (-a)) * pr_y
            if total > 0:
                p_next[x] = total
        p_n = p_next

    return p_n


def mu_hat_exact(n: int, xi: int) -> complex:
    """Exact μ̂_n(ξ) from the exact Syrac(Z/3^n) distribution."""
    dist = syrac_exact_distribution(n)
    Mn = 3 ** n
    val = 0.0 + 0.0j
    for x, p in dist.items():
        val += p * np.exp(-2j * np.pi * xi * x / Mn)
    return val


def mu_hat_mc(n: int, xi: int, n_samp: int = 200000, rng=None) -> tuple[complex, float]:
    """Monte-Carlo estimate of μ̂_n(ξ) using Tao (1.26).

    Syrac(Z/3^n) ≡ Σ_{i=1..n} 3^{i-1} · 2^{-a_{[1,i]}}  (mod 3^n)
    with a_i ~ iid Geom(2) (P(a_i = k) = 2^{-k} for k ≥ 1).
    """
    if rng is None:
        rng = np.random.default_rng(0xBEEF)
    Mn = 3 ** n
    # Sample Geom(2): a = floor(-log(U)/log(2)) + 1 with U ~ Uniform(0,1)
    # Equivalently np.random.geometric(p=0.5) gives values 1, 2, ...
    a = rng.geometric(p=0.5, size=(n_samp, n))  # shape (n_samp, n)
    # Cumulative sums: a_[1,j] = a_1 + ... + a_j
    a_cum = np.cumsum(a, axis=1)
    # We need 2^{-a_cum[i]} mod 3^n. Use Python-level modular inverse.
    # 2^{-k} mod 3^n: 2 is invertible mod 3^n, so 2^{-k} ≡ (inv(2))^k mod 3^n.
    inv2 = pow(2, -1, Mn)
    # Vectorise: compute 3^{i-1} · (inv2)^{a_cum[:, i-1]} mod Mn, sum over i.
    # Use Python objects since a_cum can be large; vectorise via pow mod via numpy uint64? not safe.
    # For modest n_samp, do a Python loop per sample. For larger, use numpy with int64 + bitmask carefully.
    # Cleanest: per sample, accumulate the sum mod Mn.
    syrac_vals = np.zeros(n_samp, dtype=np.int64)
    for s in range(n_samp):
        total = 0
        three_pow = 1
        for i in range(n):
            k = int(a_cum[s, i])
            term = (three_pow * pow(inv2, k, Mn)) % Mn
            total = (total + term) % Mn
            three_pow = (three_pow * 3) % Mn
        syrac_vals[s] = total
    # μ̂_n(ξ) = (1/n_samp) Σ exp(-2πi ξ S/Mn)
    phase = np.exp(-2j * np.pi * xi * syrac_vals.astype(np.float64) / Mn)
    val = phase.mean()
    # Standard error magnitude (rough): 1/sqrt(n_samp)
    se = 1.0 / math.sqrt(n_samp)
    return complex(val), se


def main() -> None:
    import sys
    sys.stdout.reconfigure(encoding="utf-8")

    rows = []

    # Cells: p = 3, n ∈ {2, 3, 4, 5}, r = n - 1 (so M = 3^n matches Z/3^n).
    print("=" * 88)
    print("BRIDGE Phase 3 — Candidate D1: F̂_3(ξ)/M vs μ̂_n(ξ) on Z/3^n")
    print("=" * 88)

    for n in [2, 3, 4, 5]:
        r = n - 1
        p = 3
        M = p ** (r + 1)
        assert M == 3 ** n, (M, n)

        print(f"\n--- n = {n} (r = {r}, M = 3^n = {M}) ---")

        # F̂_full on Z/M.
        F_full = F_full_on_ZM(p=p, r=r, c=1)

        # Exact μ̂_n distribution (computationally feasible for n ≤ 5).
        if n <= 4:
            dist = syrac_exact_distribution(n)
        else:
            dist = None  # MC only

        # Compare at each ξ ∈ Z/M.
        print(f"  {'xi':>4} {'in_F_supp':>10} {'F_full(xi)/M':>22} "
              f"{'mu_hat_exact':>22} {'mu_hat_mc':>22} {'mc_se':>10}")

        # Predicted F̂ support: ξ = p·a with a ∈ Z/p^r and a ≡ 1 (mod p).
        supp_F = {3 * a for a in range(0, 3 ** r) if a % 3 == 1}

        for xi in range(M):
            F_val = F_full[xi] / M
            in_supp = xi in supp_F
            # μ̂ exact (small n)
            if dist is not None:
                mu_exact = 0.0 + 0.0j
                for x, pr in dist.items():
                    mu_exact += pr * np.exp(-2j * np.pi * xi * x / M)
            else:
                mu_exact = None
            # μ̂ MC (always)
            n_samp = 50000 if n <= 3 else (30000 if n == 4 else 20000)
            mu_mc, mc_se = mu_hat_mc(n=n, xi=xi, n_samp=n_samp)

            # Only print first few xi and ones with |mu| or |F| above threshold for brevity.
            print_this = (xi < 6) or (in_supp) or (abs(mu_mc) > 0.05)

            if print_this:
                f_str = f"{F_val.real:+.4f}{F_val.imag:+.4f}j"
                me_str = f"{mu_exact.real:+.4f}{mu_exact.imag:+.4f}j" if mu_exact is not None else "—"
                mc_str = f"{mu_mc.real:+.4f}{mu_mc.imag:+.4f}j"
                print(f"  {xi:>4} {'YES' if in_supp else 'no':>10} "
                      f"{f_str:>22} {me_str:>22} {mc_str:>22} {mc_se:>10.4f}")

            rows.append({
                "n": n, "M": M, "xi": xi,
                "xi_mod_3": xi % 3,
                "in_F_supp": int(in_supp),
                "F_real": F_val.real, "F_imag": F_val.imag, "F_mag": abs(F_val),
                "mu_exact_real": mu_exact.real if mu_exact is not None else "",
                "mu_exact_imag": mu_exact.imag if mu_exact is not None else "",
                "mu_exact_mag": abs(mu_exact) if mu_exact is not None else "",
                "mu_mc_real": mu_mc.real, "mu_mc_imag": mu_mc.imag, "mu_mc_mag": abs(mu_mc),
                "mc_se": mc_se,
            })

    # ============================================================
    # Comparison summary
    # ============================================================
    print()
    print("=" * 88)
    print("SUMMARY: D1 falsification check")
    print("=" * 88)

    falsified_pairs = []
    for row in rows:
        if row["xi_mod_3"] != 0:
            # ξ with 3 ∤ ξ: F̂_3 must vanish; μ̂ may be nonzero per Tao 1.17.
            if abs(row["F_mag"]) > 1e-8 and row["mu_mc_mag"] > 3 * row["mc_se"]:
                # Both nonzero — interesting
                pass
            elif abs(row["F_mag"]) < 1e-8 and row["mu_mc_mag"] > 0.01:
                falsified_pairs.append((row["n"], row["xi"], row["F_mag"], row["mu_mc_mag"]))

    print(f"  Number of (n, ξ) pairs where 3 ∤ ξ, F̂ ≈ 0, but |μ̂_mc| > 0.01:")
    print(f"  {len(falsified_pairs)} pairs (each is a D1 falsification).")
    if falsified_pairs:
        print("  First 10:")
        for n_, xi_, fm, mm in falsified_pairs[:10]:
            print(f"    n={n_}, ξ={xi_}, |F̂|={fm:.2e}, |μ̂_mc|={mm:.4f}")

    # CSV out
    csv_path = Path("C:/Collatz/bridge_d1_results.csv")
    keys = sorted({k for r_ in rows for k in r_.keys()})
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r_ in rows:
            w.writerow(r_)
    print(f"\n[write] {csv_path}")


if __name__ == "__main__":
    main()
