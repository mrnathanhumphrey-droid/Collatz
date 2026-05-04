"""C2: Bourgain-Konyagin sum-product on multiplicative subgroup ⟨4⟩.

Empirical test: compute moments M_{2k}(r) = Σ_{c ∈ (Z/q)*} |K(r, c, 0)|^{2k}
for k = 1, 2, 3, 4 at r ∈ {6, 8, 10, 12}, where K(r, c, 0) = Σ_{u=0}^{N-1} e_q(c·4^u).

Theoretical references:
  M_2 = q · N (Plancherel, exact)
  M_4 = q · E_×({4^u}) where E_× is the multiplicative energy of {4^u : u=0..N-1}.
        Since 4 has order 3^r in (Z/q)^*, {4^u} ↔ {0..N-1} ⊂ Z/3^r via the iso, and
        E_×({4^u}) = E_+([0, N-1] in Z/3^r) = E_+([0, N-1] in Z) (since N < 3^r/3,
        no wrap). E_+([0, N-1]) = (2/3)N³ + O(N²).
  So M_4 ≈ (2/3) q N³ = 6 N^4 (using q = 9N).

  BK saving prediction (if applies): M_4 ≤ q · |H|^{4-2δ_BK} for some δ_BK > 0.
  For δ_BK = 0: M_4 ≤ q · |H|^4 (trivial bound on subgroup sum).
  Our empirical M_4 will tell us where the saving is.

Comparison:
  - Random model (uncorrelated phases, |K|~√N): M_4 = #c · 2N² = 6N · 2N² = 12N³.
    [factor 2 from 4!/(2!·2!) ways to pair (u1,u2)=(u3,u4) with u1=u3, u2=u4 OR u1=u4, u2=u3.]
    So random predicts M_4 ≈ 12N³.
  - Full multiplicative energy (saturated): M_4 = q · 2N³/3 = 6 N^4.
  - Ratio M_4_saturated / M_4_random = 6N^4 / 12N³ = N/2.
    So saturation is N/2 times higher than random — a HUGE difference.

Empirical M_4 close to 6N^4 ⟹ saturated, no BK saving.
Empirical M_4 close to 12N³ ⟹ random-like, BK could give saving.
"""
import sys, os, math, time, csv
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")

PI = math.pi


@njit(cache=True, parallel=True)
def K_for_all_c_units(r: int) -> np.ndarray:
    """Compute |K(r, c, 0)| for all c ∈ (Z/q)^*. Returns array of length 2·3^r."""
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    # Build list of c-units: c ∈ [1, q-1] with gcd(c, 3) = 1
    n_units = 2 * 3 ** r  # |(Z/q)^*|
    out = np.zeros(n_units, dtype=np.float64)
    inv_q = 2.0 * PI / q
    for idx in prange(n_units):
        # Map idx to c-unit:
        # c-units are {c ∈ [1, q-1] : c % 3 != 0}. There are 2·3^r of them.
        # Index idx = 0 → c=1, 1 → c=2, 2 → c=4, 3 → c=5, 4 → c=7, 5 → c=8, ...
        # i.e., for idx, take c = (idx // 2) * 3 + (idx % 2 + 1)
        c = (idx // 2) * 3 + (idx % 2) + 1
        sre = 0.0; sim = 0.0
        x = 1
        for u in range(N):
            phase_int = (c * x) % q
            angle = inv_q * phase_int
            sre += math.cos(angle)
            sim += math.sin(angle)
            x = (x * 4) % q
        out[idx] = math.sqrt(sre * sre + sim * sim)
    return out


def true_4th_moment_via_energy(r: int) -> float:
    """M_4 = q · E_+([0, N-1]) in Z, exact closed form.
    E_+([0, N-1]) = Σ_{t=0}^{2N-2} r²(t)
                  = 2·(1² + 2² + ... + (N-1)²) + N²
                  = N(N-1)(2N-1)/3 + N²
    """
    N = 3 ** (r - 1)
    q = 3 ** (r + 1)
    E_plus = N * (N - 1) * (2 * N - 1) // 3 + N * N
    return q * E_plus


def main():
    out_csv = r"C:\Collatz\bk_moments_data.csv"
    fieldnames = ["r", "q", "N", "n_units",
                  "M_2_emp", "M_2_predicted", "M_2_match",
                  "M_4_emp", "M_4_predicted_saturated", "M_4_random_model",
                  "M_4_emp_over_saturated", "M_4_emp_over_random",
                  "M_6_emp", "M_8_emp",
                  "max_abs_K", "median_abs_K", "p99_abs_K",
                  "elapsed_s"]

    print("="*88)
    print("# C2: Bourgain-Konyagin moment analysis on ⟨4⟩ ⊂ (Z/3^{r+1})^*")
    print("# M_{2k}(r) := Σ_{c ∈ (Z/q)^*} |K(r, c, 0)|^{2k}")
    print("="*88)
    print()

    # Warmup numba
    print("[warmup]")
    _ = K_for_all_c_units(3)
    print("[warmup done]")
    print()

    rows = []
    for r in [6, 8, 10, 12]:
        t0 = time.time()
        q = 3 ** (r + 1)
        N = 3 ** (r - 1)
        n_units = 2 * 3 ** r

        print(f"## r = {r}: q = {q}, N = {N}, |(Z/q)^*| = {n_units}")

        # Compute |K(r, c, 0)| for all c-units
        abs_K = K_for_all_c_units(r)
        elapsed_compute = time.time() - t0
        print(f"   computed |K| for {n_units} c-units in {elapsed_compute:.2f}s")

        # Moments
        K2 = abs_K ** 2
        K4 = abs_K ** 4
        K6 = abs_K ** 6
        K8 = abs_K ** 8
        M_2 = float(K2.sum())
        M_4 = float(K4.sum())
        M_6 = float(K6.sum())
        M_8 = float(K8.sum())

        # Predictions
        M_2_pred = q * N  # Plancherel exact (over all c ∈ Z/q, then restricted to units)
        # Actually Σ_{c ∈ Z/q} |K|² = q · N. The c=0 term contributes |K(c=0)|² = |Σ_u 1|² = N².
        # And Σ_{c ≡ 0 mod 3, c ≠ 0} contributes the rest. For (Z/q)^* (gcd(c,3)=1): subtract.
        # Easier: just compute predicted = q·N, expect empirical < this.
        M_4_saturated = true_4th_moment_via_energy(r)  # over all c ∈ Z/q
        M_4_random = 12 * N ** 3 * (n_units / (q))  # rescale to units only
        # Actually random model: for "typical" c-unit, |K|² ≈ N. Σ over n_units = n_units · N. M_2 = n_units·N
        # But we computed M_2 = sum over units of |K|². So M_2 ≈ n_units · N for random model.
        # Saturated M_2 = q · N (all c ∈ Z/q including c=0). For c-units only: subtract c=0 (giving N²) and c divisible by 3.
        # We'll just compare empirical to absolute counts.

        # Statistics on |K|
        max_abs = float(abs_K.max())
        median_abs = float(np.median(abs_K))
        p99_abs = float(np.percentile(abs_K, 99))

        # Ratios
        M_4_over_saturated_total = M_4 / M_4_saturated  # if computing over all Z/q gives saturated, units gives ~2/3
        M_4_over_random_units = M_4 / (n_units * 2 * N**2)  # random predicts 2N² per unit
        M_2_match = M_2 / (n_units * N)  # should be ≈1 for square-root cancellation

        print(f"   M_2 emp = {M_2:.4e},  per-unit avg = {M_2/n_units:.4f},  "
              f"prediction n_units·N = {n_units*N:.4e},  ratio = {M_2_match:.4f}")
        print(f"   M_4 emp = {M_4:.4e}")
        print(f"     vs saturated (q·E_+) = {M_4_saturated:.4e},  ratio = {M_4/M_4_saturated:.4f}")
        print(f"     vs random model (n_units · 2N²) = {n_units*2*N**2:.4e},  ratio = {M_4_over_random_units:.4f}")
        print(f"   M_6 emp = {M_6:.4e}")
        print(f"   M_8 emp = {M_8:.4e}")
        print(f"   |K| stats: max = {max_abs:.4f}, p99 = {p99_abs:.4f}, median = {median_abs:.4f}, √N = {math.sqrt(N):.4f}")
        print(f"   max/√N = {max_abs/math.sqrt(N):.4f}, median/√N = {median_abs/math.sqrt(N):.4f}")
        print()

        rows.append({
            "r": r, "q": q, "N": N, "n_units": n_units,
            "M_2_emp": f"{M_2:.6e}",
            "M_2_predicted": f"{n_units*N:.6e}",
            "M_2_match": f"{M_2_match:.6f}",
            "M_4_emp": f"{M_4:.6e}",
            "M_4_predicted_saturated": f"{M_4_saturated:.6e}",
            "M_4_random_model": f"{n_units * 2 * N**2:.6e}",
            "M_4_emp_over_saturated": f"{M_4/M_4_saturated:.6f}",
            "M_4_emp_over_random": f"{M_4_over_random_units:.6f}",
            "M_6_emp": f"{M_6:.6e}",
            "M_8_emp": f"{M_8:.6e}",
            "max_abs_K": f"{max_abs:.4f}",
            "median_abs_K": f"{median_abs:.4f}",
            "p99_abs_K": f"{p99_abs:.4f}",
            "elapsed_s": f"{elapsed_compute:.3f}",
        })

    # Save CSV
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[csv: {out_csv}]")
    print()

    # Power-law fits
    print("="*88)
    print("# Moment scaling: log(M_{2k}) = a + b · log(N)")
    print("# Saturated prediction: M_4 ∝ N^4, M_6 ∝ N^?, M_8 ∝ N^?")
    print("# Random prediction:    M_4 ∝ N^3, M_6 ∝ N^4, M_8 ∝ N^5")
    print("="*88)
    print()
    rs = [6, 8, 10, 12]
    Ns = [3**(r-1) for r in rs]
    log_Ns = [math.log(N) for N in Ns]
    M2s = [float(r["M_2_emp"]) for r in rows]
    M4s = [float(r["M_4_emp"]) for r in rows]
    M6s = [float(r["M_6_emp"]) for r in rows]
    M8s = [float(r["M_8_emp"]) for r in rows]

    def fit_slope(ys):
        log_ys = [math.log(y) for y in ys]
        x = np.array(log_Ns); y = np.array(log_ys)
        x_m = x.mean(); y_m = y.mean()
        ssxy = ((x - x_m) * (y - y_m)).sum()
        ssxx = ((x - x_m) ** 2).sum()
        b = ssxy / ssxx if ssxx > 0 else 0
        a = y_m - b * x_m
        y_pred = a + b * x
        ss_res = ((y - y_pred)**2).sum()
        ss_tot = ((y - y_m)**2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
        return b, a, r2

    b2, a2, r2_2 = fit_slope(M2s)
    b4, a4, r2_4 = fit_slope(M4s)
    b6, a6, r2_6 = fit_slope(M6s)
    b8, a8, r2_8 = fit_slope(M8s)

    print(f"  M_2:  β = {b2:.4f},  R² = {r2_2:.4f}  (Plancherel predicts β = 2 since q·N = 9N²)")
    print(f"  M_4:  β = {b4:.4f},  R² = {r2_4:.4f}  (saturated β=4, random β=3)")
    print(f"  M_6:  β = {b6:.4f},  R² = {r2_6:.4f}  (saturated β=5, random β=4)")
    print(f"  M_8:  β = {b8:.4f},  R² = {r2_8:.4f}  (saturated β=6, random β=5)")
    print()
    print("# Verdict:")
    if abs(b4 - 4) < 0.2:
        print(f"  M_4 slope {b4:.3f} ≈ 4 ⟹ MULTIPLICATIVE ENERGY SATURATED.")
        print(f"     ⟹ {{4^u : u=0..N-1}} has FULL multiplicative energy of [0,N-1]+[0,N-1].")
        print(f"     ⟹ Bourgain-Konyagin sum-product gives NO saving for this set.")
        print(f"     ⟹ C2 closure path is CLOSED.")
    elif abs(b4 - 3) < 0.2:
        print(f"  M_4 slope {b4:.3f} ≈ 3 ⟹ random-like multiplicative structure.")
        print(f"     ⟹ BK could give a polynomial saving.")
        print(f"     ⟹ C2 closure path requires further analytical follow-up.")
    else:
        print(f"  M_4 slope {b4:.3f} intermediate.")
        print(f"     ⟹ Some multiplicative structure but not full saturation.")
        print(f"     ⟹ Detailed analysis needed.")


if __name__ == "__main__":
    main()
