"""R79b side-by-side comparison: |Σ 1̂·ψ_lead| vs |Σ 1̂·G/√q| vs direct K
at r = 4, 6, 8, 10. Verifies the agent's empirical run measured the right
object (Scenario B equivalent) and quantifies the leading-vs-true phase divergence.

Definitions:
  q = 3^{r+1}, N = 3^{r-1}, period = 3^r, supp = {a ≡ 1 mod 3 in Z/3^r}, |supp| = 3^{r-1}
  D(-3a) := Σ_{u=0}^{N-1} e_q(-3au)                         (Dirichlet kernel, length N)
  ψ_lead(a) := e_q(P_a(s*(C_a)))                             (saddle leading order)
  G(a) := Σ_{s=0}^{period-1} e_q(P_a(s))                     (true inner Gauss sum)
  ψ_true(a) := G(a)/√q                                        (|ψ_true| = 1 by Th 78.3 + 78.4)

Compared:
  S_lead(r) := |Σ_{a ∈ supp} D(-3a) · ψ_lead(a)|              (Scenario A)
  S_true(r) := |Σ_{a ∈ supp} D(-3a) · ψ_true(a)|              (Scenario B)
  K_c1m0(r) := |Σ_{u=0}^{N-1} e_q(4^u)|                       (direct, c=1, m=0)
  Plancherel-reconstructed: K_c1m0_reconstr := (3/√q) · S_true (should equal K_c1m0 to numerical precision)
"""
import sys, os, math, cmath, time, csv
from fractions import Fraction
import numpy as np
from numba import njit, prange

sys.stdout.reconfigure(encoding="utf-8")
PI = math.pi

# ----------------------------------------------------------------------
# Truncated 3-adic log helpers (sympy-free, integer-mod-q)
# ----------------------------------------------------------------------

def J_for_p3(m):
    j = 1
    while True:
        x = j + 1; v = 0
        while x % 3 == 0:
            x //= 3; v += 1
        if (j + 1) - v >= m:
            return j
        j += 1


def truncated_3adic_log_mod_q(s, J, q):
    if s == 0:
        return 0
    val = 0
    for j in range(1, J + 1):
        jp = j; v3j = 0
        while jp % 3 == 0:
            jp //= 3; v3j += 1
        sign = 1 if (j - 1) % 2 == 0 else -1
        if 3 ** (j - v3j) >= q:
            continue
        coeff = (sign * 3 ** (j - v3j)) % q
        inv_jp = pow(jp, -1, q)
        s_pow = pow(s, j, q)
        val = (val + coeff * inv_jp % q * s_pow) % q
    return val


def compute_psi_lead_phase(r, a, L_tilde_inv_mod_p_mm1, J, q):
    """Phase φ such that ψ_lead(a) = e_q(φ) per Theorem 78.6 leading order."""
    p_mm1 = 3 ** r
    C_a = (a * L_tilde_inv_mod_p_mm1) % p_mm1
    s_star = ((C_a - 1) // 3) % 3
    L_at_sstar = truncated_3adic_log_mod_q(s_star, J, q)
    return (3 * s_star - C_a * L_at_sstar) % q


def precompute_P_a_coeff_table(r, q, J):
    """Precompute (3s, L(1+3s) mod q) for s = 0..period-1.
    Returns array of length period giving (3s mod q, L mod q).
    Used inside Numba."""
    period = 3 ** r
    Ls = np.zeros(period, dtype=np.int64)
    threes = np.zeros(period, dtype=np.int64)
    for s in range(period):
        Ls[s] = truncated_3adic_log_mod_q(s, J, q)
        threes[s] = (3 * s) % q
    return threes, Ls


# ----------------------------------------------------------------------
# Numba-accelerated computations
# ----------------------------------------------------------------------

@njit(cache=True, parallel=True)
def compute_dirichlet_arr(r, supp_arr):
    """For each a in supp, return D(-3a) = Σ_{u=0}^{N-1} e_q(-3au) as complex.
    supp_arr is a numpy int64 array of a-values."""
    q = 3 ** (r + 1)
    N = 3 ** (r - 1)
    n = len(supp_arr)
    out_re = np.zeros(n, dtype=np.float64)
    out_im = np.zeros(n, dtype=np.float64)
    inv_q = 2.0 * math.pi / q
    for i in prange(n):
        a = supp_arr[i]
        sre = 0.0; sim = 0.0
        # phase = (-3·a·u) mod q
        # increment per u: (-3a) mod q
        delta = (-3 * a) % q
        cur = 0
        for u in range(N):
            angle = inv_q * cur
            sre += math.cos(angle)
            sim += math.sin(angle)
            cur = (cur + delta) % q
        out_re[i] = sre
        out_im[i] = sim
    return out_re, out_im


@njit(cache=True, parallel=True)
def compute_G_arr(r, supp_arr, threes_mod_q, Ls_mod_q, L_tilde_inv, p_mm1):
    """For each a in supp, return G(a) = Σ_s e_q(P_a(s)) where P_a(s) = 3s − C_a·L(1+3s).
    Period = 3^r. Numba parallelizes over a-values."""
    q = 3 ** (r + 1)
    period = 3 ** r
    n = len(supp_arr)
    out_re = np.zeros(n, dtype=np.float64)
    out_im = np.zeros(n, dtype=np.float64)
    inv_q = 2.0 * math.pi / q
    for i in prange(n):
        a = supp_arr[i]
        C_a = (a * L_tilde_inv) % p_mm1
        sre = 0.0; sim = 0.0
        for s_idx in range(period):
            phase = (threes_mod_q[s_idx] - C_a * Ls_mod_q[s_idx]) % q
            angle = inv_q * phase
            sre += math.cos(angle)
            sim += math.sin(angle)
        out_re[i] = sre
        out_im[i] = sim
    return out_re, out_im


@njit(cache=True, parallel=True)
def compute_psi_lead_arr(r, supp_arr, L_tilde_inv, p_mm1, threes_mod_q, Ls_mod_q):
    """For each a in supp, return ψ_lead(a) phase. Vectorized in a."""
    q = 3 ** (r + 1)
    n = len(supp_arr)
    psi_re = np.zeros(n, dtype=np.float64)
    psi_im = np.zeros(n, dtype=np.float64)
    inv_q = 2.0 * math.pi / q
    for i in prange(n):
        a = supp_arr[i]
        C_a = (a * L_tilde_inv) % p_mm1
        s_star = ((C_a - 1) // 3) % 3
        # P_a(s_star) = 3·s_star − C_a · L(1+3·s_star)
        L_at_s = Ls_mod_q[s_star]
        phase = (threes_mod_q[s_star] - C_a * L_at_s) % q
        angle = inv_q * phase
        psi_re[i] = math.cos(angle)
        psi_im[i] = math.sin(angle)
    return psi_re, psi_im


# ----------------------------------------------------------------------
# Main comparison
# ----------------------------------------------------------------------

def main():
    out_csv = r"C:\Collatz\r79b_scenario_comparison.csv"
    fieldnames = ["r", "q", "N", "supp_size", "S_lead", "S_true", "K_c1m0_direct",
                  "K_reconstr_from_S_true", "ratio_lead_true", "ratio_K_to_recon",
                  "elapsed_s"]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # Warmup
        print("[warmup] Numba JIT compile...")
        _ = J_for_p3(3)
        threes_, Ls_ = precompute_P_a_coeff_table(2, 27, J_for_p3(3))
        small_supp = np.array([1, 4, 7], dtype=np.int64)
        _ = compute_dirichlet_arr(2, small_supp)
        _ = compute_G_arr(2, small_supp, threes_, Ls_, 1, 9)
        _ = compute_psi_lead_arr(2, small_supp, 1, 9, threes_, Ls_)
        print("[warmup done]")
        print()

        for r in [4, 6, 8, 10]:
            t0 = time.time()
            q = 3 ** (r + 1)
            N = 3 ** (r - 1)
            period = 3 ** r
            p_mm1 = 3 ** r
            J = J_for_p3(r + 1)
            print(f"## r={r}: q={q}, N={N}, period={period}, J={J}")

            # Build supp
            supp = np.array([a for a in range(p_mm1) if a % 3 == 1], dtype=np.int64)
            print(f"   |supp| = {len(supp)}")

            # Precompute L_tilde_inv
            L4 = truncated_3adic_log_mod_q(1, J, q)
            L_tilde = L4 // 3
            L_tilde_inv = pow(L_tilde, -1, p_mm1)

            # Precompute (3s, L_mod_q) tables
            threes, Ls = precompute_P_a_coeff_table(r, q, J)

            # Compute D(-3a) for each a
            t_d = time.time()
            d_re, d_im = compute_dirichlet_arr(r, supp)
            print(f"   D(-3a) computed in {time.time()-t_d:.2f}s")

            # Compute G(a) for each a
            t_g = time.time()
            g_re, g_im = compute_G_arr(r, supp, threes, Ls, L_tilde_inv, p_mm1)
            print(f"   G(a) computed in {time.time()-t_g:.2f}s")

            # Compute ψ_lead(a) for each a
            t_p = time.time()
            psi_re, psi_im = compute_psi_lead_arr(r, supp, L_tilde_inv, p_mm1, threes, Ls)
            print(f"   ψ_lead(a) computed in {time.time()-t_p:.2f}s")

            # Verify |G(a)| = √q for all a
            G_mag = np.sqrt(g_re ** 2 + g_im ** 2)
            sqrt_q = math.sqrt(q)
            G_mag_ok = np.allclose(G_mag, sqrt_q, atol=1e-6)
            print(f"   |G(a)| = √q = {sqrt_q:.4f} for all a: {G_mag_ok}  (max dev = {np.abs(G_mag - sqrt_q).max():.2e})")

            # Compute Scenario A: S_lead = Σ_a D(-3a) · ψ_lead(a)
            S_lead_re = (d_re * psi_re - d_im * psi_im).sum()
            S_lead_im = (d_re * psi_im + d_im * psi_re).sum()
            S_lead = math.sqrt(S_lead_re ** 2 + S_lead_im ** 2)

            # Compute Scenario B: S_true = Σ_a D(-3a) · (G(a)/√q)
            psi_true_re = g_re / sqrt_q
            psi_true_im = g_im / sqrt_q
            S_true_re = (d_re * psi_true_re - d_im * psi_true_im).sum()
            S_true_im = (d_re * psi_true_im + d_im * psi_true_re).sum()
            S_true = math.sqrt(S_true_re ** 2 + S_true_im ** 2)

            # Direct K(r, c=1, m=0)
            # Reuse the kalafatelis_sum_serial from r79b_compute_S_partial
            K_re = 0.0; K_im = 0.0
            inv_q = 2.0 * PI / q
            x = 1
            for u in range(N):
                phase = x % q
                K_re += math.cos(inv_q * phase)
                K_im += math.sin(inv_q * phase)
                x = (x * 4) % q
            K_direct = math.sqrt(K_re ** 2 + K_im ** 2)

            # Plancherel reconstruction: K = (3·e_q(1)/q) · Σ_a D(-3a)·G(a)
            # |K| = (3/q) · |Σ_a D(-3a)·G(a)| = (3/q) · √q · S_true = (3/√q) · S_true
            K_reconstr = (3.0 / sqrt_q) * S_true

            ratio_lead_true = S_lead / S_true if S_true > 0 else float("nan")
            ratio_K = K_direct / K_reconstr if K_reconstr > 0 else float("nan")

            elapsed = time.time() - t0

            print(f"   S_lead     = {S_lead:>14.4f}  (Scenario A: |Σ D·ψ_lead|)")
            print(f"   S_true     = {S_true:>14.4f}  (Scenario B: |Σ D·G/√q|)")
            print(f"   K_direct   = {K_direct:>14.4f}  (direct Kalafatelis sum)")
            print(f"   K_recon    = {K_reconstr:>14.4f}  (= (3/√q)·S_true; should = K_direct)")
            print(f"   ratio S_lead/S_true   = {ratio_lead_true:.4f}")
            print(f"   ratio K/K_recon (Plancherel check) = {ratio_K:.6f} (should be 1.0)")
            print(f"   total elapsed = {elapsed:.2f}s")
            print()

            row = {
                "r": r, "q": q, "N": N, "supp_size": len(supp),
                "S_lead": f"{S_lead:.6f}", "S_true": f"{S_true:.6f}",
                "K_c1m0_direct": f"{K_direct:.6f}",
                "K_reconstr_from_S_true": f"{K_reconstr:.6f}",
                "ratio_lead_true": f"{ratio_lead_true:.6f}",
                "ratio_K_to_recon": f"{ratio_K:.6f}",
                "elapsed_s": f"{elapsed:.3f}",
            }
            writer.writerow(row)
            f.flush()

    print(f"\n[done] csv = {out_csv}")


if __name__ == "__main__":
    main()
