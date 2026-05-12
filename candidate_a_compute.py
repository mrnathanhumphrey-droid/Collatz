"""
candidate_a_compute.py — Reading A scoping Candidate A minimum-viable test.

Pre-registered patterns (READING_A_SCOPING_MIN_VIABLE_TEST.md):
  A1: dominant-k rate-0.5 decay
  A2: phase-cancellation signed-sum rate-0.5
  F1: all c_{n,k} ≈ 0
  F2: c_{n,k} dominated by k = n-1
  F3: rate ≠ 1/2

Goal: compute c_{n,k} := <φ_n, lift_n(R_k)> over Q (exact rationals via
fractions.Fraction) for all (n, k) with 0 ≤ k < n ≤ 6.

φ_n construction:
  S_n = M_n(1) = Σ_{ξ: 3∤ξ in Z/3^n} |μ̂_n(ξ)|^2
      = Σ_r Σ_s π_n(r) π_n(s) K_n(r-s)
  where K_n(d) := Σ_{ξ: 3∤ξ in Z/3^n} e^{-2πi d ξ / 3^n}.

  K_n(d) simplifies to an exact integer-valued function:
      K_n(d) = 3^n · 1[d ≡ 0 mod 3^n] − 3^{n−1} · 1[d ≡ 0 mod 3^{n−1}].
  Specifically, on the difference set d ∈ {-(3^n-1), ..., 3^n-1}:
    K_n(0)             = 3^n − 3^{n−1} = 2 · 3^{n−1}
    K_n(±3^{n−1})      = − 3^{n−1}
    K_n(±2·3^{n−1})    = − 3^{n−1}
    K_n(d)             = 0           otherwise

  Define φ_n(r) := Σ_s π_n(s) K_n(r−s). Then ⟨φ_n, π_n⟩ = S_n.

  For the test we use the moment functional definition:
    ε_n = S_n − 7/15 = ⟨φ_n, π_n⟩ − 7/15
        = ⟨φ_n, π_n − π_∞^{(n)}⟩  +  (⟨φ_n, π_∞^{(n)}⟩ − 7/15)
  where π_∞^{(n)} is uniform on V_n (= 1/N_n on coprime states).

  We track the residual term ⟨φ_n, π_∞^{(n)}⟩ explicitly; if it equals 7/15
  exactly (or has a clean structural form), then ⟨φ_n, π_n − π_∞^{(n)}⟩ = ε_n
  exactly, and we can use Σ_k c_{n,k} = ε_n as the decomposition sanity check.
  Otherwise, c_{n,k} still measures ⟨φ_n, lift_n(R_k)⟩, the W_k contributions
  to the bilinear-pair moment, and the rate-1/2 question is whether |c_{n,k}|
  or |S_n^signed| := |Σ_k c_{n,k}| scales like (1/2)^n.

Lift construction:
  R_k ∈ V_{k+1} from R77.5 (R_k := π_{k+1} − T(π_k)).
  For k = 0: R_0 := π_1 − π_∞^{(1)} (the "level 0 residual" — π_1 minus
  uniform on V_1).
  lift_n(R_k) ∈ V_n: apply T^{n−1−k} to R_k, where T is the uniform lift.
  Algebraically: lift_n(R_k)(r') = R_k(r' mod 3^{k+1}) / 3^{n−k−1}
  for r' coprime in Z/3^n.

Outputs:
  C:\\Collatz\\candidate_a_c_nk_table.csv
  C:\\Collatz\\candidate_a_phi_n_verify.csv
  C:\\Collatz\\candidate_a_lift_orthogonality.csv
"""
import sys
import os
import time
import csv
from fractions import Fraction

sys.path.insert(0, r"C:\Collatz")
from result_77_5_compute_R_k import (
    pi_dict, lift_pi, sum_entries, squared_l2_norm
)

sys.stdout.reconfigure(encoding="utf-8")
OUT_C_NK = r"C:\Collatz\candidate_a_c_nk_table.csv"
OUT_PHI_VERIFY = r"C:\Collatz\candidate_a_phi_n_verify.csv"
OUT_LIFT_ORTHO = r"C:\Collatz\candidate_a_lift_orthogonality.csv"


def K_n_value(d, n):
    """K_n(d) = #{ξ in Z/3^n with 3∤ξ}*(closed form).
    K_n(d) = 3^n * 1[d % 3^n == 0]  -  3^{n-1} * 1[d % 3^{n-1} == 0].
    d is any integer (positive or negative).
    """
    N = 3 ** n
    Nm1 = 3 ** (n - 1)
    d_mod_N = d % N
    val = 0
    if d_mod_N == 0:
        val += N
    d_mod_Nm1 = d % Nm1
    if d_mod_Nm1 == 0:
        val -= Nm1
    return val


def compute_phi_n(pi_n, coprime_n, n):
    """φ_n(r) := Σ_s π_n(s) K_n(r-s), as dict r -> Fraction.

    pi_n: dict r -> Fraction on coprime states in Z/3^n.
    """
    N = 3 ** n
    Nm1 = 3 ** (n - 1)
    # K_n is supported on d ≡ 0 mod 3^{n-1}, i.e., d ∈ {0, ±3^{n-1}, ±2·3^{n-1}}
    # (interpreted as 5 cosets mod 3^n, with d=0 having coefficient 2·3^{n-1}
    # and the other 4 having coefficient -3^{n-1}).
    # K_n(0 mod 3^n)        = 2 * 3^{n-1}
    # K_n(3^{n-1} mod 3^n)  = -3^{n-1}     (this is d = 3^{n-1})
    # K_n(2*3^{n-1} mod 3^n)= -3^{n-1}     (this is d = -3^{n-1})
    # So K_n(d) is a function of (d mod 3^n).
    phi = {r: Fraction(0) for r in coprime_n}
    for r in coprime_n:
        acc = Fraction(0)
        for s, p_s in pi_n.items():
            d = (r - s) % N
            # d mod 3^{n-1} must be 0 for K_n to be nonzero
            if d % Nm1 != 0:
                continue
            if d == 0:
                acc += p_s * Fraction(2 * Nm1)
            else:
                # d == 3^{n-1} or d == 2*3^{n-1}
                acc -= p_s * Fraction(Nm1)
        phi[r] = acc
    return phi


def inner_product(u, v):
    """Σ_r u[r] * v[r] over the index set of u (assumed shared)."""
    s = Fraction(0)
    for r in u:
        if r in v:
            s += u[r] * v[r]
    return s


def lift_n_of_R_k(R_k, k, n):
    """Lift R_k ∈ V_{k+1} to V_n via T^{n-k-1}.
    Algebraically: lift_n(R_k)(r') = R_k(r' mod 3^{k+1}) / 3^{n-k-1}
    for r' coprime in Z/3^n. (k+1 ≤ n required.)
    """
    assert k + 1 <= n
    if n == k + 1:
        return dict(R_k)
    N_n = 3 ** n
    N_kp1 = 3 ** (k + 1)
    scale = Fraction(1, 3 ** (n - k - 1))
    out = {}
    for rp in range(N_n):
        if rp % 3 == 0:
            continue
        r_lower = rp % N_kp1
        # r_lower is automatically coprime in Z/3^{k+1} since rp coprime
        out[rp] = R_k[r_lower] * scale
    return out


def R_0_residual(pi_1, coprime_1):
    """R_0 := π_1 − π_∞^{(1)} = π_1 − uniform on V_1.
    V_1 has N_1 = 2 coprime states {1, 2}.
    π_∞^{(1)} is uniform = 1/2 on each coprime state.
    """
    N1 = len(coprime_1)
    uniform = Fraction(1, N1)
    R0 = {}
    for r in coprime_1:
        R0[r] = pi_1[r] - uniform
    return R0


def main():
    print("# Candidate A — minimum-viable test")
    print("# Computing c_{n,k} := <φ_n, lift_n(R_k)> over Q for 0 ≤ k < n ≤ 6")
    print()

    # Stage 0: π_k for k = 1..6
    pi_at = {}
    coprime_at = {}
    t0 = time.time()
    for k in [1, 2, 3, 4, 5, 6]:
        t1 = time.time()
        pi, cop = pi_dict(k)
        pi_at[k] = pi
        coprime_at[k] = cop
        print(f"  π_{k}: {len(cop)} coprime states, solved in {time.time()-t1:.2f}s")
    print(f"  total Stage 0: {time.time()-t0:.2f}s")
    print()

    # Stage 1: R_k for k = 0..5
    # R_0 := π_1 - uniform on V_1
    # R_k (k≥1) := π_{k+1} - T(π_k) lives on V_{k+1}
    R = {}
    R[0] = R_0_residual(pi_at[1], coprime_at[1])
    print(f"  R_0: lives on V_1 ({len(R[0])} entries), sum = {sum_entries(R[0])}, "
          f"‖R_0‖² = {float(squared_l2_norm(R[0])):.6e}")
    for k in [1, 2, 3, 4, 5]:
        T_pi_k = lift_pi(pi_at[k], k)
        R_k = {}
        for rp in pi_at[k + 1]:
            R_k[rp] = pi_at[k + 1][rp] - T_pi_k[rp]
        s = sum_entries(R_k)
        assert s == 0, f"R_{k} sum != 0: {s}"
        R[k] = R_k
        print(f"  R_{k}: lives on V_{k+1} ({len(R_k)} entries), "
              f"‖R_{k}‖² = {float(squared_l2_norm(R_k)):.6e}")
    print()

    # Stage 2: φ_n for n = 1..6
    phi_at = {}
    S_at = {}
    inner_phi_uniform = {}
    eps_check = {}
    for n in [1, 2, 3, 4, 5, 6]:
        t1 = time.time()
        phi_n = compute_phi_n(pi_at[n], coprime_at[n], n)
        phi_at[n] = phi_n
        # Verify <φ_n, π_n> = S_n
        S_n = inner_product(phi_n, pi_at[n])
        S_at[n] = S_n
        # Compute <φ_n, π_∞^{(n)}> = (1/N_n) * Σ_r φ_n(r)
        N_n = len(coprime_at[n])
        sum_phi = Fraction(0)
        for r in phi_n:
            sum_phi += phi_n[r]
        inner_phi_unif = sum_phi / N_n
        inner_phi_uniform[n] = inner_phi_unif
        eps_n = S_n - Fraction(7, 15)
        eps_check[n] = eps_n
        print(f"  φ_{n}: built in {time.time()-t1:.2f}s")
        print(f"     S_{n} = <φ_{n}, π_{n}> = {S_n.numerator}/{S_n.denominator} ≈ {float(S_n):.10f}")
        print(f"     <φ_{n}, π_∞^{(n)}> = {inner_phi_unif.numerator}/{inner_phi_unif.denominator} ≈ {float(inner_phi_unif):.10f}")
        print(f"     ε_{n} = S_{n} - 7/15 = {eps_n.numerator}/{eps_n.denominator} ≈ {float(eps_n):.6e}")
    print()

    # Cross-check S_n against known R76 values
    S_known = {
        1: Fraction(2, 3),
        2: Fraction(10, 21),
        3: Fraction(31370, 67963),
    }
    print("# Verify S_n against R76 known values:")
    for n in [1, 2, 3]:
        ok = S_at[n] == S_known[n]
        print(f"  n={n}: S_n = {S_at[n]}, R76 = {S_known[n]}, match = {ok}")
    print()

    # Compute c_{n,k} := <φ_n, lift_n(R_k)>
    c_nk = {}
    sum_c_n = {}
    for n in [1, 2, 3, 4, 5, 6]:
        sum_c = Fraction(0)
        for k in range(n):
            t1 = time.time()
            lift_R_k = lift_n_of_R_k(R[k], k, n)
            cval = inner_product(phi_at[n], lift_R_k)
            c_nk[(n, k)] = cval
            sum_c += cval
            print(f"  c_{{{n},{k}}} = <φ_{n}, lift_{n}(R_{k})> = "
                  f"{cval.numerator}/{cval.denominator} ≈ {float(cval):+.6e}  "
                  f"({time.time()-t1:.2f}s)")
        sum_c_n[n] = sum_c
        eps_check_n = eps_check[n]
        # Decomposition check: Σ_k c_{n,k} should equal <φ_n, π_n - π_∞^{(n)}> = S_n - inner_phi_uniform[n]
        target = S_at[n] - inner_phi_uniform[n]
        match = (sum_c == target)
        print(f"  -> Σ_k c_{{{n},k}} = {float(sum_c):+.6e}, "
              f"<φ_n, π_n - π_∞^{(n)}> = {float(target):+.6e}, match = {match}")
        print(f"     ε_n = {float(eps_check_n):+.6e}")
        print()

    # Lift orthogonality check
    print("# Lift orthogonality: <lift_n(R_k1), lift_n(R_k2)> for k1 != k2")
    lift_ortho_rows = []
    for n in [3, 4, 5, 6]:
        for k1 in range(n):
            for k2 in range(k1 + 1, n):
                l1 = lift_n_of_R_k(R[k1], k1, n)
                l2 = lift_n_of_R_k(R[k2], k2, n)
                ip = inner_product(l1, l2)
                lift_ortho_rows.append((n, k1, k2, ip))
                if ip != 0:
                    print(f"  WARN n={n} k1={k1} k2={k2}: <lift,lift> = "
                          f"{ip.numerator}/{ip.denominator} != 0")
    if all(row[3] == 0 for row in lift_ortho_rows):
        print(f"  ✓ All cross-lift inner products exactly 0 ({len(lift_ortho_rows)} pairs checked)")
    print()

    # Lift norm check: ‖lift_n(R_k)‖² = ‖R_k‖² / 3^{n-k-1}
    print("# Lift norm scaling: ‖lift_n(R_k)‖² = ‖R_k‖² / 3^{n-k-1}")
    for n in [3, 4, 5, 6]:
        for k in range(n):
            l = lift_n_of_R_k(R[k], k, n)
            ln = squared_l2_norm(l)
            rn = squared_l2_norm(R[k])
            scale = Fraction(1, 3 ** (n - k - 1)) if n > k + 1 else Fraction(1)
            target = rn * scale
            ok = ln == target
            if not ok:
                print(f"  WARN n={n} k={k}: ‖lift‖²={ln}, target={target}, match={ok}")
    print(f"  (per-lift norm scaling validated)")
    print()

    # Write CSV outputs
    with open(OUT_C_NK, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["n", "k", "c_nk_num", "c_nk_den", "c_nk_float",
                    "abs_c_nk_float", "sign"])
        for n in [1, 2, 3, 4, 5, 6]:
            for k in range(n):
                cval = c_nk[(n, k)]
                w.writerow([
                    n, k,
                    str(cval.numerator),
                    str(cval.denominator),
                    f"{float(cval):.15e}",
                    f"{abs(float(cval)):.15e}",
                    "+" if cval > 0 else ("-" if cval < 0 else "0"),
                ])
    print(f"[save] {OUT_C_NK}")

    with open(OUT_PHI_VERIFY, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["n", "S_n_num", "S_n_den", "S_n_float",
                    "phi_unif_num", "phi_unif_den", "phi_unif_float",
                    "eps_n_num", "eps_n_den", "eps_n_float",
                    "sum_c_n_num", "sum_c_n_den", "sum_c_n_float",
                    "decomp_target_num", "decomp_target_den", "decomp_match"])
        for n in [1, 2, 3, 4, 5, 6]:
            S_n = S_at[n]
            pu = inner_phi_uniform[n]
            en = eps_check[n]
            sc = sum_c_n[n]
            target = S_n - pu
            w.writerow([
                n,
                str(S_n.numerator), str(S_n.denominator), f"{float(S_n):.15e}",
                str(pu.numerator), str(pu.denominator), f"{float(pu):.15e}",
                str(en.numerator), str(en.denominator), f"{float(en):.15e}",
                str(sc.numerator), str(sc.denominator), f"{float(sc):.15e}",
                str(target.numerator), str(target.denominator),
                "TRUE" if sc == target else "FALSE",
            ])
    print(f"[save] {OUT_PHI_VERIFY}")

    with open(OUT_LIFT_ORTHO, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["n", "k1", "k2", "inner_product_num", "inner_product_den",
                    "inner_product_float"])
        for (n, k1, k2, ip) in lift_ortho_rows:
            w.writerow([
                n, k1, k2,
                str(ip.numerator), str(ip.denominator),
                f"{float(ip):.15e}",
            ])
    print(f"[save] {OUT_LIFT_ORTHO}")

    # Pattern matching
    print()
    print("=" * 72)
    print("# Pattern matching across n=2..6")
    print("=" * 72)

    # A1: dominant-k rate
    print()
    print("## A1 — dominant-k rate-0.5 decay")
    k_star = {}
    abs_c_kstar = {}
    for n in [1, 2, 3, 4, 5, 6]:
        ks = -1
        best = Fraction(0)
        for k in range(n):
            ac = abs(c_nk[(n, k)])
            if ac > best:
                best = ac
                ks = k
        k_star[n] = ks
        abs_c_kstar[n] = best
        total_abs = sum(abs(c_nk[(n, k)]) for k in range(n))
        share = float(best / total_abs) if total_abs != 0 else 0
        print(f"  n={n}: k*={ks}, |c_{{n,k*}}| = {float(best):.6e}, "
              f"share of Σ|c| = {share:.4f}")
    print()
    print("## A1 ratios |c_{n,k*(n)}| / |c_{n-1,k*(n-1)}|:")
    for n in [2, 3, 4, 5, 6]:
        if abs_c_kstar[n - 1] == 0:
            print(f"  n={n}: prev zero, undefined")
            continue
        r = abs_c_kstar[n] / abs_c_kstar[n - 1]
        print(f"  n={n}: ratio = {float(r):.6f}  (target 0.5)")
    print()

    # A2: signed sum rate
    print("## A2 — signed sum rate-0.5 decay")
    for n in [1, 2, 3, 4, 5, 6]:
        sc = sum_c_n[n]
        print(f"  n={n}: S_n^signed = Σ_k c_{{n,k}} = {float(sc):+.6e}")
    print()
    print("## A2 ratios |S_n^signed| / |S_{n-1}^signed|:")
    for n in [2, 3, 4, 5, 6]:
        if sum_c_n[n - 1] == 0:
            print(f"  n={n}: prev zero, undefined")
            continue
        r = abs(sum_c_n[n]) / abs(sum_c_n[n - 1])
        print(f"  n={n}: ratio = {float(r):.6f}  (target 0.5)")
    print()

    # F-pattern diagnostics
    print("## Sign pattern of c_{n,k} (for A2 cancellation diagnosis):")
    for n in [2, 3, 4, 5, 6]:
        signs = []
        for k in range(n):
            v = c_nk[(n, k)]
            signs.append("+" if v > 0 else ("-" if v < 0 else "0"))
        print(f"  n={n}: signs by k=0..{n-1}: {' '.join(signs)}")
    print()

    print("## F2 — is c_{n, n-1} dominant? share of |c_{n, n-1}|:")
    for n in [2, 3, 4, 5, 6]:
        last = abs(c_nk[(n, n - 1)])
        tot = sum(abs(c_nk[(n, k)]) for k in range(n))
        share = float(last / tot) if tot != 0 else 0
        print(f"  n={n}: |c_{{n,n-1}}| / Σ|c| = {share:.4f}")
    print()

    print("## Comparison: ε_n direct empirical (Σ c if decomp closes, else S_n - 7/15):")
    for n in [1, 2, 3, 4, 5, 6]:
        en = eps_check[n]
        print(f"  ε_{n} = {float(en):+.6e}")
    print("## Empirical rate ε_n / ε_{n-1}:")
    for n in [2, 3, 4, 5, 6]:
        if eps_check[n - 1] == 0:
            continue
        r = eps_check[n] / eps_check[n - 1]
        print(f"  n={n}: ε_{n}/ε_{n-1} = {float(r):+.6f}")
    print()


if __name__ == "__main__":
    main()
