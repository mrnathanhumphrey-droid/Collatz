"""
interlevel_twisted_endomorphism.py — Build the square endomorphism on V_n via
twisted fiber selection composed with the Fourier-side Tao transfer U_n.

Setup:
- U_n : V_n -> V_{n+1} (Fourier-side Tao recursion)
- Image of U_n lies entirely in W_n (3-fiber-zero-mean subspace) by phase-cube-root cancellation.
- The 3-fiber translation Z3 acts on V_{n+1} by sending xi' -> xi' + 3^n.
  Z3 has order 3; eigenspaces are at omega in {1, omega_3, omega_3^2}.
  Eigenspace at omega=1 is T(V_n) (lift subspace).
  Eigenspaces at omega=omega_3, omega_3^2 are each (1/2)*dim(W_n) = dim(V_n).
- Define T^omega: V_{n+1} -> V_n by
    T^omega(g)(xi) := (1/3) [ g(xi) + omega * g(xi + 3^n) + omega^2 * g(xi + 2*3^n) ]
  for xi in (Z/3^n)*, where the 3 lifts of xi in (Z/3^{n+1})* are {xi, xi+3^n, xi+2*3^n}.
- For omega = 1: T^omega = (1/3) * T_sum -> gives 0 on Image(U_n) since image is in W_n.
- For omega = omega_3 or omega_3^2: T^omega projects onto the non-trivial Z3 eigenspaces of g.

The endomorphism of interest:
    Phi_omega := T^omega @ U_n : V_n -> V_n

Eigenvalues of Phi_omega can be COMPLEX. We hunt for a stable complex-conjugate pair
across n=2..5 whose argument matches the empirical period-9.2 oscillation
(theta = 2*pi/9.2 ≈ 0.683 rad) at some modulus near 0.984 or 43/45.

Cross-check: spectrum of Phi_{omega=1} should be ≈ 0 (since U_n image avoids T(V_n)).
"""
import sys, os, json, time, cmath
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"


def build_U(n, V_MAX=30):
    """U_n : V_n -> V_{n+1} (Fourier side, complex)."""
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    inv2_cod = pow(2, -1, N_cod)
    pow_inv2_cod = [pow(inv2_cod, v, N_cod) for v in range(1, V_MAX + 1)]
    Z = 1.0 - 2.0 ** (-V_MAX)
    two_pi_over_3pp = 2.0 * np.pi / N_cod
    n_dom = len(coprime_dom)
    n_cod = len(coprime_cod)
    U = np.zeros((n_cod, n_dom), dtype=complex)
    for xi_prime in coprime_cod:
        i = idx_cod[xi_prime]
        for v in range(1, V_MAX + 1):
            xi_prime_inv2v = (xi_prime * pow_inv2_cod[v - 1]) % N_cod
            xi = xi_prime_inv2v % N_dom
            if xi == 0 or xi % 3 == 0:
                continue
            j = idx_dom[xi]
            U[i, j] += (2.0 ** (-v)) / Z * cmath.exp(-1j * two_pi_over_3pp * xi_prime * pow_inv2_cod[v - 1])
    return U, coprime_dom, coprime_cod


def build_T_omega(n, omega):
    """T^omega : V_{n+1} -> V_n, twisted fiber selection.
    T^omega[xi, xi'] = (1/3) * omega^a if xi' = xi + a*3^n with xi coprime, else 0.
    """
    N_dom = 3 ** n
    N_cod = 3 ** (n + 1)
    coprime_dom = [r for r in range(N_dom) if r % 3 != 0]
    coprime_cod = [r for r in range(N_cod) if r % 3 != 0]
    idx_dom = {r: i for i, r in enumerate(coprime_dom)}
    idx_cod = {r: i for i, r in enumerate(coprime_cod)}
    T = np.zeros((len(coprime_dom), len(coprime_cod)), dtype=complex)
    for xi in coprime_dom:
        for a in (0, 1, 2):
            xi_prime = xi + a * N_dom
            if xi_prime % 3 == 0:
                continue
            T[idx_dom[xi], idx_cod[xi_prime]] = (1.0 / 3.0) * (omega ** a)
    return T


def run(n, V_MAX=30):
    print(f"\n=== n={n} (V_n dim={2*3**(n-1)}) ===")
    t0 = time.time()
    U, _, _ = build_U(n, V_MAX=V_MAX)
    omega1 = 1.0
    omega3 = cmath.exp(2j * np.pi / 3)
    omega3sq = omega3 * omega3
    T1 = build_T_omega(n, omega1)
    T_w = build_T_omega(n, omega3)
    T_w2 = build_T_omega(n, omega3sq)
    print(f"  built operators in {time.time()-t0:.2f}s")

    # Sanity: T1 @ U = 0 (T sum kills U's image since image is in W)
    Phi1 = T1 @ U
    err1 = np.linalg.norm(Phi1)
    print(f"  ||T^1 @ U_n|| (should be ~0): {err1:.2e}")

    # Twisted endomorphisms
    Phi_w = T_w @ U
    Phi_w2 = T_w2 @ U
    eigs_w = np.linalg.eigvals(Phi_w)
    eigs_w2 = np.linalg.eigvals(Phi_w2)

    # Sort by modulus, descending
    eigs_w_sorted = sorted(eigs_w, key=lambda z: -abs(z))
    eigs_w2_sorted = sorted(eigs_w2, key=lambda z: -abs(z))

    print(f"  Top 10 eigenvalues of Phi_{{omega_3}} = T^omega_3 @ U_n (complex):")
    for i, z in enumerate(eigs_w_sorted[:10]):
        mod = abs(z)
        arg = np.angle(z)
        period_n = 2 * np.pi / abs(arg) if abs(arg) > 1e-12 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f}rad  period_n={period_n:.3f}")

    print(f"  Top 10 eigenvalues of Phi_{{omega_3^2}} = T^omega_3^2 @ U_n (complex):")
    for i, z in enumerate(eigs_w2_sorted[:10]):
        mod = abs(z)
        arg = np.angle(z)
        period_n = 2 * np.pi / abs(arg) if abs(arg) > 1e-12 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f}rad  period_n={period_n:.3f}")

    return {
        "n": n,
        "T1_at_U_norm": float(err1),
        "Phi_omega_eigs_top10": [[float(z.real), float(z.imag)] for z in eigs_w_sorted[:10]],
        "Phi_omega2_eigs_top10": [[float(z.real), float(z.imag)] for z in eigs_w2_sorted[:10]],
    }


def main():
    out = {}
    for n in (2, 3, 4, 5):
        out[f"n={n}"] = run(n, V_MAX=30 if n <= 3 else 25)
    with open(os.path.join(OUTDIR, "interlevel_twisted_endomorphism.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "interlevel_twisted_endomorphism.json"))

    # Cross-n summary: top-5 |eig| of Phi_omega + their args + period_n
    print("\n--- Summary: top 5 (|eig|, arg/rad, period_n) of Phi_{omega_3} across n ---")
    for n in (2, 3, 4, 5):
        r = out[f"n={n}"]
        rows = []
        for re, im in r["Phi_omega_eigs_top10"][:5]:
            z = complex(re, im)
            mod = abs(z)
            arg = np.angle(z)
            per = 2 * np.pi / abs(arg) if abs(arg) > 1e-12 else float('inf')
            rows.append((mod, arg, per))
        s = "  ".join(f"{m:.4f}@{a:+.3f}(per{p:.2f})" for m, a, p in rows)
        print(f"  n={n}: {s}")


if __name__ == "__main__":
    main()
