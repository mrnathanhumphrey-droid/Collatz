"""
T_M_class_mod9_spectrum.py — Option III probe: project T_M_tensor onto mod-9 class-pair
space (36-dim) and compute spectrum.

T_lead's 43/45 emerges from the (P_+, P_-) mod-3 class projection (2-dim, after R76
conservation kills off-diagonal classes). The mod-9 refinement gives 36 classes
(ξ_1, ξ_2) ∈ (Z/9)*², which subsumes the mod-3 structure but may carry richer
spectral structure including CC pairs.

Construction:
  V_n^class_36 = span over (c_1, c_2) ∈ (Z/9)*² of class-indicator vectors:
    E^{(c_1, c_2)} in V_n ⊗ V_n*, with (E^{(c_1,c_2)})(xi_1, xi_2) = 1 if xi_1 ≡ c_1 mod 9
                                                                   and xi_2 ≡ c_2 mod 9,
                                                                  0 otherwise.
  Normalize so each E has unit ℓ² norm.

  Bilinear T_M: V_n ⊗ V_n* -> V_{n+1} ⊗ V_{n+1}* via U_n ⊗ conj(U_n).
  Project: V_{n+1}^class_36 -> V_n^class_36 (identical class space since 9 | 3^n for n≥2).
  Endomorphism on V_n^class_36 dim 36.

Compute spectrum at n=2, 3, 4. Look for:
  1. 43/45 eigenvalue (T_lead embedded)
  2. CC pair matching empirical period 9.2 in n
  3. Sub-leading structure
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
    U = np.zeros((len(coprime_cod), len(coprime_dom)), dtype=complex)
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


def build_class_projection(coprime, mod_base=9):
    """Build the projection matrix P from V_n (coprime states mod 3^n) onto
    V_n^class_{|mod_base|} = function space on (Z/mod_base)*.

    P[c, xi] = (1 / N_c) if xi ≡ c (mod mod_base), else 0.
    Shape: (num_classes, len(coprime)).
    """
    coprime_classes = [c for c in range(mod_base) if c % 3 != 0]
    num_classes = len(coprime_classes)
    idx_class = {c: i for i, c in enumerate(coprime_classes)}
    P = np.zeros((num_classes, len(coprime)))
    count_per_class = np.zeros(num_classes)
    for j, xi in enumerate(coprime):
        c = xi % mod_base
        if c not in idx_class:
            continue
        i = idx_class[c]
        P[i, j] = 1.0
        count_per_class[i] += 1
    # Normalize each row by class count
    for i in range(num_classes):
        if count_per_class[i] > 0:
            P[i, :] /= count_per_class[i]
    return P, coprime_classes


def build_class_injection(coprime, mod_base=9):
    """Build the injection matrix I from V_n^class to V_n.

    I[xi, c] = 1 if xi ≡ c (mod mod_base), else 0.
    Shape: (len(coprime), num_classes).
    """
    coprime_classes = [c for c in range(mod_base) if c % 3 != 0]
    idx_class = {c: i for i, c in enumerate(coprime_classes)}
    I = np.zeros((len(coprime), len(coprime_classes)))
    for j, xi in enumerate(coprime):
        c = xi % mod_base
        if c in idx_class:
            I[j, idx_class[c]] = 1.0
    return I


def run(n, mod_base=9, V_MAX=30):
    print(f"\n=== n={n}, mod_base={mod_base} ===")
    t0 = time.time()
    U, coprime_dom, coprime_cod = build_U(n, V_MAX=V_MAX)
    P_dom, classes_dom = build_class_projection(coprime_dom, mod_base=mod_base)
    P_cod, classes_cod = build_class_projection(coprime_cod, mod_base=mod_base)
    I_dom = build_class_injection(coprime_dom, mod_base=mod_base)
    I_cod = build_class_injection(coprime_cod, mod_base=mod_base)
    num_classes = len(classes_dom)
    print(f"  num_classes = {num_classes}, dim V_n = {len(coprime_dom)}, dim V_{n+1} = {len(coprime_cod)}")
    print(f"  classes = {classes_dom}")

    # Method A: project U_n onto class space (V_n^class -> V_{n+1}^class)
    # Class-resolved U_n: lift class vector via I_dom, apply U_n, project via P_cod.
    U_class = P_cod @ U @ I_dom  # shape (num_classes, num_classes)

    # Both U_class and conj(U_class) compose to give bilinear T_M restricted to class-pair.
    # For the bilinear pair-correlation operator on V_n^class_pair (dim num_classes^2):
    # T_M_class_pair[(c1', c2'), (c1, c2)] := U_class[c1', c1] * conj(U_class)[c2', c2]
    # (Kronecker product)
    n_cl = num_classes
    T_M_kron = np.kron(U_class, np.conj(U_class))
    print(f"  built T_M class kron, shape {T_M_kron.shape} in {time.time()-t0:.2f}s")

    # Compute spectrum
    eigs = np.linalg.eigvals(T_M_kron)
    eigs_sorted = sorted(eigs, key=lambda z: -abs(z))
    print(f"  top 15 |eigs| of T_M_class_pair (mod {mod_base} = {n_cl} classes -> {n_cl**2} pair-classes):")
    for i, z in enumerate(eigs_sorted[:15]):
        mod = abs(z)
        arg = np.angle(z)
        per = 2 * np.pi / abs(arg) if abs(arg) > 1e-9 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f}  period_n={per:.3f}")

    # Find closest to 43/45 and to 0.984
    target_43_45 = 43.0 / 45.0
    target_0_984 = 0.984
    closest_43 = min(eigs, key=lambda z: abs(abs(z) - target_43_45))
    closest_984 = min(eigs, key=lambda z: abs(abs(z) - target_0_984))
    print(f"  closest |z| to 43/45 = {target_43_45:.6f}: {abs(closest_43):.6f} (eigvalue {closest_43.real:+.4f}{closest_43.imag:+.4f}j, arg {np.angle(closest_43):+.4f})")
    print(f"  closest |z| to 0.984: {abs(closest_984):.6f} (eigvalue {closest_984.real:+.4f}{closest_984.imag:+.4f}j, arg {np.angle(closest_984):+.4f})")

    # Look for stable angle near 2π/9.2 ≈ 0.683 rad in eigvals with substantial modulus
    period_target = 9.2
    arg_target = 2 * np.pi / period_target
    candidates = [(z, abs(np.angle(z)) - arg_target) for z in eigs if abs(z) > 0.1]
    candidates.sort(key=lambda t: abs(t[1]))
    print(f"  eigvalues with |z| > 0.1 and arg closest to ±{arg_target:.4f} rad (period 9.2):")
    for i, (z, gap) in enumerate(candidates[:5]):
        mod = abs(z)
        arg = np.angle(z)
        per = 2 * np.pi / abs(arg) if abs(arg) > 1e-9 else float('inf')
        print(f"    [{i}] {z.real:+.6f}{z.imag:+.6f}j  |z|={mod:.6f}  arg={arg:+.4f} period_n={per:.3f} (Δ_arg={gap:+.4f})")

    return {
        "n": n,
        "mod_base": mod_base,
        "num_classes": int(num_classes),
        "U_class_top_singular": [float(s) for s in np.linalg.svd(U_class, compute_uv=False)[:10]],
        "top10_complex_eigs": [[float(z.real), float(z.imag)] for z in eigs_sorted[:10]],
        "closest_to_43_45": [float(closest_43.real), float(closest_43.imag), float(abs(closest_43))],
        "closest_to_0_984": [float(closest_984.real), float(closest_984.imag), float(abs(closest_984))],
    }


def main():
    out = {}
    for n in (2, 3, 4):
        out[f"n={n}, mod=9"] = run(n, mod_base=9, V_MAX=30 if n <= 3 else 25)
    # Also try mod=27 at n=3 to refine
    out["n=3, mod=27"] = run(3, mod_base=27, V_MAX=30)
    with open(os.path.join(OUTDIR, "T_M_class_mod9_spectrum.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "T_M_class_mod9_spectrum.json"))


if __name__ == "__main__":
    main()
