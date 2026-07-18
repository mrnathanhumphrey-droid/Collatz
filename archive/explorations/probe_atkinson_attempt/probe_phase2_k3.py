"""
probe_phase2_k3.py — Phase 2: Rota-Baxter axiom verification at k=3.

Algebra R = C^18 (functions on (Z/27)* under pointwise multiplication).
This is the natural commutative semisimple algebra at k=3.

For each candidate T : R -> R, check the Rota-Baxter axiom of weight θ:
    T(u)·T(v) = T( T(u)·v + u·T(v) - θ·u·v )
on a sample of (u, v) pairs. The axiom must hold for ALL u, v in R; we
test on 8 pairs of basis elements + 2 random elements as a numerical
sanity check (not a proof, but exposes obvious failures).

Candidates:
  T_a: projection onto subset A_a = {residues r mod 27 in {1,2,4,5,7,8,10,11,13}}
       (= half the coprime residues, "lower half" by integer order)
  T_c: projection onto characters with positive real part of dominant eigenvalue
       (here: characters with chi(2) having Re > 0)
  T_d: projection onto principal characters (= chi_0 only) — trivial RB
  T_K: K_3 itself (the Tao Markov kernel) — non-projector candidate
  T_resK: (Id - K_3)^(-1) when defined — homomorphism-resolvent attempt

Output:
  construction_attempts_k3.csv - per-candidate axiom verification
  spectrum_comparison.csv - spectra for candidates that pass
"""
from __future__ import annotations

import csv
import os
import sys
from itertools import product

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT_DIR = r"C:\Collatz\probe_atkinson_attempt"
K = 3
N_MOD = 3 ** K       # 27
COPRIME = [r for r in range(N_MOD) if r % 3 != 0]  # 18 elements
N = len(COPRIME)
STATE_IDX = {r: i for i, r in enumerate(COPRIME)}


def build_K_dense(k):
    """Tao-Syracuse Markov kernel K_k on coprime classes of Z/3^k."""
    Nm = 3 ** k
    M = 1
    v = 2 % Nm
    while v != 1:
        v = (v * 2) % Nm
        M += 1
    inv2 = pow(2, -1, Nm)
    powers_inv2 = np.empty(M, dtype=np.int64)
    pi = inv2
    for vv in range(M):
        powers_inv2[vv] = pi
        pi = (pi * inv2) % Nm
    coprime = [r for r in range(Nm) if r % 3 != 0]
    n = len(coprime)
    state_idx = {r: i for i, r in enumerate(coprime)}
    weights = np.zeros(M, dtype=np.float64)
    for vv in range(min(M, 1074)):
        weights[vv] = 2.0 ** -(vv + 1)
    weights /= weights.sum()
    K_mat = np.zeros((n, n), dtype=np.float64)
    for i_r in range(n):
        r = int(coprime[i_r])
        base = (3 * r + 1) % Nm
        for vv in range(M):
            target = (base * powers_inv2[vv]) % Nm
            K_mat[i_r, state_idx[target]] += weights[vv]
    return K_mat


def axiom_residual(T, theta, u, v):
    """Compute || T(u)·T(v) - T( T(u)·v + u·T(v) - theta·u·v ) ||_inf
    where · is pointwise (Hadamard) product."""
    Tu = T @ u
    Tv = T @ v
    lhs = Tu * Tv
    rhs = T @ (Tu * v + u * Tv - theta * u * v)
    return float(np.max(np.abs(lhs - rhs)))


def test_candidate(name, T, theta=1.0):
    """Test the Rota-Baxter axiom on a sample of (u, v) pairs.
    Returns max residual and per-pair details."""
    rng = np.random.default_rng(seed=12345)
    pairs = []
    # 8 standard basis pairs
    for i in [0, 1, 2, 5, 9, 13]:
        for j in [0, 3, 7, 11, 17]:
            if (i, j) in [(0, 0), (5, 11), (9, 17)]:
                e_i = np.zeros(N); e_i[i] = 1.0
                e_j = np.zeros(N); e_j[j] = 1.0
                pairs.append((f"e_{i}, e_{j}", e_i, e_j))
    # 4 random pairs
    for s in range(4):
        u = rng.standard_normal(N)
        v = rng.standard_normal(N)
        pairs.append((f"random_{s}", u, v))
    # All-ones (identity-like) test
    one = np.ones(N)
    pairs.append(("ones, ones", one, one))

    residuals = []
    for label, u, v in pairs:
        r = axiom_residual(T, theta, u, v)
        residuals.append((label, r))

    max_res = max(r for _, r in residuals)
    return max_res, residuals


def main():
    print(f"Phase 2: RB axiom check at k={K}, N={N}")
    print()

    # Build framework operator K_3
    K_mat = build_K_dense(K)
    print(f"K_3 built: shape {K_mat.shape}, row sums {K_mat.sum(axis=1)[:3]} ... 1.0")
    print()

    # === Candidates ===
    candidates = []

    # T_a: projection onto subset A = first 9 residues by index
    A_a = list(range(9))
    T_a = np.zeros((N, N))
    for i in A_a:
        T_a[i, i] = 1.0
    candidates.append(("T_a (proj on first 9 residues by COPRIME index)", T_a))

    # T_c: Atkinson-style "positive half" projector — by integer ordering
    # of residues. COPRIME = [1,2,4,5,7,8,...,26]. Pick r < 14: that's roughly half.
    A_c = [i for i, r in enumerate(COPRIME) if r < 14]  # r in {1,2,4,5,7,8,10,11,13} = 9 residues
    T_c = np.zeros((N, N))
    for i in A_c:
        T_c[i, i] = 1.0
    candidates.append((f"T_c (proj on r < 14: |A|={len(A_c)})", T_c))

    # T_d: pick a different half: r mod 3 == 1
    A_d = [i for i, r in enumerate(COPRIME) if r % 3 == 1]  # 9 residues
    T_d = np.zeros((N, N))
    for i in A_d:
        T_d[i, i] = 1.0
    candidates.append((f"T_d (proj on r ≡ 1 mod 3: |A|={len(A_d)})", T_d))

    # T_K: K_3 itself
    candidates.append(("T_K (K_3 Markov kernel, not projector)", K_mat))

    # T_resK: (I - K)^-1 IF nonsingular. K has eigenvalue 1 (stationary), so
    # I - K is singular. Use pseudo-inverse on (I - K).
    IK = np.eye(N) - K_mat
    rank = np.linalg.matrix_rank(IK, tol=1e-10)
    print(f"rank(I - K_3) = {rank} (= N-1 = 17 expected since K stochastic, "
          f"shows pseudoinverse needed)")
    T_resK = np.linalg.pinv(IK)
    candidates.append(("T_resK (pseudoinverse of I - K_3)", T_resK))

    # T_K_squared: K^2 (still stochastic, equally non-RB-natural)
    candidates.append(("T_K2 (K_3^2)", K_mat @ K_mat))

    # T_chi: Fourier-projector — keep only characters of order dividing 9
    # The character group of (Z/27)* is Z/18 (cyclic). Characters of order |9
    # are characters of the quotient Z/9 ⊂ Z/18: 9 of them. Their orthogonal
    # complement is the 9 characters of order 18 / non-trivial-square.
    # Build character table:
    coprime = COPRIME
    # (Z/27)* is cyclic of order 18; find a primitive root g.
    # g=2 is a primitive root mod 27 (well-known).
    g = 2
    # Index each coprime residue by its discrete log in Z/18.
    discrete_log = {}
    cur = 1
    for j in range(18):
        if cur in STATE_IDX:
            discrete_log[cur] = j
        cur = (cur * g) % 27
    # Character chi_a(r) = exp(2πi · a · log_g(r) / 18) for a ∈ Z/18.
    # In real basis: real and imaginary parts of chi_a give an
    # 18-dim representation. Use complex characters directly for the projector.
    omega = np.exp(2j * np.pi / 18)
    char_table = np.zeros((18, N), dtype=complex)
    for a in range(18):
        for i, r in enumerate(coprime):
            char_table[a, i] = omega ** (a * discrete_log[r])
    # T_chi: project onto characters with a in {0, 2, 4, 6, 8, 10, 12, 14, 16}
    # (= even-order elements of Z/18 = 9 characters)
    a_set = list(range(0, 18, 2))
    P = np.zeros((N, N), dtype=complex)
    for a in a_set:
        v = char_table[a] / np.sqrt(N)  # orthonormalize per-character
        P += np.outer(v, v.conj())
    # Take real part (should be exactly real if a_set is closed under conjugation
    # a -> -a mod 18 = 18-a, and our set {0,2,4,...,16} is closed since
    # 18-2=16, 18-4=14, etc. For a=0: 18-0 = 18 = 0. ✓)
    T_chi = P.real
    candidates.append(("T_chi (Fourier proj on a in {0,2,...,16})", T_chi))

    # Now test all candidates
    rows = []
    for name, T in candidates:
        # Theta = 1 for projectors and Wiener-Hopf-ish operators; try other theta too
        for theta in [0.0, 1.0]:
            max_res, details = test_candidate(name, T, theta=theta)
            # Spectrum
            try:
                eigs = np.linalg.eigvals(T)
                eigs_sorted = sorted(eigs, key=lambda x: -abs(x))
                lam_top = abs(eigs_sorted[0])
                lam_2 = abs(eigs_sorted[1]) if len(eigs_sorted) > 1 else float("nan")
                # Count distinct |lambda| values
                eig_mags = sorted(set(round(abs(e), 6) for e in eigs), reverse=True)
                n_distinct = len(eig_mags)
            except Exception as e:
                lam_top = float("nan"); lam_2 = float("nan")
                n_distinct = -1
                eig_mags = []
            # Is it a projector? T^2 - T should be ~0
            proj_resid = float(np.max(np.abs(T @ T - T)))
            print(f"  {name} | theta={theta} | max_res = {max_res:.4e} | "
                  f"lam_top = {lam_top:.4f} | proj_resid = {proj_resid:.2e}")
            rows.append({
                "candidate": name,
                "theta": theta,
                "max_residual": max_res,
                "passes_RB": max_res < 1e-8,
                "is_projector": proj_resid < 1e-10,
                "lambda_top": lam_top,
                "lambda_2": lam_2,
                "n_distinct_eigvals": n_distinct,
                "spectrum_top5": ",".join(f"{m:.6f}" for m in eig_mags[:5]),
            })

    # Write CSV
    csv_path = os.path.join(OUT_DIR, "construction_attempts_k3.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        cols = ["candidate", "theta", "max_residual", "passes_RB",
                "is_projector", "lambda_top", "lambda_2",
                "n_distinct_eigvals", "spectrum_top5"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: (f"{v:.10e}" if isinstance(v, float) else v)
                        for c, v in r.items()})
    print(f"\n[csv: {csv_path}]")

    # Find passing candidates and analyze spectrum further
    print()
    print("=" * 78)
    print("SPECTRA OF RB-PASSING CANDIDATES")
    print("=" * 78)

    spec_rows = []
    for r in rows:
        if r["passes_RB"]:
            print(f"\nPasser: {r['candidate']} (θ={r['theta']})")
            # For projectors, spectrum is {0, 1}; nothing to see.
            # For non-projectors, list full spectrum.
            T_passer = next(T for nm, T in candidates if nm == r["candidate"])
            eigs = np.linalg.eigvals(T_passer)
            sorted_eigs = sorted(eigs, key=lambda x: -abs(x))
            for j, e in enumerate(sorted_eigs):
                spec_rows.append({
                    "candidate": r["candidate"],
                    "theta": r["theta"],
                    "rank": j + 1,
                    "lambda_real": float(np.real(e)),
                    "lambda_imag": float(np.imag(e)),
                    "abs_lambda": float(abs(e)),
                })
            print(f"  unique |lambda| values:")
            mags = sorted(set(round(abs(e), 6) for e in eigs), reverse=True)
            for m in mags:
                count = sum(1 for e in eigs if round(abs(e), 6) == m)
                print(f"    |λ| = {m:.6f}: multiplicity {count}")

            # Check 1-q^n match
            non_zero_one = [m for m in mags if 1e-8 < m < 1 - 1e-8]
            if non_zero_one:
                print(f"  spectrum has nontrivial values: {non_zero_one}")
                # Check geometric pattern 1-q^n
                if len(non_zero_one) >= 2:
                    qs = []
                    for n_test in [1, 2, 3, 4]:
                        for m in non_zero_one:
                            if 0 < m < 1:
                                q = (1 - m) ** (1.0 / n_test)
                                qs.append((m, n_test, q))
                    print(f"  candidate q from 1-q^n inversion: {qs[:10]}")

    spec_csv = os.path.join(OUT_DIR, "spectrum_comparison.csv")
    if spec_rows:
        with open(spec_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(spec_rows[0].keys()))
            w.writeheader()
            for r in spec_rows:
                w.writerow({c: (f"{v:.10e}" if isinstance(v, float) else v)
                            for c, v in r.items()})
        print(f"\n[csv: {spec_csv}]")

    # Repunit position analysis for k=3
    print()
    print("=" * 78)
    print("REPUNIT POSITION ANALYSIS at k=3")
    print("=" * 78)
    # Repunit at k=3: r = (4^3 - 1)/3 = 21
    print(f"Repunit residue (4^3-1)/3 = 21, mod 27 = 21.")
    if 21 in STATE_IDX:
        idx_21 = STATE_IDX[21]
        print(f"  21 ∈ COPRIME at index {idx_21}")
        print(f"  21 mod 3 = {21 % 3}")
        print(f"  discrete_log_2(21) mod 18 = {discrete_log[21]}")
        # Position in each candidate projector's image set
        for name, T in candidates:
            if T.shape == (N, N):
                # is e_{idx_21} in image of T? (T @ e_{21} mostly e_{21}?)
                e21 = np.zeros(N); e21[idx_21] = 1.0
                Te21 = T @ e21
                # If projector, T @ e21 = e21 (in image) or 0 (not in image)
                in_image = np.allclose(Te21, e21, atol=1e-10)
                not_in_image = np.allclose(Te21, np.zeros(N), atol=1e-10)
                if in_image:
                    print(f"  {name}: 21 IN T's image")
                elif not_in_image:
                    print(f"  {name}: 21 NOT in T's image")
                else:
                    print(f"  {name}: 21 mixed (T·e_21 = {Te21[idx_21]:.4f} on diag)")
    else:
        print(f"  21 not in COPRIME (shouldn't happen since 21 mod 3 = 0)")
        # Wait — 21 mod 3 = 0, so 21 is NOT in (Z/27)* coprime classes!
        # Let me reconsider: (4^k-1)/3 for k=3 is 21, but 21 = 3·7 so NOT coprime to 3!
        # The repunit at k=k is 1, 5, 21, 85, 341 — these are the residues mod 2^k,
        # not mod 3^k. So the "repunit at k=3 mod 3^3" framing in the brief
        # is conflating two structures.

    print(f"\n*** NOTE: the repunit (4^k-1)/3 = 21 (at k=3) is mod 2^k, "
          f"NOT mod 3^k. ***")
    print(f"  21 mod 2^5 = 21 (the repunit at modulus 2^5 = 32)")
    print(f"  21 mod 3^3 = 21, but 21 = 3·7 is NOT coprime to 3, so 21 is NOT "
          f"a state in the K_k Markov chain on (Z/3^k)*.")
    print(f"  The 'repunit connection' from the framework is a 2-adic feature,")
    print(f"  not a 3-adic feature, so it doesn't directly index into the")
    print(f"  character algebra of (Z/3^k)* used in this Atkinson construction.")


if __name__ == "__main__":
    main()
