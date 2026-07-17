"""
phase2_approx_dark_probe.py — Approximate-darkness leakage ratios for candidate
structural subspaces of Syracuse's adaptive Kraus family.

Per DWM_DARK_SUBSPACE_ATTACK_PLAN.md Phase 2 and PHASE1_DARK_SUBSPACE_RESULT.md
spec: for each candidate D ⊂ H_n, compute
  alpha_D(M) := ||P_D · M · P_{D^perp}|| / ||M||  (operator norm)
across the Kraus family.

Candidates at level n:
  - D_W: 3-fiber-zero-mean subspace (dim = n_dim - 2). Identified this session as K_k's
    kernel and the natural inverse-limit dark-subspace candidate.
  - D_class: complement; class-resolved functions (constant within each ξ mod 3 class). Dim 2.
  - D_T_diag: 1-dim (1, 4) eigendirection in class space, carrying T_lead's 43/45.

Structural prediction (chat reasoning):
  x_j(b_prior) = 3^{2j-2} · 2^{-b_prior} mod 3^n.
  - j = 1: x_1 ≡ ±1 mod 3, so phase mixes within 3-fibers → D_W NOT preserved → α > 0.
  - j ≥ 2: x_j ≡ 0 mod 3 (since 3^{2j-2} | x_j for j ≥ 2), so phase is constant within each
    3-fiber → D_W exactly preserved → α = 0.

Bonus: compute commutant dim restricted to {j ≥ 2}. Prediction = 2 (block-diagonal on
(D_W, D_class) at j ≥ 2; first step j=1 is the only mixing event).
"""
import sys, os, json, cmath, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_kraus_family(n, V_MAX=12, J_MAX=3):
    N = 3 ** n
    inv2 = pow(2, -1, N)
    coprime = [r for r in range(N) if r % 3 != 0]
    idx = {r: i for i, r in enumerate(coprime)}
    n_dim = len(coprime)
    two_pi_over_N = 2.0 * np.pi / N

    family = []
    seen_keys = set()
    for j in range(1, J_MAX + 1):
        b_max = 3 * V_MAX
        for b_prior in range(0, b_max + 1):
            x_j = pow(3, 2 * j - 2, N) * pow(inv2, b_prior, N) % N
            for v in range(1, V_MAX + 1):
                pow_inv2_v = pow(inv2, v, N)
                x_phase = (x_j * pow_inv2_v) % N
                key = (x_phase, v)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                M = np.zeros((n_dim, n_dim), dtype=complex)
                weight = 2.0 ** (-v / 2)
                for xi in coprime:
                    target_xi = (xi * pow_inv2_v) % N
                    if target_xi % 3 == 0 or target_xi not in idx:
                        continue
                    phase = cmath.exp(-1j * two_pi_over_N * xi * x_phase)
                    M[idx[xi], idx[target_xi]] = weight * phase
                family.append({"label": f"v={v}_j={j}_b={b_prior}_xphase={x_phase}",
                               "j": j, "v": v, "b_prior": b_prior, "x_phase": x_phase,
                               "M": M})
    return family, coprime


def build_projections(coprime):
    n_dim = len(coprime)
    class_plus = [i for i, xi in enumerate(coprime) if xi % 3 == 1]
    class_minus = [i for i, xi in enumerate(coprime) if xi % 3 == 2]
    e_plus = np.zeros(n_dim, dtype=complex)
    e_plus[class_plus] = 1.0 / np.sqrt(len(class_plus))
    e_minus = np.zeros(n_dim, dtype=complex)
    e_minus[class_minus] = 1.0 / np.sqrt(len(class_minus))
    P_class = np.outer(e_plus, e_plus.conj()) + np.outer(e_minus, e_minus.conj())
    P_W = np.eye(n_dim, dtype=complex) - P_class
    v_14 = e_plus + 4.0 * e_minus
    v_14 /= np.linalg.norm(v_14)
    P_T_diag = np.outer(v_14, v_14.conj())
    return {"P_W": P_W, "P_class": P_class, "P_T_diag": P_T_diag}


def leakage_ratio(P_D, M):
    n_dim = P_D.shape[0]
    P_perp = np.eye(n_dim, dtype=complex) - P_D
    A = P_D @ M @ P_perp
    return np.linalg.norm(A, ord=2) / np.linalg.norm(M, ord=2)


def restricted_commutant_dim(family_subset, n_dim, tol=1e-9):
    if not family_subset:
        return 0, np.array([])
    T_rows = []
    I = np.eye(n_dim, dtype=complex)
    for op in family_subset:
        M = op["M"]
        T_M = np.kron(M.T, I) - np.kron(I, M)
        T_rows.append(T_M)
    T = np.vstack(T_rows)
    U, sigma, Vh = np.linalg.svd(T, full_matrices=False)
    null_mask = sigma < tol * max(sigma.max(), 1.0)
    return int(np.sum(null_mask)), sigma


def stats(vals):
    if len(vals) == 0:
        return None
    a = np.array(vals)
    return {"n": int(len(a)), "max": float(a.max()), "mean": float(a.mean()),
            "median": float(np.median(a)), "min": float(a.min())}


def main():
    out = {}
    for n in (2, 3):
        print(f"\n=== n={n} (V_n = L²((Z/{3**n})*), dim {2*3**(n-1)}) ===")
        t0 = time.time()
        family, coprime = build_kraus_family(n, V_MAX=12, J_MAX=3)
        n_dim = len(coprime)
        print(f"  built {len(family)} distinct Kraus operators in {time.time()-t0:.2f}s")
        projs = build_projections(coprime)

        leakages = {name: {"j1": [], "j_ge_2": [], "all": []} for name in projs}
        for op in family:
            M = op["M"]
            j = op["j"]
            for name, P_D in projs.items():
                alpha = leakage_ratio(P_D, M)
                leakages[name]["all"].append(alpha)
                (leakages[name]["j1"] if j == 1 else leakages[name]["j_ge_2"]).append(alpha)

        family_j1 = [op for op in family if op["j"] == 1]
        family_j_ge_2 = [op for op in family if op["j"] >= 2]
        print(f"  family split: {len(family_j1)} j=1 operators, {len(family_j_ge_2)} j>=2 operators")

        dim_j1, _ = restricted_commutant_dim(family_j1, n_dim)
        dim_j_ge_2, _ = restricted_commutant_dim(family_j_ge_2, n_dim)
        print(f"  dim(A') restricted to j=1   : {dim_j1}")
        print(f"  dim(A') restricted to j>=2  : {dim_j_ge_2}")

        for name in ["P_W", "P_class", "P_T_diag"]:
            print(f"\n  Subspace {name}:")
            for restr in ["all", "j1", "j_ge_2"]:
                s = stats(leakages[name][restr])
                if s is None:
                    continue
                print(f"    {restr:8s}: n={s['n']:4d}, max={s['max']:.6e}, mean={s['mean']:.6e}, median={s['median']:.6e}, min={s['min']:.6e}")

        out[f"n={n}"] = {
            "n": n,
            "n_dim": int(n_dim),
            "n_kraus": len(family),
            "n_j1": len(family_j1),
            "n_j_ge_2": len(family_j_ge_2),
            "commutant_dim_j1": dim_j1,
            "commutant_dim_j_ge_2": dim_j_ge_2,
            "leakages": {name: {restr: stats(leakages[name][restr]) for restr in ["all", "j1", "j_ge_2"]} for name in projs},
        }

    with open(os.path.join(OUTDIR, "phase2_approx_dark_probe.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "phase2_approx_dark_probe.json"))


if __name__ == "__main__":
    main()
