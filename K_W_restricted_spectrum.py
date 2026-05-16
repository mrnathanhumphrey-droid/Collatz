"""
K_W_restricted_spectrum.py — Compute the spectrum of the Syracuse Markov chain K_k
restricted to the orthogonal complement subspace W_{k-1} of the lift T(V_{k-1}) in V_k.

Setup (per R77.5):
- V_k = R^{N_k}, N_k = 2*3^{k-1} (coprime states mod 3^k)
- Lift T : V_{k-1} -> V_k is isometric up to sqrt(3): T(u)(r') = u(r' mod 3^{k-1}) / 3
- W_{k-1} = T(V_{k-1})^perp = { f in V_k : sum over each 3-fiber r' lifts of r is 0 }
- dim W_{k-1} = 4*3^{k-2}

Question: K_k preserves T(V_{k-1}) by marginal consistency. Does K_k's action on W_{k-1}
have a stable eigenvalue near 1/2?

If yes -> finite-truncation T_M = K_k|_W carries the rate-1/2 mode (R77.2 conjecture closes).
If no  -> rate-1/2 is structurally absent at finite truncation; lives in inverse-limit (per R77.5 reframing).

We compute eigenvalues of:
  K_W := P_W^T @ K_k @ P_W
where P_W is the (N_k x dim(W)) orthonormal basis of W_{k-1}, and K_k acts on functions
via (K f)(r) = sum_{r'} K_k(r, r') f(r') (so left-multiplication on column vectors).

Note: K_k is row-stochastic; eigenvalue 1 corresponds to the constant function 1 (not in W).

Output: top 5 absolute-eigenvalues at each k = 2, 3, 4, 5.
"""
import sys, os, json, time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
OUTDIR = r"C:\Collatz\experiments_output"
os.makedirs(OUTDIR, exist_ok=True)


def build_K(k, V_MAX=40):
    """Row-stochastic Syracuse Markov transition K_k on coprime states mod 3^k.

    State r -> ((3r+1) * 2^{-v}) mod 3^k with weight 2^{-v} / Z, v = 1..V_MAX, Z = 1 - 2^{-V_MAX}.
    """
    N = 3 ** k
    inv2 = pow(2, -1, N)
    pow_inv2 = [pow(inv2, v, N) for v in range(1, V_MAX + 1)]
    coprime = [r for r in range(N) if r % 3 != 0]
    state_idx = {r: i for i, r in enumerate(coprime)}
    n = len(coprime)
    Z = 1.0 - 2.0 ** (-V_MAX)
    K = np.zeros((n, n))
    for r in coprime:
        i = state_idx[r]
        for v in range(1, V_MAX + 1):
            weight = (2.0 ** (-v)) / Z
            target = ((3 * r + 1) * pow_inv2[v - 1]) % N
            j = state_idx.get(target)
            if j is not None:
                K[i, j] += weight
    return K, coprime, state_idx


def build_W_basis(k):
    """Return an N_k x dim(W) orthonormal basis matrix for W_{k-1} inside V_k.

    W_{k-1} = { f in V_k : for each r in (Z/3^{k-1})*, sum over 3 lifts r' -> r equals 0 }.

    For each level-(k-1) coprime residue r, the three lifts in (Z/3^k)* are
        r, r + 3^{k-1}, r + 2*3^{k-1}
    (all coprime to 3 since gcd(r, 3) = 1 and 3^{k-1} is divisible by 3).
    Two orthonormal mean-zero vectors per fiber:
        v1 = (1, -1, 0) / sqrt(2)
        v2 = (1, 1, -2) / sqrt(6)
    """
    assert k >= 2
    N_kp = 3 ** k
    N_km = 3 ** (k - 1)
    coprime_kp = [r for r in range(N_kp) if r % 3 != 0]
    state_idx_kp = {r: i for i, r in enumerate(coprime_kp)}
    coprime_km = [r for r in range(N_km) if r % 3 != 0]
    n_kp = len(coprime_kp)
    cols = []
    for r in coprime_km:
        lifts = [r, r + N_km, r + 2 * N_km]
        idxs = [state_idx_kp[l] for l in lifts]
        v1 = np.zeros(n_kp)
        v1[idxs[0]] = 1.0 / np.sqrt(2.0)
        v1[idxs[1]] = -1.0 / np.sqrt(2.0)
        v2 = np.zeros(n_kp)
        v2[idxs[0]] = 1.0 / np.sqrt(6.0)
        v2[idxs[1]] = 1.0 / np.sqrt(6.0)
        v2[idxs[2]] = -2.0 / np.sqrt(6.0)
        cols.append(v1)
        cols.append(v2)
    B = np.column_stack(cols)
    # Sanity: B^T B should be identity, and B has dim 4*3^{k-2} columns
    assert B.shape == (n_kp, 4 * 3 ** (k - 2)), f"got shape {B.shape}"
    err = np.linalg.norm(B.T @ B - np.eye(B.shape[1]))
    assert err < 1e-10, f"W basis not orthonormal: err {err}"
    return B


def project_onto_W(K, B):
    """Compute B^T K B. K acts on column vectors as (K f)(r) = sum_r' K[r, r'] f(r'),
    so K is applied via left-multiplication; the restriction is B^T @ K @ B."""
    return B.T @ K @ B


def run(k, V_MAX=40):
    t0 = time.time()
    K, coprime, _ = build_K(k, V_MAX=V_MAX)
    B = build_W_basis(k)
    K_W = project_onto_W(K, B)
    t_build = time.time() - t0
    eigs = np.linalg.eigvals(K_W)
    abs_eigs = np.sort(np.abs(eigs))[::-1]
    print(f"k={k}: N_k={len(coprime)}, dim(W)={B.shape[1]}, t_build={t_build:.2f}s")
    print(f"  top 8 |eigs|: {abs_eigs[:8]}")
    # Full K spectrum for cross-check
    eigs_full = np.linalg.eigvals(K)
    abs_eigs_full = np.sort(np.abs(eigs_full))[::-1]
    print(f"  full K |eigs| top 8: {abs_eigs_full[:8]}")
    return {
        "k": k,
        "N_k": len(coprime),
        "dim_W": int(B.shape[1]),
        "V_MAX": V_MAX,
        "K_W_abs_eigs_top": [float(x) for x in abs_eigs[:10]],
        "K_full_abs_eigs_top": [float(x) for x in abs_eigs_full[:10]],
        "K_W_eigs": [[float(x.real), float(x.imag)] for x in eigs],
    }


def main():
    out = {}
    for k in (2, 3, 4, 5):
        out[f"k={k}"] = run(k)
    with open(os.path.join(OUTDIR, "K_W_restricted_spectrum.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote", os.path.join(OUTDIR, "K_W_restricted_spectrum.json"))
    print("\n--- Summary table: top 5 |eigs| of K_k restricted to W_{k-1} ---")
    print(f"{'k':>3} {'N_k':>5} {'dim_W':>6} {'|e1|':>9} {'|e2|':>9} {'|e3|':>9} {'|e4|':>9} {'|e5|':>9}")
    for k in (2, 3, 4, 5):
        r = out[f"k={k}"]
        e = r["K_W_abs_eigs_top"]
        print(f"{k:>3} {r['N_k']:>5} {r['dim_W']:>6} " + " ".join(f"{x:>9.5f}" for x in e[:5]))


if __name__ == "__main__":
    main()
