"""
Sensitivity: how stable is lambda_max(M_closed) under perturbations?
If the "match" depends on a precise eigenvalue, robustness matters.
"""
import sys
import io
import numpy as np
from pathlib import Path

sys.stdout = io.TextIOWrapper(open(1, 'wb', 0), encoding='utf-8', write_through=True)


def build_M_variants(k, variant='default'):
    Mod = 1 << k
    inv3 = pow(3, -1, Mod)
    M = np.zeros((Mod, Mod), dtype=np.float64)

    if variant == 'default':
        # Standard Result 23 M_closed
        for r in range(Mod):
            M[r, (2*r) % Mod] += 1.0
            if r % 2 == 0:
                for j in range(3):
                    lift = r + j*Mod
                    if lift % 6 == 4:
                        r_right = ((lift - 1) * inv3) % Mod
                        M[r, r_right] += 1.0/3.0
                        break
        return M

    elif variant == 'no_inv3':
        # Pure doubling: trivially lambda = 1
        for r in range(Mod):
            M[r, (2*r) % Mod] += 1.0
        return M

    elif variant == 'full_inv3_all':
        # Inverse-3 from EVERY residue (not just EVEN): hypothetical
        for r in range(Mod):
            M[r, (2*r) % Mod] += 1.0
            for j in range(3):
                lift = r + j*Mod
                if lift % 6 == 4:
                    r_right = ((lift - 1) * inv3) % Mod
                    M[r, r_right] += 1.0/3.0
                    break
        return M

    elif variant == 'inv3_weight_1':
        # Inverse-3 weighted at 1 instead of 1/3
        for r in range(Mod):
            M[r, (2*r) % Mod] += 1.0
            if r % 2 == 0:
                for j in range(3):
                    lift = r + j*Mod
                    if lift % 6 == 4:
                        r_right = ((lift - 1) * inv3) % Mod
                        M[r, r_right] += 1.0
                        break
        return M

    elif variant == 'inv3_weight_half':
        # Inverse-3 weighted at 1/2 instead of 1/3
        for r in range(Mod):
            M[r, (2*r) % Mod] += 1.0
            if r % 2 == 0:
                for j in range(3):
                    lift = r + j*Mod
                    if lift % 6 == 4:
                        r_right = ((lift - 1) * inv3) % Mod
                        M[r, r_right] += 0.5
                        break
        return M

    raise ValueError(variant)


def lam_max(M):
    eigvals = np.linalg.eigvals(M.T)
    return float(np.max(np.abs(eigvals)))


def main():
    print("=== Sensitivity: lambda_max(M_closed) under perturbations ===\n", flush=True)
    print(f"  {'variant':<22}  {'lam_max':>10}  {'log/log2':>10}  {'2log/log2':>10}", flush=True)

    for var in ['default', 'no_inv3', 'full_inv3_all', 'inv3_weight_1', 'inv3_weight_half']:
        M = build_M_variants(8, var)
        lam = lam_max(M)
        d_natural = np.log(lam) / np.log(2)
        d_doubled = 2 * np.log(lam) / np.log(2)
        print(f"  {var:<22}  {lam:>10.6f}  {d_natural:>10.6f}  {d_doubled:>10.6f}", flush=True)

    # The big question: how much does lam_max have to shift to land at 0.6942?
    print(f"\n  Target: 2*log(lam)/log(2) = 0.6942 ⟹ lam = 2^(0.6942/2) = {2**(0.6942/2):.6f}", flush=True)
    print(f"  Target: log(lam)/log(2) = 0.6942 ⟹ lam = 2^0.6942 = {2**0.6942:.6f}", flush=True)
    print(f"  Default M_closed lam = 1.263763", flush=True)
    print(f"  Gap to first target: {1.263763 - 2**(0.6942/2):.6f}", flush=True)

    print(f"\nNote: lam_max is structurally determined by the 1+1/3 = 4/3 row sum on EVEN residues", flush=True)
    print(f"(odd-row sums are 1, even-row sums are 4/3). The Perron eigenvalue sits between", flush=True)
    print(f"these bounds: 1 < lam_max = 1.2638 < 4/3 = 1.3333. The exact value depends on the", flush=True)
    print(f"residue-graph topology which determines how mass cycles between odd/even classes.", flush=True)


if __name__ == "__main__":
    main()
