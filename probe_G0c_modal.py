"""
G0c modal side -- the INDEPENDENT (M) coordinate done right: the AMPLITUDE-CARRYING subdominant
eigenvalue of build_M (R26/R32), not the raw |l2|/l1 (which is the within-cell/tower mode ~0.98).

||pi_k||^2 = 1^T M^k v0 = Sum_i A_i lambda_i^k,  A_i = (1^T r_i)(l_i^T v0)/(l_i^T r_i),
r_i/l_i right/left eigenvectors. r_q = |lambda_i|/lambda_1 for the SUBDOMINANT i with |A_i| NOT
negligible (tower/within modes have A_i ~ 0). This is the modal coordinate that must match the
gate-validated cross-rho (probe_27): r_5~0.62, r_7~0.39.

Dense eig where feasible (L=1 all q; L=2 q=3). Report the amplitude-ranked spectrum.
"""
import numpy as np
from probe_25_transfer_operator_Aprime import build_M
from probe_6_conservation_generalize import order_of_two

LOG = []


def log(m=""):
    print(m); LOG.append(str(m))


def modal_amplitudes(q, L):
    M, idx, n = build_M(q, L)
    A = M.toarray()
    w, R = np.linalg.eig(A)                 # right eigvecs (columns)
    Rinv = np.linalg.pinv(R)                # pinv: build_M is DEFECTIVE at q=3 (R39 Jordan)
    v0 = np.zeros(n); v0[idx[(1, 1, 0)]] = 1.0
    ones = np.ones(n)
    # A_i = (1^T r_i) * (l_i^T v0)   [l_i^T r_i = 1 by construction of Rinv]
    coeff_right = ones @ R                  # 1^T r_i  (length n)
    coeff_left = Rinv @ v0                  # l_i^T v0 (length n)
    Amp = coeff_right * coeff_left          # A_i
    idxs = np.argsort(-np.abs(w))
    return w[idxs], Amp[idxs]


def main():
    log("# G0c modal -- amplitude-carrying subdominant of build_M vs cross-rho (probe_27) r_5~0.62 r_7~0.39")
    log("")
    PLAN = {3: [1, 2], 5: [1], 7: [1]}      # dense-feasible
    for q in PLAN:
        log(f"## q={q}")
        for L in PLAN[q]:
            w, Amp = modal_amplitudes(q, L)
            lam1 = abs(w[0])
            Atot = np.abs(Amp).sum()
            log(f"   L={L}: lambda_1={lam1:.5f} (A_1={Amp[0].real:+.4f})")
            log(f"   {'|lambda|':>10} {'|lam|/lam1':>10} {'A_i':>12} {'|A_i|/sumA':>11}  amplitude-carrying?")
            shown = 0
            r_amp = None
            for i in range(1, len(w)):
                a = abs(Amp[i])
                frac = a / Atot
                carrying = frac > 1e-4
                if carrying and r_amp is None:
                    r_amp = abs(w[i]) / lam1
                if shown < 8 or carrying:
                    log(f"   {abs(w[i]):>10.5f} {abs(w[i])/lam1:>10.5f} {Amp[i].real:>+12.5f} "
                        f"{frac:>11.2e}  {'YES <== r_q' if (carrying and abs(abs(w[i])/lam1 - (r_amp or 0))<1e-9) else ('yes' if carrying else '')}")
                    shown += 1
                if shown >= 14:
                    break
            log(f"   => amplitude-carrying subdominant r_q(L={L}) = {r_amp:.4f}" if r_amp else "   => none found")
            log("")
    log("## READ: amplitude-carrying subdominant should approach cross-rho r_q (0.62/0.39) as L grows.")
    log("   Raw |l2|/l1 (~0.98) is a ZERO-AMPLITUDE within-cell/tower mode -- correctly excluded here.")
    with open("result_G0c_modal_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
