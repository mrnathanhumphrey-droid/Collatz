"""
PROBE 26 -- ROUTE 1 (option 1): AMPLITUDE-RESOLVED extraction of r_q from the R25 operator.

PRE-REGISTRATION (written before running).
------------------------------------------------------------------
R25 built a GATE-VALIDATED operator M_L (sum(M^k v0)=||pi_k||^2 exact for k<=L) but the
raw lambda_2/lambda_1 did NOT give r_q -- the subdominant spectrum is a tower/within
cluster near 1/3, and the cross rate is buried. FIX: use AMPLITUDES.

  g(k) = sum(M^k v0) = ||pi_k||^2 = sum_i A_i lambda_i^k   (A_i = modal amplitude in v0).
  h(k) = g(k) / lambda_1^k = sum_i A_i mu_i^k,   mu_i = lambda_i/lambda_1,  mu_1 = 1.
  c_k = h(k) - h(k-1) = sum_{i>=2} A_i mu_i^{k-1}(mu_i - 1).
  rho_k = c_{k+1}/c_k  ->  the LARGEST |mu_i| (i>=2) with A_i != 0  =:  r_q(L).
The increment kills the mu_1=1 constant; modes with ZERO amplitude never appear -> this
AUTOMATICALLY skips tower/within modes that v0 does not excite. r_q = lim rho_k.

Two readouts:
  (a) ITERATED rho_k(L): iterate M to large k (cheap matvec), rho_k -> r_q(L). Robust,
      no eigenvectors. For k<=L this MUST reproduce R22/R23 exact rho (gate).
  (b) MODAL: top eigenpairs (right r_i via eigs, left l_i via eigs(M^T)); amplitude
      A_i = (1.r_i)(l_i.v0)/(l_i.r_i); mu_i = lambda_i/lambda_1. Print modes by mu_i with
      A_i so we SEE which eigenvalue carries cross (largest mu_i<1 with |A_i| non-negligible).

HYPOTHESES / GATES:
  G_MATCH (gate): iterated rho_k for k<=L equals R22/R23 exact rho to <1e-6. If FALSE the
      normalization/iteration is wrong -> STOP.
  R_Q(L): the settled rho_k = r_q(L). Report; check convergence L=1->2 (q>=5), 1->2->3 (q=3);
      r_3(L)->1; consistency with R22 exact rho (q=5 climbing 0.53->0.62->...).
  AMP: does the modal readout show the near-1/3 cluster carrying SMALL amplitude (tower)
      while the true cross mode (larger |mu-1| weight) sits lower? PRIOR: yes (that is why
      raw lambda_2 failed). Reported as measurement.

  NOT predicted: the value r_q (priors 0-for-8). Only the structure (r_3->1, r_q<1 q>=5,
  L-convergence) is committed.

BUDGET: reuse build_M (R25). Iteration is free. eigs top-~30 for modal. q=3 L<=3, q=5/7
  L<=2, q=11 L=1 (dims 12..21609). No new heavy compute.

NOT AT STAKE: R10-R25, R5's rate, R6, R7, R12, THEOREM_C_745.
"""
import sys
import numpy as np
import scipy.sparse.linalg as spla

from probe_6_conservation_generalize import order_of_two, stationary
from probe_8_selfsimilar_overlap import sum_p2_exact
from probe_15_tower_k_count import ratio_within
from probe_25_transfer_operator_Aprime import build_M

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def perron(M):
    n = M.shape[0]
    if n <= 400:
        w, V = np.linalg.eig(M.toarray())
        i = int(np.argmax(w.real))
        return float(w[i].real)
    val = spla.eigs(M, k=1, which='LM', return_eigenvectors=False)
    return float(np.real(val[0]))


def iterate_g(M, idx, K):
    n = M.shape[0]
    v = np.zeros(n)
    v[idx[(1, 1, 0)]] = 1.0
    g = []
    for _ in range(K):
        v = M.dot(v)
        g.append(v.sum())
    return g


def modal(M, v0idx, howmany=30):
    """Top eigenpairs + amplitudes A_i and mu_i. Dense if small, else sparse eigs."""
    n = M.shape[0]
    v0 = np.zeros(n); v0[v0idx] = 1.0
    one = np.ones(n)
    if n <= 400:
        w, R = np.linalg.eig(M.toarray())
        wl, Lft = np.linalg.eig(M.toarray().T)
    else:
        m = min(howmany, n - 2)
        w, R = spla.eigs(M, k=m, which='LM')
        wl, Lft = spla.eigs(M.T.tocsr(), k=m, which='LM')
    l1 = w[int(np.argmax(w.real))].real
    # match left to right by nearest eigenvalue
    out = []
    used = set()
    for i in range(len(w)):
        j = min((jj for jj in range(len(wl)) if jj not in used),
                key=lambda jj: abs(wl[jj] - w[i]), default=None)
        if j is None:
            continue
        used.add(j)
        ri = R[:, i]; li = Lft[:, j]
        denom = li.conj().dot(ri)
        if abs(denom) < 1e-14:
            continue
        Ai = (one.dot(ri)) * (li.conj().dot(v0)) / denom
        out.append((w[i], Ai))
    out.sort(key=lambda t: -abs(t[0]))
    return l1, out


def exact_rho_ref(q, L):
    """Real rho_k for k<=L+1 from stationary, to gate the iterated rho."""
    P2 = lambda k: float(sum_p2_exact(order_of_two(q ** k)))
    cr = {}
    for k in range(1, L + 3):
        if (q - 1) * q ** (k - 1) > 200_000:
            break
        pi, cp, _ = stationary(q, k)
        cr[k] = float(np.dot(pi, pi)) / (P2(k) ** k) - 1 - ratio_within(q, k)
    ck = {k: cr[k] - cr.get(k - 1, 0.0) for k in cr if k >= 2}
    rho = {k: ck[k + 1] / ck[k] for k in ck if (k + 1) in ck and ck[k] != 0}
    return rho


def main():
    log("# PROBE 26 -- amplitude-resolved r_q from the R25 operator")
    log("# Pre-reg: G_MATCH(gate) / R_Q(L) iterated / AMP modal readout")
    log("")

    PLAN = {3: [1, 2, 3], 5: [1, 2], 7: [1, 2], 11: [1]}
    DIMCAP = 200_000
    K = 60
    rq_by_L = {}

    for q in [3, 5, 7, 11]:
        d = order_of_two(q)
        ref = exact_rho_ref(q, max(PLAN[q]))
        for L in PLAN[q]:
            dim = (d * q ** (L - 1)) ** 2 * q ** L
            if dim > DIMCAP:
                log(f"## q={q} L={L}: dim={dim} > cap -- SKIPPED (SAID)")
                continue
            M, idx, n = build_M(q, L)
            l1 = perron(M)
            g = iterate_g(M, idx, K)
            h = [g[k] / (l1 ** (k + 1)) for k in range(K)]      # h[k-1] = h(k), k=1..K
            # c_k = h(k)-h(k-1); rho_k = c_{k+1}/c_k
            c = {k: h[k - 1] - (h[k - 2] if k >= 2 else 0.0) for k in range(1, K + 1)}
            rho = {}
            for k in range(1, K):
                if abs(c[k]) > 1e-13 and abs(c[k + 1]) > 1e-13:
                    rho[k] = c[k + 1] / c[k]
            log(f"## q={q}  L={L}  (dim={n}, lambda_1={l1:.8f})")
            # gate vs exact rho (k<=L)
            gate_ok = True
            gk = [k for k in sorted(rho) if k in ref and k <= max(PLAN[q])]
            for k in gk[:4]:
                rel = abs(rho[k] - ref[k])
                if k <= L and rel >= 1e-5:
                    gate_ok = False
                log(f"   rho_{k}(iter)={rho[k]:.6f}  rho_{k}(exact)={ref[k]:.6f}  |diff|={rel:.2e}  {'k<=L' if k<=L else ''}")
            log(f"   G_MATCH(k<=L): {'PASS' if gate_ok else '*** FAIL ***'}")
            # settled tail
            tail_ks = [k for k in sorted(rho) if k >= 3][-8:]
            log(f"   rho_k tail: {['%.5f' % rho[k] for k in tail_ks]}")
            rq = rho[tail_ks[-1]] if tail_ks else float('nan')
            rq_by_L[(q, L)] = rq
            log(f"   => r_q(L={L}) ~ {rq:.6f}")
            # modal readout
            l1m, modes = modal(M, idx[(1, 1, 0)])
            log(f"   modal (mu_i=lambda_i/lambda_1, A_i=amplitude):")
            shown = 0
            for lam, A in modes:
                mu = lam / l1m
                if abs(mu) > 1.0001:
                    continue
                tag = "  <- Perron(mu=1)" if abs(mu - 1) < 1e-3 else ""
                log(f"      mu={float(np.real(mu)):+.5f}{'' if abs(np.imag(mu))<1e-9 else '%+.5fi'%np.imag(mu)}"
                    f"   |A|={abs(A):.3e}   |A*(mu-1)|={abs(A*(mu-1)):.3e}{tag}")
                shown += 1
                if shown >= 10:
                    break
            log("")

    log("## R_Q convergence in L")
    log(f"   {'q':>4} {'L=1':>10} {'L=2':>10} {'L=3':>10}")
    for q in [3, 5, 7, 11]:
        row = " ".join(f"{rq_by_L.get((q,L), float('nan')):>10.5f}" for L in [1, 2, 3])
        log(f"   {q:>4} {row}")
    log("")
    log("## READ: r_q(L) = largest |mu_i|<1 with amplitude. |A*(mu-1)| ranks a mode's")
    log("   contribution to cross's growth. r_3->1 (degenerate), r_q<1 q>=5. Watch L-conv.")
    flush()


def flush():
    with open("result_26_amplitude_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
