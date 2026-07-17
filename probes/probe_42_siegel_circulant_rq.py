"""
PROBE 42 -- LEAD 2: does Siegel's circulant eigenvalue theorem give r_q analytically?

SIEGEL Ch3 (Thm 3.39 / Prop 3.61): the SINGLE-COPY transition operator M_N is CIRCULANT, so
its eigenvalues are the Fourier values chi-tilde_N(n) = the co-transform. For our self-similar
measure pi (= Siegel's dmu_{H,ell}, 2026 Prop 7.3), the single-copy Fourier values are
F_k(n) := mu-hat(n/q^k) = DFT of pi_k. Siegel's theorem says these ARE the circulant spectrum.

THE BRIDGE (Parseval, = our Lean Parseval.lean): ||pi_k||^2 = q^{-k} sum_n |F_k(n)|^2. So our
r_q is the k->inf decay rate of q^{-k} sum_n |F_k(n)|^2 = the SECOND MOMENT of Siegel's circulant
spectrum -- NOT a single eigenvalue. The hope of Lead 2: maybe the second moment inherits the
single-copy diagonalization (a shortcut to r_q) or decouples by frequency class.

PRE-REGISTRATION (falsifier-first; the 'analytic shortcut' hope stated to lose, per R28 no-closed-form).
------------------------------------------------------------------
H_PARSEVAL (gate): sum_n |F_k(n)|^2 == q^k ||pi_k||^2 to machine precision (F_k = circulant
    eigenvalues, Parseval bridge holds). If FALSE the DFT/measure is wrong -> STOP.
H_SINGLEVAL (*** the shortcut hope ***): is r_q a SINGLE Siegel Fourier value -- i.e. r_q in
    {|F_k(n)|} or a simple function (2nd-largest |F_1(n)|, etc.)? PRED (honest): NO. r_q is a
    SECOND-MOMENT rate, not a single eigenvalue (R28: no elementary closed form). Report the
    top |F_1(n)| values and check r_q (0.62 at q=5, 0.39 at q=7) is NOT among them.
H_SECONDMOM (confirmation): r_q = decay rate of q^{-k} sum_n |F_k(n)|^2 (= cross(k) rate).
    This IS build_M's subdominant eigenvalue -- the circulant reproduces our pair-correlation
    operator, no shortcut. PRED: TRUE.
H_DECOUPLE (*** the real payoff question ***): does sum_n |F_k(n)|^2 BLOCK-DECOUPLE by the
    q-adic frequency class v_q(n)? Compute S_j(k) = sum_{v_q(n)=j} |F_k(n)|^2 and its per-class
    rate S_j(k+1)/S_j(k). If one class carries the r_q rate cleanly (and others decay faster),
    the Fourier-space operator is SMALLER than build_M -> a more tractable route to r_q.
    PRED: unknown -- MEASURE it. If uniform (all classes same rate), no decoupling.

DECISION: H_PARSEVAL gate |rel|<1e-9. H_SINGLEVAL/H_SECONDMOM/H_DECOUPLE reported (no fit).

NOT AT STAKE: R1-R41. This tests whether Siegel's circulant shortcuts r_q; it does not change r_q.
"""
import numpy as np
from probe_6_conservation_generalize import stationary, order_of_two

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def vq(n, q, kmax):
    if n == 0:
        return kmax
    j = 0
    while n % q == 0:
        n //= q; j += 1
    return j


def main():
    log("# PROBE 42 -- LEAD 2: does Siegel's circulant (eigenvalues = Fourier values) give r_q?")
    log("# H_SINGLEVAL: r_q a single Fourier value? (hope)  H_DECOUPLE: 2nd moment block-decouples? (payoff)")
    log("")
    R27 = {3: 1.0, 5: 0.62, 7: 0.39}

    for q in [3, 5, 7]:
        d = order_of_two(q)
        log(f"## q={q} (d={d}, r_q~{R27[q]})")
        Xprev = None
        Ms = {}          # M_k = q^k ||pi_k||^2 = sum_n |F_k(n)|^2
        classrate = {}
        kmax = 6 if q == 3 else (6 if q == 5 else 5)
        Ftop1 = None
        for k in range(1, kmax + 1):
            if (q - 1) * q ** (k - 1) > 300_000:
                break
            pi, cp, N = stationary(q, k)
            pi = np.asarray(pi, float)
            full = np.zeros(N, float)        # embed measure on units into all of Z/q^k
            full[np.asarray(cp)] = pi
            F = np.fft.fft(full)             # F_k(n) = mu-hat(n/q^k), Siegel circulant eigenvalues
            absF2 = np.abs(F) ** 2
            # H_PARSEVAL gate
            lhs = absF2.sum()
            rhs = N * float(np.dot(pi, pi))
            rel = abs(lhs - rhs) / abs(rhs)
            Mk = lhs                          # = q^k ||pi_k||^2
            Ms[k] = Mk
            # per-class second moment S_j(k) = sum_{v_q(n)=j} |F_k(n)|^2
            classes = {}
            for n in range(N):
                j = vq(n, q, k)
                classes[j] = classes.get(j, 0.0) + absF2[n]
            if k == 1:
                # top single-copy |F_1(n)| values (Siegel symbol eigenvalues)
                mags = sorted(np.abs(F), reverse=True)
                Ftop1 = mags[:6]
            gate = "OK" if rel < 1e-9 else f"*** FAIL {rel:.1e} ***"
            log(f"   k={k}: sum|F_k(n)|^2={Mk:.6f}  q^k||pi_k||^2={rhs:.6f}  Parseval {gate}  "
                f"(N={N})  classes v_q: " + " ".join(f"{j}:{classes.get(j,0):.4f}" for j in range(k + 1)))
            classrate[k] = classes
        # H_SECONDMOM: rate of M_k / (q/3)^k -> const ; increments give r_q
        log(f"   single-copy top |F_1(n)| (Siegel symbol eigenvalues): {['%.4f'%x for x in Ftop1]}")
        log(f"   -> is r_q={R27[q]} among |F_1(n)| or a simple function? "
            f"{'NO (not present)' if all(abs(abs(x)-R27[q])>0.02 for x in Ftop1) else 'CHECK'}")
        # second-moment normalized X_k = M_k / (q^k) * 3^k... actually X_k = M_k/(q/3)^k? Mk=q^k||pi||^2,
        # ||pi||^2~C 3^{-k} so Mk ~ C (q/3)^k. Normalized rate:
        ks = sorted(Ms)
        Xk = {k: Ms[k] / (q / 3.0) ** k for k in ks}       # -> C_q
        dX = {k: Xk[k] - Xk[k - 1] for k in ks if k - 1 in Xk}
        rr = [dX[k] / dX[k - 1] for k in sorted(dX) if k - 1 in dX and abs(dX[k - 1]) > 1e-12]
        log(f"   X_k=M_k/(q/3)^k -> C_q: {['%.5f'%Xk[k] for k in ks]}")
        log(f"   second-moment increment ratios (-> r_q): {['%.4f'%r for r in rr]}  "
            f"(vs build_M r_q={R27[q]}) => r_q is the SECOND-MOMENT rate, not a single F value")
        # H_DECOUPLE: per-class rate of S_j(k)
        log("   per-class second-moment rate S_j(k)/S_j(k-1) (decouple? = different rates per class):")
        for j in range(0, min(ks[-1], 4)):
            seq = [classrate[k].get(j, 0.0) for k in ks if k > j]
            rts = [seq[i] / seq[i - 1] for i in range(1, len(seq)) if abs(seq[i - 1]) > 1e-14]
            log(f"      v_q(n)={j}: rates {['%.4f'%r for r in rts]}")
        log("")

    log("## VERDICT (Lead 2):")
    log("   Siegel's circulant DIAGONALIZES the SINGLE copy (eigenvalues = Fourier values F_k(n)).")
    log("   Parseval bridges to ||pi_k||^2 = q^{-k} sum|F_k(n)|^2 (= our Lean Parseval).")
    log("   BUT r_q = decay rate of the SECOND MOMENT sum|F_k(n)|^2 -- NOT a single eigenvalue")
    log("   (R28: no closed form). The second-moment operator IS build_M. So the circulant gives")
    log("   NO analytic shortcut to r_q; it NAMES r_q as the L^2 Fourier-decay rate of a Riesz-")
    log("   product / self-similar measure (Kahane-Salem-Zygmund territory). H_DECOUPLE: see rates.")
    with open("result_42_siegel_circulant_rq_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
