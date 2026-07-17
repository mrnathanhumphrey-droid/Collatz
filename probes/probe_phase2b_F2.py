"""
PROBE F2 -- gate LEMMA E-FORM (the judge) + LTE cascade + QSD braid + 2x2 effective-pair dump.
INSTRUMENT LAW: dense/direct at q=3. No rate fit; deviations reported AS deviations.

F2-1  E-FORM judge:  entry[(e,g)->(e',g')] = R(e'-e) * N(e',g,g') / D
        R(s)          = sum_da w_da w_{da-s}   (circular autocorr = Real-T1's R at k=0)
        N(e',g,g')    = #{units u : (g + u*(1-2^{e'})) %3==0  AND  ((g+u*(1-2^{e'}))//3)%q^L == g'}
        D             = 2*3^{L-1}
      Recompute EVERY entry of the frozen compressed chain (= build_compressed) from the formula;
      pre-registered EXACT match (<=1e-12) + zero pattern (entry==0 <=> N==0, since R(s)>0 always).

F2-2  LTE ladder:  v3(T) = v3(2^{e'}-1) = 1 + v3(e'/2) for even e'.  Pre-reg: at L=3 the valuation-2
      targets are EXACTLY e' in {6,12}; all other even e' are valuation 1; odd e' die from g=0.
      Theorem-shape readouts: (a) exact integer ladder, (b) from g=0 the carry valuation into target
      e' is DETERMINISTICALLY t(e')-1, (c) reproduce C2's per-class cascade deviation and attribute.

F2-3  QSD braid (admissibility condition 2): QSD-compressed partner ABOVE c0 @L=2, BELOW @L=3 (the
      crossing uniform provably missed). Read out the L=2 side. PASS -> QSD carries proof weight.

F2-4  2x2 effective-pair dump: spectrally project onto span{c0-mode, partner} (direct, EP-robust),
      report the 2x2 in the (kinematic=gamma0-aligned, tower) basis: detunings (diag), couplings
      (offdiag), trace/det/discriminant (basis-invariant). L=2,3 from the full operator. NO fit.
"""
import numpy as np
from fractions import Fraction

from probe_phase2a_q2b_q6 import build_M_gen, subgroup
from probe_phase2b_E import build_compressed

LOG = []
def log(m=""):
    try: print(m, flush=True)
    except UnicodeEncodeError: print(str(m).encode("ascii", "replace").decode("ascii"), flush=True)
    LOG.append(str(m))


# ----------------------------- shared -----------------------------
def setup_real(q, L, lam=0.5):
    qL = q ** L
    sub = subgroup(2 % qL, qL); D = len(sub)
    dl = {}; x = 1 % qL
    for e in range(D): dl[x] = e; x = (x * 2) % qL
    raw = np.array([lam ** d for d in range(1, D + 1)]); w = raw / raw.sum()
    return qL, sub, D, dl, w


def v3(n):
    if n == 0: return None
    v = 0
    while n % 3 == 0: n //= 3; v += 1
    return v


# ============================ F2-1 : E-FORM JUDGE ============================
def eform_lmat(q, L, lam=0.5):
    """Build the compressed chain PURELY from the closed-form E-FORM (R * N / D)."""
    qL, sub, D, dl, w = setup_real(q, L, lam)
    R = np.array([np.sum(w * np.roll(w, s)) for s in range(D)])         # circular autocorr
    units = np.array(sub, dtype=np.int64)                              # the D units (=<2>)
    two_e = np.array([pow(2, e, qL) for e in range(D)], dtype=np.int64)
    m = D * qL
    Lm = np.zeros((m, m))
    Ncount = {}                                                        # (e',g,g') -> N (for audit)
    for e in range(D):
        for ep in range(D):
            s = (ep - e) % D
            fac = R[s] / D
            coef = (1 - two_e[ep]) % qL
            T = (units * coef) % qL                                    # over all units u
            for g in range(qL):
                S = g + T
                mask = (S % q) == 0
                if not mask.any():
                    continue
                gp = (S // q) % qL
                gpm = gp[mask]
                np.add.at(Lm, (e * qL + g, ep * qL + gpm), fac)        # each u adds fac -> fac*N
    return Lm, R, D, qL


def f2_1(L):
    log(f"\n## F2-1 -- GATE LEMMA E-FORM (q=3 L={L}): entry = R(e'-e) * N / D")
    Lm_true, D, qL, w, _ = build_compressed(3, L)                      # the actual frozen chain
    Lm_form, R, _, _ = eform_lmat(3, L)
    m = D * qL
    diff = np.abs(Lm_true - Lm_form)
    maxd = diff.max()
    # zero pattern
    sup_true = np.abs(Lm_true) > 1e-14
    sup_form = np.abs(Lm_form) > 1e-14
    sup_ok = np.array_equal(sup_true, sup_form)
    nnz = int(sup_true.sum())
    # R(s) positivity -> zero pattern is entirely N==0
    Rmin = R.min()
    log(f"   dim={m}  nonzeros={nnz}  max|entry_true - entry_form| = {maxd:.3e}  "
        f"({'PASS <=1e-12' if maxd <= 1e-12 else 'FAIL'})")
    log(f"   support identical (zero pattern): {sup_ok}   min R(s)={Rmin:.4e} (>0 => every zero is N==0)")
    # localize any mismatch
    if maxd > 1e-12:
        bad = np.argwhere(diff > 1e-12)
        for r, c in bad[:12]:
            log(f"      MISMATCH ({r//qL},{r%qL})->({c//qL},{c%qL})  true={Lm_true[r,c]:.10g} form={Lm_form[r,c]:.10g}")
    # dump per-entry comparison at L=2 (the localized artifact)
    if L == 2:
        lines = [f"# E-FORM gate q=3 L=2: every nonzero entry, actual (build_compressed) vs formula R*N/D",
                 "# src(e,g)\tdst(e',g')\tactual\tformula\ts=e'-e\tR(s)\tN=actual*D/R(s)"]
        for r in range(m):
            for c in range(m):
                if sup_true[r, c] or sup_form[r, c]:
                    e, g = r // qL, r % qL; ep, gp = c // qL, c % qL
                    s = (ep - e) % D
                    Nest = Lm_true[r, c] * D / R[s] if R[s] > 0 else float('nan')
                    lines.append(f"({e},{g})\t({ep},{gp})\t{Fraction(Lm_true[r,c]).limit_denominator(10**7)}"
                                 f"\t{Fraction(Lm_form[r,c]).limit_denominator(10**7)}\t{s}"
                                 f"\t{Fraction(R[s]).limit_denominator(10**7)}\t{round(Nest)}")
        with open("outputs/eform_gate_q3_L2.tsv", "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        log("   dumped outputs/eform_gate_q3_L2.tsv (every nonzero entry, actual vs formula, with s,R,N)")
    return maxd <= 1e-12 and sup_ok


# ============================ F2-2 : LTE CASCADE ============================
def f2_2(L):
    log(f"\n## F2-2 -- LTE ladder  t(e')=v3(2^{{e'}}-1)=1+v3(e'/2) (q=3 L={L})")
    qL, sub, D, dl, w = setup_real(3, L)
    # (a) exact integer ladder
    ladder = {ep: v3(2 ** ep - 1) for ep in range(1, D)}
    two_targets = sorted([ep for ep, t in ladder.items() if t == 2])
    one_targets = sorted([ep for ep, t in ladder.items() if t == 1])
    zero_targets = sorted([ep for ep, t in ladder.items() if t == 0])
    log(f"   e'=0 -> T=0 (self-phase, no shift).  valuation-2 targets e' = {two_targets}")
    log(f"   valuation-1 (even) e' = {one_targets}")
    log(f"   valuation-0 (odd, die from g=0) e' = {zero_targets}")
    if L == 3:
        pred = (two_targets == [6, 12] and
                one_targets == [2, 4, 8, 10, 14, 16] and
                zero_targets == [1, 3, 5, 7, 9, 11, 13, 15, 17])
        log(f"   PRE-REG (L=3): valuation-2 == {{6,12}} AND even->1 AND odd->0 :  {'CONFIRMED' if pred else 'DEVIATION'}")
    # (b) deterministic g=0 carry law: from g=0, carry valuation into target e' == t(e')-1
    units = np.array(sub, dtype=np.int64)
    two_e = np.array([pow(2, e, qL) for e in range(D)], dtype=np.int64)
    ok_det = True; examples = []
    for ep in range(1, D):
        coef = (1 - two_e[ep]) % qL
        T = (units * coef) % qL
        gate = (T % 3) == 0
        if not gate.any():
            examples.append((ep, "odd: no g=0 transition", ladder[ep]))
            if ladder[ep] != 0: ok_det = False
            continue
        carry = (T[gate] // 3) % qL
        cv = np.array([v3(int(c)) if c != 0 else L for c in carry])     # v3(carry); 0 -> L (max)
        # predicted: t(e')-1 (uniform); carry==0 has valuation >= L (only when T divisible high)
        predv = ladder[ep] - 1
        uniform = np.all(cv == cv[0])
        match = (cv[0] == predv) if predv < L else (cv[0] >= L)
        if not (uniform and match): ok_det = False
        examples.append((ep, f"t={ladder[ep]} -> carry v3={cv[0]} (pred {predv}), uniform={uniform}", ladder[ep]))
    log(f"   from g=0: carry valuation into target e' is DETERMINISTIC = t(e')-1 :  "
        f"{'CONFIRMED (all even targets)' if ok_det else 'DEVIATION (see examples)'}")
    for ep, msg, t in examples[:8]:
        log(f"      e'={ep:2d}: {msg}")
    # (c) reproduce C2's per-(theta=e mod3, gamma) cascade deviation + attribute to LTE
    M, idx, n = build_M_gen(3, L, 2, [0.5 ** d for d in range(1, D + 1)])
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    erho = np.array([dl[(b * pow(a, -1, qL)) % qL] for (a, b, g) in states])
    gam = np.array([s[2] for s in states])
    cls_of = list(zip((erho % 3).tolist(), gam.tolist()))
    Mc = M.tocoo()
    from collections import defaultdict
    dist = defaultdict(lambda: np.zeros(L + 1)); tot = defaultdict(float)
    for r, c, val in zip(Mc.row, Mc.col, Mc.data):
        dist[cls_of[c]][v3(states[r][2]) if states[r][2] != 0 else L] += val
        tot[cls_of[c]] += val
    law = np.array([2 * 3 ** (-(j + 1)) for j in range(L)] + [0.0]); law[L] = 1 - law[:L].sum()
    devs = np.array([np.max(np.abs(dist[sc] / tot[sc] - law)) for sc in sorted(dist)])
    log(f"   C2 reproduction: per-class max deviation from uniform cascade  min={devs.min():.2e} "
        f"med={np.median(devs):.2e} max={devs.max():.2e}  (class-DEPENDENT; the target-side driver is the LTE ladder)")
    if L == 3:
        return None  # ladder object dumped once (below at L=3 call site)
    return ladder


# ============================ F2-3 : QSD BRAID ============================
def qsd_partner_side(L, partner_val):
    """Q-weighted (quasi-stationary) compressed partner and its side vs c0."""
    qL, sub, D, dl, w = setup_real(3, L)
    M, idx, n = build_M_gen(3, L, 2, [0.5 ** d for d in range(1, D + 1)])
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    erho = np.array([dl[(b * pow(a, -1, qL)) % qL] for (a, b, g) in states])
    gam = np.array([s[2] for s in states])
    cls = list(zip(erho.tolist(), gam.tolist()))
    classes = sorted(set(cls)); ci = {c: j for j, c in enumerate(classes)}; m = len(classes)
    cls_idx = np.array([ci[c] for c in cls])
    Md = M.toarray()
    evR, VR = np.linalg.eig(Md)
    muQ = np.abs(VR[:, int(np.argmax(np.abs(evR)))])
    # weighted lump
    num = np.zeros((m, m)); den = np.zeros(m)
    Mc = M.tocoo()
    np.add.at(num, (cls_idx[Mc.col], cls_idx[Mc.row]), muQ[Mc.col] * Mc.data)
    np.add.at(den, cls_idx, muQ); den[den == 0] = 1.0
    LQ = num / den[:, None]
    c0 = float(np.sum(w ** 2))
    cks = np.array([complex(np.sum(w ** 2 * np.exp(2j * np.pi * k * (np.arange(D) + 1) / D))) for k in range(D)])
    ev = np.linalg.eigvals(LQ)
    i0 = int(np.argmin(np.abs(ev - c0)))
    cand = [i for i in range(len(ev)) if i != i0 and min(abs(ev[i] - cks)) > 1e-9]
    ip = min(cand, key=lambda i: abs(ev[i] - partner_val))
    p = ev[ip].real
    return p, c0, ("above" if p > c0 else "below"), partner_val, ("above" if partner_val > c0 else "below")


def f2_3():
    log("\n## F2-3 -- QSD ADMISSIBILITY (the braid): QSD partner ABOVE c0 @L=2, BELOW @L=3?")
    p2, c02, s2, t2, ts2 = qsd_partner_side(2, 0.346827)
    p3, c03, s3, t3, ts3 = qsd_partner_side(3, 0.333236)
    log(f"   L=2: c0={c02:.6f}  QSD partner={p2:.6f} ({s2})   [true partner {t2} {ts2}]")
    log(f"   L=3: c0={c03:.6f}  QSD partner={p3:.6f} ({s3})   [true partner {t3} {ts3}]")
    braid = (s2 == "above" and s3 == "below")
    matches_true = (s2 == ts2 and s3 == ts3)
    log(f"   BRAID (above@L2, below@L3): {'CONFIRMED' if braid else 'NOT SEEN'}   "
        f"matches TRUE partner's braid: {matches_true}")
    log(f"   => {'QSD carries proof weight (reproduces the crossing uniform missed).' if (braid and matches_true) else 'QSD demoted to exploratory; proceed on projection (F2-4).'}")
    return braid and matches_true


# ============================ F2-4 : 2x2 EFFECTIVE PAIR ============================
def eff_2x2(L, partner_val):
    qL, sub, D, dl, w = setup_real(3, L)
    M, idx, n = build_M_gen(3, L, 2, [0.5 ** d for d in range(1, D + 1)])
    states = [None] * n
    for s, i in idx.items(): states[i] = s
    gam = np.array([s[2] for s in states])
    Md = M.toarray()
    c0 = float(np.sum(w ** 2))
    evR, VR = np.linalg.eig(Md)
    i0 = int(np.argmin(np.abs(evR - c0)))
    ip = int(np.argmin(np.abs(evR - partner_val)))
    lam0, lamp = evR[i0], evR[ip]
    # orthonormal basis of the invariant subspace span{r0, rp}
    U = np.column_stack([VR[:, i0], VR[:, ip]])
    Q, _ = np.linalg.qr(U)
    # orient to (kinematic = gamma=0 aligned, tower) via gamma=0 projector restricted to subspace
    p0 = (gam == 0).astype(float)
    G0 = Q.conj().T @ (p0[:, None] * Q)                                # 2x2 Hermitian
    wv, Vv = np.linalg.eigh(G0)                                        # ascending
    kin = Vv[:, 1]; tow = Vv[:, 0]                                     # max gamma0-weight = kinematic
    Qp = Q @ np.column_stack([kin, tow])
    B = Qp.conj().T @ (Md @ Qp)                                       # 2x2 effective (kin,tow) basis
    tr = np.trace(B); det = np.linalg.det(B)
    disc = tr ** 2 - 4 * det                                          # = (lam0 - lamp)^2
    gap = lam0 - lamp
    kin_wt = float(wv[1]); tow_wt = float(wv[0])                       # gamma=0 weight of each basis vec
    return dict(B=B, c0=c0, lam0=lam0.real, lamp=lamp.real, tr=tr, det=det,
                disc=disc, gap=gap.real, kin_wt=kin_wt, tow_wt=tow_wt)


def f2_4():
    log("\n## F2-4 -- 2x2 EFFECTIVE-PAIR (basis: kinematic=gamma0-aligned, tower); direct/EP-robust")
    rows = []
    for L, pv in [(2, 0.346827), (3, 0.333236)]:
        d = eff_2x2(L, pv)
        B = d['B']
        log(f"\n   L={L}: c0(kinematic anchor)={d['c0']:.8f}  partner={d['lamp']:.8f}  gap(c0-partner)={d['gap']:.3e}")
        log(f"      B[kin,kin]={B[0,0].real:+.8f}{B[0,0].imag:+.2e}j   B[kin,tow]={B[0,1].real:+.8f}{B[0,1].imag:+.2e}j")
        log(f"      B[tow,kin]={B[1,0].real:+.8f}{B[1,0].imag:+.2e}j   B[tow,tow]={B[1,1].real:+.8f}{B[1,1].imag:+.2e}j")
        log(f"      basis-invariants: tr={d['tr'].real:+.8f}  det={d['det'].real:+.8f}  "
            f"discriminant=(l0-lp)^2={d['disc'].real:.6e}  sqrt(disc)={np.sqrt(abs(d['disc'].real)):.6e}")
        log(f"      gamma=0 weight of basis vecs: kinematic={d['kin_wt']:.4f}  tower={d['tow_wt']:.4f} "
            f"(kinematic localizes on gamma=0; tower off it -- P-consistent)")
        rows.append((L, d))
    # dump
    lines = ["# 2x2 effective-pair matrix, q=3, basis (kinematic=gamma0-aligned, tower). Coalescence-derivation JUDGE.",
             "# L\tc0\tpartner\tB_kk\tB_kt\tB_tk\tB_tt\ttrace\tdet\tdiscriminant=(l0-lp)^2\tsqrt_disc\tkin_gamma0wt\ttow_gamma0wt"]
    for L, d in rows:
        B = d['B']
        lines.append("\t".join(str(x) for x in [
            L, f"{d['c0']:.10f}", f"{d['lamp']:.10f}",
            f"{B[0,0].real:.10f}", f"{B[0,1].real:.10f}", f"{B[1,0].real:.10f}", f"{B[1,1].real:.10f}",
            f"{d['tr'].real:.10f}", f"{d['det'].real:.10f}", f"{d['disc'].real:.6e}",
            f"{np.sqrt(abs(d['disc'].real)):.6e}", f"{d['kin_wt']:.6f}", f"{d['tow_wt']:.6f}"]))
    with open("outputs/effective_2x2_q3.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log("\n   dumped outputs/effective_2x2_q3.tsv (L=2,3).  L=4 SIZED+DEFERRED (see below).")
    log("   L=4 note: the QSD-weighted L=4 2x2 needs mu=|dominant right eigvec| of the 236,196-state")
    log("             operator (nnz ~2.3e8, ~5.5GB CSR). mu = LARGEST mode -> power iteration (pure sparse")
    log("             matvec, GPU-friendly, NOT the banned interior LU), BUT near-degenerate top (gap<1e-4)")
    log("             => slow subspace separation; heavy compute -> Lambda/greenlight. CONFIRMATORY only")
    log("             (rate-fit banned); L=2,3 are the load-bearing points for Wilson's L-trend derivation.")


def emit_lte_dump(L=3):
    qL, sub, D, dl, w = setup_real(3, L)
    lines = [f"# LTE ladder q=3 L={L}: t(e') = v3(2^e' - 1) = 1 + v3(e'/2) for even e'. Governs carry valuation.",
             "# e'\tt(e')=v3(2^e'-1)\t2^e'-1\tnote"]
    for ep in range(D):
        val = 2 ** ep - 1
        t = v3(val)
        note = "self (T=0)" if ep == 0 else ("odd: dies from g=0" if t == 0 else f"even: carry v3 = t-1 = {t-1} from g=0")
        lines.append(f"{ep}\t{t if t is not None else 'inf'}\t{val}\t{note}")
    with open("outputs/lte_ladder_q3.tsv", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    log("# PROBE F2 -- E-FORM judge + LTE cascade + QSD braid + 2x2 effective-pair. Dense/direct at q=3.")
    g1_2 = f2_1(2); g1_3 = f2_1(3)
    log(f"\n   >> F2-1 VERDICT: L=2 {'PASS' if g1_2 else 'FAIL'}, L=3 {'PASS' if g1_3 else 'FAIL'}  "
        f"=> LEMMA E-FORM {'GATE-CONFIRMED' if (g1_2 and g1_3) else 'NOT confirmed'}")
    f2_2(2); f2_2(3); emit_lte_dump(3)
    braid = f2_3()
    f2_4()
    log(f"\n## SUMMARY: E-FORM {'confirmed' if (g1_2 and g1_3) else 'FAILED'}; QSD braid {'PASS' if braid else 'FAIL'}; "
        f"2x2 dumped (L=2,3), L=4 deferred.")
    with open("logs/probe_phase2b_F2_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
