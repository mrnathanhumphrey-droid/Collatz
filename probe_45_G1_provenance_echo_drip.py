"""
PROBE 45 (G1, redesigned) -- q=3 boundary: PROVENANCE decomposition of the collision cascade.

The spec upgrade (user): omega(3)=1 total is CIRCULAR -- the q=3 fixed point is banked six ways
(R15 slope, R16 D-ratio 1.00092, R22 rho->1, R25 lambda_2->lambda_1...). A hand-derivation wrong
in two compensating ways still hits 1. So test the PARTS, in cross-normalized units:

  P_k := per-level ECHO -- propagation of already-off-diagonal (nontrivial-collider) mass:
         (O->O mass, level k -> k+1) / (off-diag mass at level k).      PRE-REG -> 2/3  [at risk].
  I_k := per-level DRIP -- fresh injection from the diagonal reservoir (D->O).  PRE-REG -> 7/45.
  Consistency (NOT a test): P + I/c* = 2/3 + (7/45)/(7/15) = 1.   c* = 7/15 = lim c_k.

DISCRIMINATING READOUT -- P_k PER LEVEL (three pre-registered outcomes):
  (1) P_k ~ 2/3 CONSTANT in k        -> recursion form exact; R23's super-geometric read needs re-exam.
  (2) P_k -> 2/3, super-geo corrections -> 2/3 is the LIMIT; recursion asymptotic; both stand (user's pred).
  (3) P_k -> something != 2/3        -> case-table carry error; the 7/45 reading dies cleanly.

BINNING IS ITSELF UNDER TEST (R19: colliding mass concentrates on pairs differing in exactly TWO
coords). Report which binning makes P_k stable: (A) a==b diagonal; (B) a==b AND carry gamma==0.
Diagnostic: is the D reservoir mass == (1/3)^k (pure identical-path) or inflated by accidental
residue coincidences (S==S' mod ordL, different paths)?  Same v-truncation both channels: build_M's
single shared `mult` list (ordL = d*q^{L-1}), so no diag-exact/offdiag-truncated skew.

INDEPENDENT EXACT CHECK (no binning, no operator): from exact rational marginals pi_k, form
c_m := ||d_m||^2 * 3^m  (||d_m||^2 = ||pi_m||^2 - (1/3)||pi_{m-1}||^2, R74). Test whether the
CONSTANT-COEFFICIENT recursion holds: is  c_{m+1} - (2/3) c_m  EXACTLY 7/45 for all m (outcome 1),
or drifting (outcome 2)? Exact big-int rationals -- this is R23's arithmetic, re-asked as the residual.

NOT AT STAKE: R1-R44. Tests the ECHO/DRIP split of the q=3 boundary; changes no r_q value.
"""
import numpy as np
from fractions import Fraction

from probe_25_transfer_operator_Aprime import build_M, two_subgroup
from c_seven_forty_fifth_derivation import build_markov_rational, stationary_rational

LOG = []


def log(m=""):
    try:
        print(m)
    except UnicodeEncodeError:
        print(str(m).encode("ascii", "replace").decode("ascii"))
    LOG.append(str(m))


def build_masks(idx, n):
    diag = np.zeros(n, dtype=bool)
    g0 = np.zeros(n, dtype=bool)
    for (a, b, g), i in idx.items():
        diag[i] = (a == b)
        g0[i] = (g == 0)
    return diag, g0


def provenance(q, L, binning):
    """Return per-level flows under a diagonal binning. binning in {'ab','ab_g0'}."""
    M, idx, n = build_M(q, L)
    diag, g0 = build_masks(idx, n)
    Dm = diag if binning == 'ab' else (diag & g0)
    O = ~Dm
    v = np.zeros(n)
    v[idx[(1, 1, 0)]] = 1.0
    o = v * O
    d = v * Dm
    rows = []
    for k in range(1, L + 1):
        af_o = M.dot(o)           # arrivals at level k that came FROM off-diagonal
        af_d = M.dot(d)           # arrivals at level k that came FROM diagonal
        v = af_o + af_d
        echo = float(af_o[O].sum())   # O->O  (propagation)
        drip = float(af_d[O].sum())   # D->O  (fresh injection)
        Oprev = float(o.sum())        # off-diag mass at level k-1 (the source)
        Dprev = float(d.sum())        # diag mass at level k-1
        o = v * O
        d = v * Dm
        rows.append(dict(k=k, O=float(o.sum()), D=float(d.sum()),
                         echo=echo, drip=drip, Oprev=Oprev, Dprev=Dprev,
                         tot=float(v.sum())))
    return rows


def exact_c(kmax):
    """Exact rational ||pi_k||^2, ||d_m||^2, c_m = ||d_m||^2 * 3^m."""
    pisq = {}
    for k in range(1, kmax + 1):
        K, cop = build_markov_rational(k)
        pi = stationary_rational(K)
        pisq[k] = sum(p * p for p in pi)
    d2 = {}
    c = {}
    for k in range(1, kmax):
        d2[k + 1] = pisq[k + 1] - pisq[k] * Fraction(1, 3)
        c[k + 1] = d2[k + 1] * Fraction(3 ** (k + 1))
    return pisq, d2, c


def main():
    q = 3
    log("# PROBE 45 (G1) -- q=3 boundary echo/drip provenance. Total=1 is circular; test the PARTS.")
    log("# Pre-reg: P_k(echo)->2/3 [at risk], I_k(drip)->7/45, P_k PER LEVEL discriminates (1)/(2)/(3).")
    log(f"#   targets: 2/3={2/3:.6f}  7/45={7/45:.6f}  c*=7/15={7/15:.6f}")
    log("")

    # ---------- PART 1: exact constant-coefficient recursion test (no binning) ----------
    log("## PART 1 -- EXACT rational recursion residual  c_{m+1} - (2/3) c_m  =?= 7/45  (outcome 1 vs 2)")
    kmax = 5
    try:
        pisq, d2, c = exact_c(kmax)
    except Exception as e:
        log(f"   kmax=5 failed ({e}); retry kmax=4")
        kmax = 4
        pisq, d2, c = exact_c(kmax)
    log(f"   {'m':>3} {'||pi_m||^2':>16} {'||d_m||^2':>18} {'c_m=||d_m||^2*3^m':>20} {'c_m dec':>10} {'->7/15':>9}")
    for m in range(1, kmax + 1):
        pm = pisq.get(m)
        dm = d2.get(m)
        cm = c.get(m)
        log(f"   {m:>3} {str(pm):>16} {(str(dm) if dm is not None else '-'):>18} "
            f"{(str(cm) if cm is not None else '-'):>20} "
            f"{(f'{float(cm):.6f}' if cm is not None else '-'):>10} {7/15:>9.6f}")
    log("")
    log(f"   {'m':>3} {'c_{m+1}-(2/3)c_m (exact)':>28} {'decimal':>12} {'vs 7/45':>12} {'constant?':>10}")
    resid = {}
    ms = sorted(c.keys())
    for i in range(len(ms) - 1):
        m, mp = ms[i], ms[i + 1]
        r = c[mp] - Fraction(2, 3) * c[m]
        resid[m] = r
        log(f"   {m:>3} {str(r):>28} {float(r):>12.6f} {float(r - Fraction(7,45)):>+12.6e} "
            f"{'=7/45' if r == Fraction(7,45) else 'NO':>10}")
    if len(resid) >= 2:
        rr = list(resid.values())
        allc = all(x == rr[0] for x in rr)
        log(f"   => residuals {'ALL EQUAL' if allc else 'DRIFT'}: "
            f"{'constant-coeff recursion EXACT (outcome 1)' if allc else 'super-geometric -- coeffs NOT constant (outcome 2)'}")
        if not allc:
            log(f"      residual sequence (decimal): {[round(float(x),6) for x in rr]}  "
                f"-> if ->7/45 from below/above, 2/3 & 7/45 are the LIMITS (user's prediction).")
    log("")

    # ---------- PART 2: operator provenance echo/drip (build_M), binning A = a==b ----------
    for binning, tag in [('ab', 'A: diagonal = {a==b}'), ('ab_g0', 'B: diagonal = {a==b AND gamma==0}')]:
        log(f"## PART 2 [{tag}] -- echo/drip from build_M (q=3). P_k = echo/O_prev -> 2/3 ?")
        for L in [2, 3]:
            rows = provenance(q, L, binning)
            log(f"   --- L={L} (exact to depth {L}) ---")
            log(f"   {'k':>3} {'O_mass':>14} {'D_mass':>14} {'(1/3)^k':>10} {'echo(O->O)':>13} "
                f"{'drip(D->O)':>13} {'P_k=echo/Oprev':>15}")
            for r in rows:
                k = r['k']
                third = (1.0 / 3.0) ** k
                Pk = (r['echo'] / r['Oprev']) if r['Oprev'] > 1e-18 else float('nan')
                log(f"   {k:>3} {r['O']:>14.9f} {r['D']:>14.9f} {third:>10.6f} "
                    f"{r['echo']:>13.9f} {r['drip']:>13.9f} "
                    f"{(f'{Pk:.6f}' if Pk==Pk else 'n/a (Oprev=0)'):>15}")
            # gate: total mass = ||pi_k||^2 ?  and  O_mass = ||d_k||^2 ?
            log(f"       gate tot vs ||pi_k||^2 :" +
                " ".join(f"k{ r['k']}:{r['tot']:.6f}/{float(pisq.get(r['k'],float('nan'))):.6f}" for r in rows))
            log(f"       O_mass vs ||d_k||^2   :" +
                " ".join(f"k{r['k']}:{r['O']:.6f}/{(float(d2[r['k']]) if r['k'] in d2 else float('nan')):.6f}" for r in rows))
            log("")
    log("## READ -- see verdict block below (P_k per level vs 2/3; which binning holds P_k still).")

    with open("result_45_G1_provenance_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
