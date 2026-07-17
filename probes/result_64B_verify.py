"""
result_64B_verify.py -- reconstruct + verify R64.B (the load-bearing 1:4 class-mass
ratio for c=7/45), which is cited in THEOREM_C_745 + D3_DERIVATION_AUDIT + R77 §1 but
has NO result file. Closes the provenance gap flagged by the 2026-07-14 thread audit.

Class definition (result_77_T_lead_spectrum.md line 14, R66 chain rule):
    class + = v EVEN,  class - = v ODD   (v = 2-adic valuation of the halving step)
The 1:4 ratio has two equivalent elementary origins, both parity splits of Geom(1/2):
  (a) class PROBABILITIES:  P(v even)=Σ_{v even>=1}2^{-v}=1/3,  P(v odd)=2/3
      -> squared-class-mass ratio (1/3)^2:(2/3)^2 = 1:4  (the (1,4) eigenvector, R77 l.119)
  (b) Plancherel |mu_hat|^2 weights: each step contributes 2^{-v} to mu_hat, 4^{-v} to |.|^2:
      Sum_{v even>=1}4^{-v}=1/15,  Sum_{v odd>=1}4^{-v}=4/15  -> 1:4  (R77 l.19-23, T_diag prefactors)
Both give 1:4 exactly. R64.B is therefore ELEMENTARY (parity split of Geom(1/2)); the gap
was provenance-only. This script proves it exactly (Fractions) + confirms numerically on
real Collatz v-values, and checks T_diag's (1,4) eigenstructure.
"""
import sys
from fractions import Fraction as F
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

def parity_split(base_inv, Vmax=4000):
    """Sum_{v>=1} base_inv^v split by parity of v, as exact Fractions (geometric, closed)."""
    # closed forms: Sum_{v even>=2} x^v = x^2/(1-x^2); Sum_{v odd>=1} x^v = x/(1-x^2)
    x = base_inv
    even = x * x / (1 - x * x)
    odd = x / (1 - x * x)
    return even, odd

def main():
    log("# R64.B verification -- the 1:4 class-mass ratio, reconstructed")
    log("# class + = v even, class - = v odd (R77 l.14); v ~ Geom(1/2), P(v=k)=2^{-k}, k>=1")
    log("")

    # (a) class probabilities: parity split of Sum 2^{-v}
    e2, o2 = parity_split(F(1, 2))
    log("## (a) class probabilities = parity split of Geom(1/2)  (Sum 2^{-v})")
    log(f"   P(v even) = {e2} = {float(e2):.6f}   (target 1/3)")
    log(f"   P(v odd)  = {o2} = {float(o2):.6f}   (target 2/3)")
    log(f"   ratio P(even):P(odd) = {e2}:{o2} = 1:{o2/e2}")
    sq_ratio = (o2 * o2) / (e2 * e2)
    log(f"   SQUARED-class-mass ratio (1/3)^2:(2/3)^2 = 1:{sq_ratio}  <- the (1,4) eigenvector")
    log(f"   exact check: P(even)=1/3? {e2 == F(1,3)}   P(odd)=2/3? {o2 == F(2,3)}   "
        f"squared ratio=4? {sq_ratio == 4}")
    log("")

    # (b) Plancherel weight split: parity split of Sum 4^{-v}
    e4, o4 = parity_split(F(1, 4))
    log("## (b) Plancherel |mu_hat|^2 weight = parity split of Sum 4^{-v}  (R77 l.19-23)")
    log(f"   Sum_{{v even}} 4^{{-v}} = {e4}   (target 1/15)")
    log(f"   Sum_{{v odd}}  4^{{-v}} = {o4}   (target 4/15)")
    log(f"   ratio = {e4}:{o4} = 1:{o4/e4}")
    log(f"   exact check: 1/15? {e4 == F(1,15)}   4/15? {o4 == F(4,15)}   ratio=4? {o4/e4 == 4}")
    log(f"   T_diag prefactors: P^++ coeff = 3*(1/15)=1/5={3*e4}, P^-- coeff = 3*(4/15)=4/5={3*o4}")
    log("")

    # T_diag eigenstructure
    log("## T_diag = (1/5)[[1,1],[4,4]] eigenstructure")
    T = np.array([[1, 1], [4, 4]], float) / 5
    w, V = np.linalg.eig(T)
    order = np.argsort(-w)
    log(f"   eigenvalues = {sorted(w.round(10), reverse=True)}  (target {{1, 0}})")
    v1 = V[:, order[0]]; v1 = v1 / v1[0]
    log(f"   lambda=1 eigenvector (normalized) = ({v1[0]:.4f}, {v1[1]:.4f})  (target (1, 4))")
    log(f"   check: eigenvalue 1? {abs(sorted(w,reverse=True)[0]-1)<1e-12}  "
        f"eigenvector (1,4)? {abs(v1[1]-4)<1e-9}")
    log("")

    # (c) numeric confirmation on REAL Collatz v-values (invariant marginal, fresh uniform sample)
    log("## (c) numeric confirmation on real Collatz v = v2(3n+1), n uniform odd 3-nmid")
    t = np.arange(20_000_000, dtype=np.int64)
    cnt = np.zeros(64, dtype=np.int64)
    for base in (1, 5):
        n = 6 * t + base
        x = 3 * n + 1
        v = np.log2((x & (-x)).astype(np.float64)).round().astype(np.int64)
        v = np.minimum(v, 63)
        cnt += np.bincount(v, minlength=64)[:64]
    N = cnt.sum(); vs = np.arange(64)
    p_even = cnt[vs % 2 == 0].sum() / N
    p_odd = cnt[vs % 2 == 1].sum() / N
    # squared-class-mass ratio (the (1,4) eigenvector): (p_even)^2 : (p_odd)^2
    sq_num = (p_odd ** 2) / (p_even ** 2)
    log(f"   N={N:,}  empirical P(v even)={p_even:.5f} (1/3=0.33333), P(v odd)={p_odd:.5f} (2/3)")
    log(f"   empirical SQUARED-class-mass ratio (p_odd)^2:(p_even)^2 = {sq_num:.4f}:1  (target 4:1)")
    log("")

    ok = (e2 == F(1,3) and o2 == F(2,3) and sq_ratio == 4 and e4 == F(1,15) and o4 == F(4,15)
          and abs(v1[1]-4) < 1e-9 and abs(p_even-1/3) < 2e-3 and abs(sq_num-4) < 0.02)
    log("## VERDICT")
    if ok:
        log("   CONFIRMED. R64.B's 1:4 class-mass ratio is ELEMENTARY and EXACT: the parity")
        log("   split of Geom(1/2) — P(v even)=1/3, P(v odd)=2/3 (squared -> 1:4), equivalently")
        log("   Sum_{even}4^{-v}:Sum_{odd}4^{-v}=1/15:4/15=1:4 (the T_diag prefactors). The value")
        log("   is correct; the provenance gap (no result_64*.md) was write-up only, now closed.")
        log("   c=7/45's dependence on R64.B is on an elementary Geom(1/2) parity identity.")
    else:
        log("   MISMATCH -- see checks above; R64.B value is NOT the naive parity split.")
    with open(r"C:\Collatz\result_64B_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("")
    log("[wrote] result_64B_log.txt")

if __name__ == "__main__":
    main()
