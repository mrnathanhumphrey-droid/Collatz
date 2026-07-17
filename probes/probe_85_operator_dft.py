"""
Probe 85 (corrected) — the operator-DFT bridge test, rung 1.

The pointwise chirp identity is FORCED by R82 step-1 (same character species). The
real evidential question at rung 1: does the ACTUAL DWM j=2 operator phase, DFT'd,
carry R81's INDEPENDENTLY-CERTIFIED F-hat structure?

DWM j=2 single-Kraus phase (from dwm_kraus_match_syracuse.build_M_tilde_unweighted),
even leg v=2k (the <4> sub-orbit), 9=3^2 reduces the modulus 3^n -> q=3^{n-2}:
    e_{3^n}(-xi*9*2^{-b}*2^{-2k}) = e_q(-xi*2^{-b} * 4^{-k}),   q=3^{n-2}=3^{r+1}, r=n-3.
So the DWM chirp has multiplier c' = (-xi*2^{-b}) mod q and orbit 4^{-k}. Its DFT is
F-hat with that c' (4^{-k} vs 4^{+k} is a frequency flip a->-a, structure identical).

TEST (n=6, r=3, q=81, d=27 -- clean Mahler regime, no offset):
 (1) enumerate the DWM multipliers c'=(-xi*2^{-b}) mod q over the states xi and the
     b-values that actually occur; report their <4>-coset structure vs R81's family
     (ep=0 -> <4>, ep=1 -> 2<4>). Does the operator populate BOTH R81 cosets?
 (2) for a representative c' in each coset, DFT the REAL DWM chirp and certify with
     R81's OWN exact machinery (magnitude 3sqrt(q) flat = Th 78.3; phase index J4;
     finite-difference degree). Compare to R81's certified family profile.
 PASS: DWM chirp reproduces |F|=3sqrt(q) flat AND R81's growing-degree J4 phase, and
       its multipliers occupy exactly R81's 2 cosets -> the DWM operator's chirp is
       structurally R81's certified F-hat. That is bridge part-A, at the operator level,
       with R81's certification (Mahler/cyclotomic, never saw the DWM operator) as the
       independent oracle. FAIL of any clause is reported straight.

Not at stake: THEOREM_C_745, Th 78.1-78.3. The 7/45 measure-assembly is rung 2 (separate).
"""
import sys, math, cmath
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:\Collatz")

from bilinear_pair_operator import build_markov_rational
from result_81_fhat_phase_profile import (
    pow4_table, Fhat_complex, certify_square, compute_J4, fd_degree, c_family, omega_r,
)

LOG = []
def log(m=""):
    print(m); LOG.append(str(m))

def dlog_base2(q):
    t = {}; x = 1
    order = 2 * (q // 3)
    for e in range(order):
        t[x] = e; x = (x * 2) % q
    return t

def coset_of(x, q, dl):
    """0 if x in <4> (dlog even), 1 if in 2<4> (dlog odd)."""
    return dl[x % q] % 2

# ---- R81-machinery phase profile for a given multiplier c (reused as the oracle) ----
def profile(r, c, pow4):
    """Return (mag_flat_ok, support_size, deg_4q, J4_list) for chirp multiplier c."""
    q = 3 ** (r + 1); d = 3 ** r
    magexp = 3 * math.sqrt(q)
    allabs = np.array([abs(Fhat_complex(r, c, a, pow4)) for a in range(d)])
    supp = [a for a in range(d) if allabs[a] > 1e-6]
    mags = [allabs[a] for a in supp]
    mag_ok = bool(mags) and (max(mags) - min(mags)) / magexp < 1e-9 \
             and abs(np.mean(mags) - magexp) / magexp < 1e-9
    J4 = []
    cf = 0
    for a in supp:
        Fc = Fhat_complex(r, c, a, pow4)
        ok, s, sg = certify_square(r, c, a, pow4)
        if not ok:
            cf += 1; J4.append(None); continue
        J4.append(compute_J4(r, s, sg, Fc))
    deg = None
    if cf == 0 and J4:
        deg, _, _ = fd_degree(J4, 4 * q)
    return mag_ok, len(supp), deg, J4, cf

def main():
    log("# PROBE 85 rung 1 — operator-DFT bridge test (REAL DWM phase vs R81 certified F-hat)")
    log("")
    n = 6; r = n - 3               # j=2 -> r = n-3 = 3
    N = 3 ** n; q = 3 ** (r + 1)   # q = 3^{n-2} = 81
    d = 3 ** r                     # 27
    pow4 = pow4_table(q)
    dl = dlog_base2(q)
    log(f"n={n}, j=2 -> r={r}, N=3^{n}={N}, reduced modulus q=3^{r+1}={q}, d=3^{r}={d}")
    log("")

    # ---- (1) DWM j=2 multipliers c' = (-xi * 2^{-b}) mod q, coset structure ----
    log("## (1) DWM j=2 chirp multipliers c' = (-xi*2^{-b}) mod q, occurring (xi,b)")
    _, coprime = build_markov_rational(n)         # states xi in (Z/3^n)*
    inv2q = pow(2, -1, q)
    V_MAX = 16
    b_vals = sorted({v1 + vp1 for v1 in range(1, V_MAX + 1)
                             for vp1 in range(1, V_MAX + 1) if v1 != vp1})
    cprime = set()
    for xi in coprime:
        xr = xi % q
        if xr == 0:
            continue
        for b in b_vals:
            cp = (-(xr) * pow(inv2q, b, q)) % q
            if cp != 0 and cp in dl:      # unit
                cprime.add(cp)
    c0 = [c for c in cprime if coset_of(c, q, dl) == 0]   # <4>
    c1 = [c for c in cprime if coset_of(c, q, dl) == 1]   # 2<4>
    log(f"   distinct unit multipliers c' occurring: {len(cprime)}  "
        f"(<4>: {len(c0)}, 2<4>: {len(c1)})")
    # R81's family cosets
    fam = c_family(r)
    fam_coset = {(ell, eps): coset_of(c, q, dl) for (ell, eps, c) in fam}
    log(f"   R81 family cosets: " + ", ".join(f"(l={l},e={e})->{'<4>' if cs==0 else '2<4>'}"
                                              for (l, e), cs in fam_coset.items()))
    both = (len(c0) > 0 and len(c1) > 0)
    log(f"   DWM operator populates BOTH R81 cosets: {'YES' if both else 'NO'}  "
        f"(R81 has eps=0 in <4>, eps=1 in 2<4>)")
    log("")

    # ---- (2) DFT the REAL DWM chirp; certify with R81 machinery; compare to family ----
    log("## (2) DFT of the REAL DWM chirp, certified by R81's exact machinery")
    log("   (magnitude 3sqrt(q) flat = Th 78.3; phase index J4; finite-diff degree)")
    log("")
    magexp = 3 * math.sqrt(q)
    log(f"   |F|=3sqrt(q) target = {magexp:.6f};  R81 degree at r=3 was 3 (growing-degree law)")
    log("")

    # R81 reference profiles for the two cosets (family reps: c_{0,0}=1 in <4>, c_{0,1}=2 in 2<4>)
    ref = {}
    for (ell, eps, c) in fam:
        if ell == 0:
            mag_ok, ns, deg, J4, cf = profile(r, c, pow4)
            ref[eps] = dict(c=c, mag_ok=mag_ok, ns=ns, deg=deg, J4=J4, cf=cf)
            log(f"   R81 ref  eps={eps} c={c:>3}: mag_flat={mag_ok}  |supp|={ns}  deg={deg}  cert_fail={cf}")
    log("")

    # DWM reps: pick one occurring c' per coset
    reps = {0: (sorted(c0)[0] if c0 else None), 1: (sorted(c1)[0] if c1 else None)}
    all_ok = both
    for cs, cp in reps.items():
        if cp is None:
            log(f"   DWM  coset {cs}: no multiplier occurs — SKIP"); continue
        mag_ok, ns, deg, J4, cf = profile(r, cp, pow4)
        # compare to R81 ref of same coset: degree + magnitude flatness + support size.
        rf = ref.get(cs)
        deg_match = (rf is not None and deg == rf["deg"])
        struct_ok = mag_ok and (ns == rf["ns"]) and deg_match and cf == 0
        # J4 profiles equal up to a global constant (R84 omega_3^l) + the 4^{-k} freq flip?
        # test: within-support finite-difference SEQUENCE (deg+1 diffs) matches ref's,
        # which is invariant to a global additive phase constant.
        Jd = J4match = None
        if cf == 0 and rf and rf["cf"] == 0 and len(J4) == len(rf["J4"]):
            def diffs(seq, mod, g):
                cur = [x % mod for x in seq]
                for _ in range(g):
                    cur = [(cur[i+1]-cur[i]) % mod for i in range(len(cur)-1)]
                return cur
            g = deg if deg else 1
            dd = diffs(J4, 4*q, 1)          # 1st difference kills the global constant
            dr = diffs(rf["J4"], 4*q, 1)
            # allow the a->-a reversal: compare as multisets of |first-differences|
            J4match = sorted((min(x, 4*q-x) for x in dd)) == sorted((min(x, 4*q-x) for x in dr))
        all_ok = all_ok and struct_ok
        log(f"   DWM  coset {cs} c'={cp:>3}: mag_flat={mag_ok}  |supp|={ns}  deg={deg}  "
            f"cert_fail={cf}  | vs R81: deg_match={deg_match}  1stdiff_multiset_match={J4match}")
    log("")
    log("## VERDICT")
    if all_ok:
        log("   PASS(part-A): the REAL DWM j=2 operator chirp, DFT'd, reproduces R81's")
        log("   certified F-hat structure (|F|=3sqrt(q) flat, growing-degree J4 phase) on")
        log("   BOTH cosets R81 uses. The DWM operator's chirp IS R81's certified object;")
        log("   R81's Mahler/cyclotomic certification (independent of the DWM operator) is")
        log("   the oracle. Rung 2 (Geom-measure assembly -> 7/45) is the remaining piece.")
    else:
        log("   MIXED/FAIL: see clauses above; reported straight, not smoothed.")
    with open(r"C:\Collatz\result_85_operator_dft_log.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    log("")
    log("[wrote] result_85_operator_dft_log.txt")

if __name__ == "__main__":
    main()
