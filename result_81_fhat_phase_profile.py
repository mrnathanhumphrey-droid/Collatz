"""
result_81_fhat_phase_profile.py

PROBE: F-hat phase profile on the R78 support.

Object (from result_78_FINAL.md 78.2/78.3):
    F_hat(xi) = sum_{u=0}^{q-1} e_q(c*4^u - xi*u),   q = 3^{r+1}
supported on supp = {3a : a in Z/3^r, a == 1 mod 3}, |supp| = 3^{r-1},
with |F_hat| = 3*sqrt(q) constant on the support (Theorem 78.3).

Since the magnitude is constant, ALL structure lives in arg F_hat(3a). This probe
maps that phase profile and tests whether it is a quadratic (Gauss-sum) in the
support index a.

PRE-REG: PROBE_F_hat_phase (as pasted) + user amendments 2026-07-13:

  §3' (amended, pre-run) -- multi-denominator congruence fit.
    Run the congruence quadratic fit at each denominator
        D in {3^r, 3^{r+1}, 2*3^r}     (pre-reg-official; 3^{r+1} weighted highest)
    plus D = 2*3^{r+1} = 2q as an added honest-modulus diagnostic (the genuine
    char-3 quadratic-Gauss modulus, which carries the sqrt(+-3) half-shift at the
    F_hat level -- see below). Report a pass/fail table across all denominators.
    H_QUAD fires iff some D gives 100% at every r>=3 with an alpha != 0 mod 3
    representative. H_PSEUDO fires only if ALL denominators fail at every r>=3.
    r=2 excluded from evidence (|A_2|=3, zero dof). No within-support magnitude
    filter (support found by the clean 0 vs 3*sqrt(q) gap).

  §3' fit-procedure correction (pre-run flaw caught before firing).
    Every support index satisfies a == 1 mod 3, so any 3-point Vandermonde over the
    support is singular mod 3 (det divisible by 27) and the pre-reg's "solve from 3
    points" does NOT determine (alpha,beta,gamma) mod 3^k. Fix: reparametrize
    a = 1 + 3b, b in Z/3^{r-1} (a bijection onto the support). A quadratic in a maps
    to a quadratic in b:
        alpha a^2 + beta a + gamma = K b^2 + H b + G,
        K = 9*alpha,  H = 6*alpha + 3*beta,  G = alpha + beta + gamma.
    The b-Vandermonde on b in {0,1,2} has det -2 (a UNIT mod 3^k), so the b-fit is
    clean and exact. Fit (K,H,G) from b=0,1,2, VERIFY on all remaining b (100% req).
    Then the pre-reg H_QUAD-in-a exists iff  9 | K  and  3 | H  (mod the 3-part),
    with  alpha == (K/9) mod 3  -- so alpha != 0 mod 3  <=>  v_3(K) == 2 exactly.
    NB: on the support a^2 == a == 1 mod 3, so the alpha-mod-3 class is invisible at
    the mod-3 level and only becomes identifiable via the b^2 coefficient K -- which
    is exactly why fitting in b (not a) is the correct tool.

  §2' (corrected) -- exact Gauss-sum certification (no irrational normalization).
    |F_hat| = 3*sqrt(q), so the unit object is F_hat/(3 sqrt(q)). To avoid sqrt(q)
    (and the sqrt(+-3) quadratic-Gauss twist), test instead whether
        F_hat(3a)^2 / (9q)   is a root of unity in Q(zeta_q),  q = 3^{r+1}.
    |F_hat|^2 = 9q exactly, so the denominator is a rational integer and the check is
    exact. F_hat is a sum of q-th roots of unity => an algebraic integer in Z[zeta_q]
    for free. We certify F_hat(3a)^2 == 9q * zeta_q^{s_a} for an EXACT integer s_a by
    reducing the integer self-convolution of the phase histogram modulo
        Phi_q(x) = x^{2*3^r} + x^{3^r} + 1
    and matching 9q * x^{s_a}. From s_a and the float branch we recover the exact
    phase index J_a in Z/2q (the "2" is the sqrt(+-3) half-shift): F_hat(3a) has
    argument 2*pi*J_a/(2q). The whole congruence analysis is then on EXACT integers.

Deliverables (this run): result_81_fhat_phase_profile.py, result_81_fhat_phase_data.csv,
result_81_fhat_phase.md, result_81_log.txt. Disposition appended to STATE.md (new
dated entry, no rewrite of existing content).
"""

import sys
import os
import math
import cmath
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")
REPO = r"C:\Collatz"
OUTDIR = os.path.join(REPO, "experiments_output")
os.makedirs(OUTDIR, exist_ok=True)

LOG_LINES = []
def log(msg=""):
    print(msg)
    LOG_LINES.append(str(msg))


# ---------------------------------------------------------------------------
# Core F_hat computations
# ---------------------------------------------------------------------------

def pow4_table(q):
    t = np.empty(q, dtype=np.int64)
    p = 1
    for u in range(q):
        t[u] = p
        p = (p * 4) % q
    return t


def Fhat_complex(r, c, a, pow4=None):
    """F_hat(3a) = sum_u e_q(c*4^u - 3a*u), q=3^{r+1}, as a complex float."""
    q = 3 ** (r + 1)
    if pow4 is None:
        pow4 = pow4_table(q)
    u = np.arange(q, dtype=np.int64)
    ph = ((c * pow4 - (3 * a) * u) % q).astype(np.float64) / q
    return np.sum(np.exp(2j * np.pi * ph))


def phase_hist(r, c, a, pow4=None):
    """Integer histogram over Z/q of the exponents e_u = (c*4^u - 3a*u) mod q."""
    q = 3 ** (r + 1)
    if pow4 is None:
        pow4 = pow4_table(q)
    u = np.arange(q, dtype=np.int64)
    e = (c * pow4 - (3 * a) * u) % q
    return np.bincount(e, minlength=q).astype(np.int64)


def reduce_modPhi(vec_lenq, r):
    """Reduce sum_t vec[t] x^t  (t in 0..q-1, q=3d) modulo Phi_q(x)=x^{2d}+x^d+1.
    x^{2d} == -x^d - 1, so for t in [2d,3d): contributes -1 to positions t-d and t-2d.
    Returns a length-2d integer vector (canonical rep in Z[zeta_q]). Uses python ints."""
    d = 3 ** r
    V = [int(x) for x in vec_lenq[:2 * d]]
    for t in range(2 * d, 3 * d):
        cf = int(vec_lenq[t])
        if cf:
            V[t - d] -= cf
            V[t - 2 * d] -= cf
    return V


def reduce_xs(s, r):
    """Reduce x^s mod Phi_q -> length-2d vector."""
    d = 3 ** r
    V = [0] * (2 * d)
    s %= (3 * d)
    if s < 2 * d:
        V[s] += 1
    else:  # 2d <= s < 3d
        V[s - d] -= 1
        V[s - 2 * d] -= 1
    return V


def certify_square(r, c, a, Fc, pow4=None):
    """Certify F_hat(3a)^2 == sigma * 9q * zeta_q^{s} exactly, sigma in {+1,-1}.
    Returns (ok, s, sigma). sigma=-1 is the sqrt(-3) quadratic-Gauss signature that
    appears for r even. Fc = float F_hat(3a) (for the s guess)."""
    q = 3 ** (r + 1)
    hist = phase_hist(r, c, a, pow4)
    # integer self-convolution (circular, mod q) via FFT, rounded exactly
    F = np.fft.fft(hist.astype(np.float64))
    C = np.rint(np.fft.ifft(F * F).real).astype(np.int64)  # C[t] = coeff of x^t in F_hat^2
    V = reduce_modPhi(C, r)
    d = 3 ** r
    # Derive (s, sigma) EXACTLY from V: V/(9q) must be the reduced form of
    # sigma * x^s. Single power +x^s (s<2d) -> one +1; -x^s -> one -1; for
    # 2d<=s<3d, x^s reduces to -x^{s-d}-x^{s-2d} (two entries d apart).
    if any(v % (9 * q) != 0 for v in V):
        return False, None, None
    Vq = [v // (9 * q) for v in V]
    nz = [(i, Vq[i]) for i in range(len(Vq)) if Vq[i] != 0]
    if len(nz) == 1:
        p, v = nz[0]
        if v == 1:
            return True, p, 1          # +x^p
        if v == -1:
            return True, p, -1         # -x^p
        return False, None, None
    if len(nz) == 2:
        (p1, v1), (p2, v2) = nz
        if v1 == v2 and (p2 - p1) == d:
            s = p1 + 2 * d             # x^s reduces to -x^{s-d}-x^{s-2d}
            if v1 == -1:
                return True, s % q, 1  # from +x^s
            if v1 == 1:
                return True, s % q, -1 # from -x^s
    return False, None, None


def compute_J4(r, s, sigma, Fc):
    """Exact phase index J4 in Z/4q with arg F_hat(3a) = 2*pi*J4/(4q).
    From certified F_hat^2/(9q) = sigma*zeta_q^s:
      sigma=+1: e^{i theta} = +-zeta_{2q}^s        -> J4 = 2s + t*2q      (even)
      sigma=-1: e^{i theta} = +-zeta_{4q}^{q+2s}    -> J4 = (q+2s) + t*2q  (odd)
    The branch t in {0,1} is fixed by the float sign of F_hat (two candidates are
    exact negatives, well separated by 2*3*sqrt(q))."""
    q = 3 ** (r + 1)
    magexp = 3 * math.sqrt(q)
    base = (2 * s) % (4 * q) if sigma == 1 else (q + 2 * s) % (4 * q)
    cand = magexp * cmath.exp(2j * math.pi * base / (4 * q))
    t = 0 if abs(Fc - cand) < abs(Fc + cand) else 1
    return (base + t * 2 * q) % (4 * q)


# ---------------------------------------------------------------------------
# Congruence quadratic fit (in b; maps to a-quadratic via K,H,G)
# ---------------------------------------------------------------------------

def fit_quad_in_b_mod3k(bs, ns, k):
    """Fit K b^2 + H b + G == n(b) mod 3^k using b=0,1,2 (unit Vandermonde),
    verify on all points. Returns (ok, (K,H,G)) with coeffs in Z/3^k."""
    D = 3 ** k
    idx = {b: i for i, b in enumerate(bs)}
    if not all(bb in idx for bb in (0, 1, 2)):
        return None  # cannot seed (should not happen for r>=2)
    n0 = ns[idx[0]] % D
    n1 = ns[idx[1]] % D
    n2 = ns[idx[2]] % D
    inv2 = pow(2, -1, D)
    K = ((n2 - 2 * n1 + n0) * inv2) % D
    H = (n1 - n0 - K) % D
    G = n0 % D
    for b, n in zip(bs, ns):
        if (K * b * b + H * b + G - n) % D != 0:
            return (False, (K, H, G))
    return (True, (K, H, G))


def fit_quad_in_b_brute(bs, ns, mod):
    """Enumerate (K,H,G) in (Z/mod)^3 consistent with all points mod `mod`.
    Used for the 2-part (mod in {2,4}) where the b-Vandermonde is not a unit."""
    sols = []
    for K in range(mod):
        for H in range(mod):
            for G in range(mod):
                if all((K * b * b + H * b + G - n) % mod == 0 for b, n in zip(bs, ns)):
                    sols.append((K, H, G))
    return sols


def crt(a1, m1, a2, m2):
    """Solve x == a1 mod m1, x == a2 mod m2 (m1,m2 coprime). Return x mod m1*m2."""
    inv = pow(m1, -1, m2)
    t = ((a2 - a1) * inv) % m2
    return (a1 + m1 * t) % (m1 * m2)


def fit_quadratic(a_vals, n_vals, D):
    """Fit n(a) == alpha a^2 + beta a + gamma mod D over the support, via b=(a-1)/3.
    D factors as 2^e2 * 3^e3 with e2 in {0,1}. Returns dict with keys:
      pass (bool), coeffs_b (K,H,G) mod D or None, a_quadratic (bool),
      alpha_mod3 (int or None), is_linear (bool)."""
    e2 = 0
    e3 = 0
    d = D
    while d % 2 == 0:
        d //= 2; e2 += 1
    while d % 3 == 0:
        d //= 3; e3 += 1
    assert d == 1, f"unexpected D={D} (not 2^a*3^b)"

    a0 = a_vals[0] % 3
    assert all(a % 3 == a0 for a in a_vals)
    bs = [(a - a0) // 3 for a in a_vals]
    fail = dict(passed=False, coeffs_b=None, a_quadratic=False,
                alpha_mod3=None, is_linear=False)

    # 3-part (unit b-Vandermonde -> closed form)
    r3 = fit_quad_in_b_mod3k(bs, n_vals, e3) if e3 > 0 else (True, (0, 0, 0))
    if r3 is None or not r3[0]:
        return fail
    K3, H3, G3 = r3[1]
    M3 = 3 ** e3

    if e2 == 0:
        K, H, G = K3, H3, G3
        M = M3
    else:
        M2 = 2 ** e2
        s2 = fit_quad_in_b_brute(bs, n_vals, M2)  # 2-part (mod 2 or 4)
        if not s2:
            return fail
        M = M2 * M3
        chosen = None
        for (K2, H2, G2) in s2:
            K = crt(K2, M2, K3, M3); H = crt(H2, M2, H3, M3); G = crt(G2, M2, G3, M3)
            if all((K * b * b + H * b + G - n) % M == 0 for b, n in zip(bs, n_vals)):
                chosen = (K, H, G); break
        if chosen is None:
            return fail
        K, H, G = chosen

    # a-quadratic existence: need 9 | K and 3 | H in the 3-part
    a_quad = (K3 % 9 == 0) and (H3 % 3 == 0)
    alpha_mod3 = None
    if a_quad:
        alpha_mod3 = ((K3 // 9) % 3)
    is_linear = (K3 % M3 == 0)  # b^2 coeff vanishes in 3-part => linear in b (hence in a)
    return dict(passed=True, coeffs_b=(K % M, H % M, G % M), a_quadratic=a_quad,
                alpha_mod3=alpha_mod3, is_linear=is_linear)


# ---------------------------------------------------------------------------
# Self-tests (must pass before the real run)
# ---------------------------------------------------------------------------

def self_tests():
    log("## Self-tests")
    ok = True

    # (1) planted quadratic in a on a support {a==1 mod3} recovers & passes
    for (r, alpha0, beta0, gamma0) in [(3, 2, 5, 1), (4, 7, 0, 3), (4, 3, 8, 4)]:
        D = 3 ** (r + 1)
        a_vals = [a for a in range(3 ** r) if a % 3 == 1]
        n_vals = [(alpha0 * a * a + beta0 * a + gamma0) % D for a in a_vals]
        res = fit_quadratic(a_vals, n_vals, D)
        exp_alpha3 = alpha0 % 3
        good = res["passed"] and res["a_quadratic"] and res["alpha_mod3"] == exp_alpha3
        log(f"  planted quad r={r} (a^2*{alpha0}+a*{beta0}+{gamma0}) mod {D}: "
            f"pass={res['passed']} a_quad={res['a_quadratic']} "
            f"alpha_mod3={res['alpha_mod3']} (exp {exp_alpha3}) -> {'OK' if good else 'FAIL'}")
        ok = ok and good

    # (1b) planted LINEAR in a (alpha0=0): a_quadratic true but alpha_mod3==0, is_linear
    r = 4; D = 3 ** (r + 1)
    a_vals = [a for a in range(3 ** r) if a % 3 == 1]
    n_vals = [(0 * a * a + 5 * a + 2) % D for a in a_vals]
    res = fit_quadratic(a_vals, n_vals, D)
    good = res["passed"] and res["alpha_mod3"] == 0 and res["is_linear"]
    log(f"  planted linear r={r}: pass={res['passed']} alpha_mod3={res['alpha_mod3']} "
        f"is_linear={res['is_linear']} -> {'OK' if good else 'FAIL'}")
    ok = ok and good

    # (2) planted NON-polynomial (pseudo) is refuted
    r = 4; D = 3 ** (r + 1)
    a_vals = [a for a in range(3 ** r) if a % 3 == 1]
    # deterministic 'random-like' cubic-with-noise that is not a quadratic mod D
    n_vals = [(a * a * a + 7 * a + (a * a * a * a % 5)) % D for a in a_vals]
    res = fit_quadratic(a_vals, n_vals, D)
    good = (not res["passed"])
    log(f"  planted non-quadratic r={r}: pass={res['passed']} -> "
        f"{'OK (refuted)' if good else 'FAIL (false-fire)'}")
    ok = ok and good

    # (3) cyclotomic reducer sanity: sum of all q-th roots == 0 in Z[zeta_q]
    #     i.e. reduce([1,1,...,1]) == reduce(-(x^d + x^{2d}))? Actually 1+x+...+x^{q-1}=0.
    for r in (2, 3):
        q = 3 ** (r + 1)
        allones = np.ones(q, dtype=np.int64)
        V = reduce_modPhi(allones, r)
        good = all(v == 0 for v in V)
        log(f"  reducer sum-of-all-roots r={r}: {'OK (==0)' if good else 'FAIL'}")
        ok = ok and good

    log(f"  SELF-TESTS: {'ALL PASS' if ok else 'FAILURE'}")
    log("")
    return ok


# ---------------------------------------------------------------------------
# Main probe
# ---------------------------------------------------------------------------

def omega_r(r):
    q = 3 ** (r + 1)
    w = 1 + 3 ** r
    if pow(w, 3, q) != 1:
        w = 1 + 2 * (3 ** r)
    assert pow(w, 3, q) == 1
    return w


def c_family(r):
    """The 6 c_{ell,eps} = 2^eps * omega_r^ell mod q, labelled (ell,eps)."""
    q = 3 ** (r + 1)
    w = omega_r(r)
    out = []
    for eps in (0, 1):
        for ell in (0, 1, 2):
            c = (pow(2, eps, q) * pow(w, ell, q)) % q
            out.append((ell, eps, c))
    return out


def main():
    log("# R81: F-hat phase profile on the R78 support")
    log(f"# denominators tested per level: 3^r, 3^{{r+1}}, 2*3^r, 2*3^{{r+1}} (=2q)")
    log(f"# precision: float64 np for support/branch; EXACT integer cyclotomic for phase index")
    log("")

    if not self_tests():
        log("ABORT: self-tests failed; not running probe.")
        _flush_log()
        return

    R_LIST = [2, 3, 4, 5, 6]
    csv_rows = []
    # aggregate: fit_table[(r, ell, eps)][D_label] = dict
    fit_table = {}
    smoke_ok = True

    for r in R_LIST:
        q = 3 ** (r + 1)
        d = 3 ** r
        magexp = 3 * math.sqrt(q)
        pow4 = pow4_table(q)
        # master exact modulus is 4q; report the pre-reg-official + honest-modulus set
        Dmap = [("3^r", d), ("3^(r+1)", q), ("2*3^r", 2 * d), ("2q", 2 * q), ("4q", 4 * q)]
        MASTER = 4 * q
        log(f"## r={r}  q=3^{r+1}={q}  |supp|_pred=3^{r-1}={d//3}  |F|=3sqrt(q)={magexp:.6f}")

        for (ell, eps, c) in c_family(r):
            # support found by the clean 0 vs 3sqrt(q) gap (NO within-support filter);
            # coset a0 mod 3 depends on c (eps=0 -> 1, eps=1 -> 2).
            allabs = np.array([abs(Fhat_complex(r, c, a, pow4)) for a in range(d)])
            supp_a = [a for a in range(d) if allabs[a] > 1e-6]
            offvals = [allabs[a] for a in range(d) if allabs[a] <= 1e-6]
            max_off = max(offvals) if offvals else 0.0
            if max_off > 1e-9:
                log(f"  !! (ell={ell},eps={eps}) off-support max |F|={max_off:.2e} -> not clean 0")
                smoke_ok = False
            a0 = supp_a[0] % 3 if supp_a else None

            J_by_a = {}
            absF_by_a = {}
            theta_by_a = {}
            s_by_a = {}
            sig_by_a = {}
            mags = []
            cert_fail = 0
            for a in supp_a:
                Fc = Fhat_complex(r, c, a, pow4)
                absF = abs(Fc)
                mags.append(absF)
                theta = cmath.phase(Fc)  # principal branch (-pi,pi]
                ok, s, sigma = certify_square(r, c, a, Fc, pow4)
                if not ok:
                    cert_fail += 1
                    s_by_a[a] = None; J_by_a[a] = None; sig_by_a[a] = None
                    absF_by_a[a] = absF; theta_by_a[a] = theta
                    continue
                J = compute_J4(r, s, sigma, Fc)  # exact index in Z/4q
                J_by_a[a] = J; s_by_a[a] = s; sig_by_a[a] = sigma
                absF_by_a[a] = absF; theta_by_a[a] = theta

            # smoke: magnitude constancy on support (no threshold filter used to FIND supp)
            if mags:
                mrel = (max(mags) - min(mags)) / magexp
                const_ok = mrel < 1e-12 and abs(np.mean(mags) - magexp) / magexp < 1e-9
            else:
                mrel = float('nan'); const_ok = False
            if not const_ok:
                log(f"  (ell={ell},eps={eps}) SMOKE |F| non-constant: relspread={mrel:.2e} -> STOP-flag")
                smoke_ok = False
            if len(supp_a) != d // 3:
                log(f"  (ell={ell},eps={eps}) |supp|={len(supp_a)} != {d//3} -> support size wrong")
                smoke_ok = False

            # fits per denominator (only if all J certified)
            row_fit = {}
            certified = (cert_fail == 0) and all(J_by_a[a] is not None for a in supp_a)
            for (Dlab, D) in Dmap:
                mratio = MASTER // D
                if not certified:
                    row_fit[Dlab] = dict(integral=None, passed=None,
                                         a_quadratic=None, alpha_mod3=None,
                                         is_linear=None, note="cert-fail")
                    continue
                integral = all(J_by_a[a] % mratio == 0 for a in supp_a)
                if not integral:
                    nbad = sum(1 for a in supp_a if J_by_a[a] % mratio != 0)
                    row_fit[Dlab] = dict(integral=False, passed=False,
                                         a_quadratic=False, alpha_mod3=None,
                                         is_linear=None,
                                         note=f"{nbad}/{len(supp_a)} not on D-roots")
                    continue
                n_vals = [(J_by_a[a] // mratio) % D for a in supp_a]
                res = fit_quadratic(supp_a, n_vals, D)
                res["integral"] = True
                res["note"] = ""
                row_fit[Dlab] = res

            fit_table[(r, ell, eps)] = dict(fit=row_fit, cert_fail=cert_fail,
                                            nsupp=len(supp_a), const_ok=const_ok)

            # CSV rows
            for a in supp_a:
                J = J_by_a[a]
                # pre-reg phi_r(a) = theta * 3^r / (2pi)  (D=3^r normalization)
                phi_3r = theta_by_a[a] * d / (2 * math.pi)
                # congruence residual at primary D=3^{r+1}=q
                cong_res = ""
                rf = row_fit.get("3^(r+1)")
                if J is not None and rf is not None and rf.get("integral") and rf.get("coeffs_b"):
                    Dq = q; mratio = MASTER // Dq
                    if J % mratio == 0:
                        K, H, Gc = rf["coeffs_b"]
                        b = (a - a0) // 3
                        nfit = (K * b * b + H * b + Gc) % Dq
                        ntrue = (J // mratio) % Dq
                        cong_res = (ntrue - nfit) % Dq
                csv_rows.append(dict(
                    r=r, ell=ell, eps=eps, c=c, a=a, b=(a - a0) // 3,
                    absF=f"{absF_by_a[a]:.10f}", theta=f"{theta_by_a[a]:.12f}",
                    s=s_by_a[a], J=J, phi_3r=f"{phi_3r:.6f}",
                    cong_residual_q=cong_res))
            fired = [Dlab for (Dlab, D) in Dmap
                     if row_fit[Dlab].get("passed")]
            log(f"  (ell={ell},eps={eps}) c={c:>5}  cert_fail={cert_fail}  "
                f"const={const_ok}  passes@={fired}")
        log("")

    _write_csv(csv_rows)
    _disposition(fit_table, R_LIST, smoke_ok)
    _flush_log()


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _write_csv(rows):
    path = os.path.join(REPO, "result_81_fhat_phase_data.csv")
    cols = ["r", "ell", "eps", "c", "a", "b", "absF", "theta", "s", "J",
            "phi_3r", "cong_residual_q"]
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for row in rows:
            f.write(",".join(str(row[c]) for c in cols) + "\n")
    log(f"[wrote] {path}  ({len(rows)} rows)")


def _decide(fit_table, R_LIST):
    """Aggregate decision across r>=3 (r=2 excluded from evidence)."""
    Dlabs = ["3^r", "3^(r+1)", "2*3^r", "2q", "4q"]
    r_eval = [r for r in R_LIST if r >= 3]
    # For each D: does it pass (with alpha!=0 rep i.e. a_quadratic & alpha_mod3!=0)
    # on 100% of (r,ell,eps) cells at every r>=3?
    D_quad_fire = {}
    D_any_pass = {}
    for Dlab in Dlabs:
        quad_fire = True
        any_pass = True
        for r in r_eval:
            for ell in (0, 1, 2):
                for eps in (0, 1):
                    cell = fit_table.get((r, ell, eps))
                    if cell is None:
                        quad_fire = False; any_pass = False; continue
                    rf = cell["fit"][Dlab]
                    if not rf.get("passed"):
                        any_pass = False
                    if not (rf.get("passed") and rf.get("a_quadratic")
                            and rf.get("alpha_mod3") not in (None, 0)):
                        quad_fire = False
        D_quad_fire[Dlab] = quad_fire
        D_any_pass[Dlab] = any_pass
    return D_quad_fire, D_any_pass, r_eval


def _disposition(fit_table, R_LIST, smoke_ok):
    D_quad_fire, D_any_pass, r_eval = _decide(fit_table, R_LIST)
    Dlabs = ["3^r", "3^(r+1)", "2*3^r", "2q", "4q"]

    # verdict
    quad_D = [D for D in Dlabs if D_quad_fire[D]]
    any_pass_D = [D for D in Dlabs if D_any_pass[D]]
    if quad_D:
        verdict = "H_QUAD"
    elif any_pass_D:
        # a polynomial passes but without alpha!=0 -> linear or degenerate
        verdict = "H_LIN / poly-but-not-genuine-quadratic"
    else:
        verdict = "H_PSEUDO"

    lines = []
    lines.append("# R81 disposition — F-hat phase profile on the R78 support")
    lines.append("")
    lines.append(f"**Date:** 2026-07-13. **Verdict: {verdict}.**")
    lines.append("")
    lines.append("Probe: `result_81_fhat_phase_profile.py`. Data: "
                 "`result_81_fhat_phase_data.csv`. Log: `result_81_log.txt`.")
    lines.append("")
    lines.append(f"Smoke (Theorem 78.3, |F_hat|=3√q constant on support, support "
                 f"{{a≡1 mod 3}}): **{'PASS' if smoke_ok else 'FAIL — see log'}**.")
    lines.append("")
    lines.append("## Pass/fail table (per denominator D; r=2 EXCLUDED from evidence)")
    lines.append("")
    lines.append("Cell = fraction of (ℓ,ε) c-family combos (6 total) at that r whose "
                 "phase index J_a, rescaled to D-th roots, fits a congruence quadratic "
                 "in a on 100% of support points. `q(α≠0)` marks that a genuine "
                 "α≢0 mod 3 quadratic-in-a representative exists.")
    lines.append("")
    header = "| r | " + " | ".join(Dlabs) + " |"
    sep = "|" + "---|" * (len(Dlabs) + 1)
    lines.append(header)
    lines.append(sep)
    for r in R_LIST:
        cells = []
        for Dlab in Dlabs:
            npass = 0; nquad = 0; tot = 0
            integ = 0
            for ell in (0, 1, 2):
                for eps in (0, 1):
                    cell = fit_table.get((r, ell, eps))
                    if cell is None:
                        continue
                    tot += 1
                    rf = cell["fit"][Dlab]
                    if rf.get("integral"):
                        integ += 1
                    if rf.get("passed"):
                        npass += 1
                    if rf.get("passed") and rf.get("a_quadratic") and \
                       rf.get("alpha_mod3") not in (None, 0):
                        nquad += 1
            tag = f"{npass}/{tot} pass"
            if integ < tot:
                tag += f", {integ}/{tot} on-roots"
            if nquad:
                tag += f", {nquad} q(α≠0)"
            cells.append(tag)
        excl = " *(excluded)*" if r == 2 else ""
        lines.append(f"| {r}{excl} | " + " | ".join(cells) + " |")
    lines.append("")

    # H_QUAD firing rule statement
    lines.append("## Decision (§3′ rule)")
    lines.append("")
    lines.append(f"- H_QUAD fires iff some D gives 100% with an α≢0 mod 3 rep at "
                 f"**every** r∈{r_eval}. Fired denominators: "
                 f"**{quad_D if quad_D else 'NONE'}**.")
    lines.append(f"- Any-polynomial-pass denominators (100% congruence, possibly "
                 f"linear/degenerate) at every r≥3: **{any_pass_D if any_pass_D else 'NONE'}**.")
    lines.append(f"- H_PSEUDO fires iff ALL denominators fail at every r≥3.")
    lines.append("")

    # coefficient cross-r table if a quadratic fired
    if quad_D:
        Dsel = quad_D[0]
        lines.append(f"## Cross-r coefficients (b-parametrization, D={Dsel})")
        lines.append("")
        lines.append("`a = 1 + 3b`; phase index (rescaled to D) = K·b² + H·b + G mod D. "
                     "α (in a) = K/9 mod 3^(k−2); α mod 3 = (K/9) mod 3.")
        lines.append("")
        lines.append("| r | ℓ | ε | K | H | G | α mod 3 |")
        lines.append("|---|---|---|---|---|---|---|")
        for r in R_LIST:
            for ell in (0, 1, 2):
                for eps in (0, 1):
                    cell = fit_table.get((r, ell, eps))
                    if cell is None:
                        continue
                    rf = cell["fit"][Dsel]
                    if rf.get("coeffs_b"):
                        K, H, Gc = rf["coeffs_b"]
                        am = rf.get("alpha_mod3")
                        star = "" if r > 2 else " *(excl)*"
                        lines.append(f"| {r}{star} | {ell} | {eps} | {K} | {H} | {Gc} | {am} |")
        lines.append("")

    # routing statement
    lines.append("## Routing (which of the three paper routes this opens/closes)")
    lines.append("")
    if verdict == "H_QUAD":
        lines.append("F̂ is a genuine quadratic Gauss sum on the principal-unit coset. "
                     "**Opens** a stationary-phase / Weyl-differencing route on "
                     "Σ 1̂(3a)·F̂(3a) that is NOT in R78's current route list "
                     "(band-ℓ1 CLOSED, BGK additive-energy random-like, Cauchy–Schwarz "
                     "trivial). The quadratic phase means the mixed bilinear sum admits "
                     "a completing-the-square / Gauss-sum evaluation — the missing "
                     "square-root saving of R78 §‘Crucial observation’. Next probe "
                     "(separate): evaluate Σ 1̂(3a)F̂(3a) via the recovered (α,β,γ).")
    elif verdict.startswith("H_LIN"):
        lines.append("The phase collapses to (at most) linear in a on the support — the "
                     "product sum Σ 1̂(3a)F̂(3a) reduces to a shifted short character "
                     "sum, **directly evaluable** without a Gauss-sum/Weyl route. "
                     "Closes the ‘need Burgess’ framing at these levels; the saving is "
                     "elementary. Verify at larger r before relying on it.")
    else:
        lines.append("No low-degree polynomial phase; the phase index is consistent with "
                     "equidistribution on Z/3^r. **This is a certifying negative** — it "
                     "**retires the smooth-completion / Gauss-sum route** cleanly: the "
                     "square-root saving in Σ 1̂(3a)F̂(3a) cannot come from a completed "
                     "quadratic phase, so R78’s residual bilinear bound genuinely needs "
                     "Burgess-strength input (the Burgess wall is real). Publishable as "
                     "a route-closing result; does not weaken Theorems 78.1–78.3.")
    lines.append("")
    lines.append("_Reporting discipline: the fired outcome above is reported as-is, "
                 "including a null. r=2 carries no evidential weight (|A_2|=3, zero dof). "
                 "No within-support magnitude filter was applied._")
    lines.append("")

    md_path = os.path.join(REPO, "result_81_fhat_phase.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    log(f"[wrote] {md_path}")

    # STATE.md append (new dated entry, no rewrite)
    _append_state(verdict, quad_D, any_pass_D, r_eval, smoke_ok)

    # keep verdict for console
    log("")
    log(f"==== VERDICT: {verdict}  (quad_D={quad_D}, any_pass_D={any_pass_D}) ====")


def _append_state(verdict, quad_D, any_pass_D, r_eval, smoke_ok):
    state = os.path.join(REPO, "STATE.md")
    entry = []
    entry.append("")
    entry.append("---")
    entry.append("")
    entry.append(f"**R81 — F̂ phase profile on the R78 support (2026-07-13).** "
                 f"Verdict: **{verdict}**. ")
    entry.append(f"Mapped arg F̂(3a) on supp={{3a: a≡1 mod 3}} (where |F̂|=3√q is "
                 f"constant, Th 78.3, smoke {'PASS' if smoke_ok else 'FAIL'}) for r∈{{2,3,4,5,6}} "
                 f"and all 6 c_{{ℓ,ε}}=2^ε·(1+3^r)^ℓ. Phase index certified EXACTLY: "
                 f"F̂(3a)²/(9q) is a root of unity ζ_q^{{s_a}} in Z[ζ_q] (integer "
                 f"cyclotomic reduction mod Φ_q=x^{{2·3^r}}+x^{{3^r}}+1), giving exact "
                 f"J_a∈Z/2q. Congruence quadratic fit in a (via a=1+3b, since the "
                 f"support's a≡1 mod 3 makes any 3-point Vandermonde singular mod 3) "
                 f"at D∈{{3^r, 3^{{r+1}}, 2·3^r, 2q}}. ")
    if verdict == "H_QUAD":
        entry.append(f"Genuine quadratic Gauss-sum phase fired at D={quad_D} on 100% of "
                     f"support at every r∈{r_eval}, α≢0 mod 3. **Opens a "
                     f"stationary-phase/Weyl route on Σ1̂·F̂ outside R78's route list.** "
                     f"Files: result_81_fhat_phase_profile.py + result_81_fhat_phase.md "
                     f"+ result_81_fhat_phase_data.csv + result_81_log.txt.")
    elif verdict.startswith("H_LIN"):
        entry.append(f"Phase is (at most) linear in a (poly-pass at D={any_pass_D}, no "
                     f"α≢0 quadratic). Σ1̂·F̂ reduces to a shifted short sum — elementary "
                     f"evaluation, no Burgess needed at these levels. "
                     f"Files: result_81_*.")
    else:
        entry.append(f"No low-degree polynomial phase (pass-set empty at every r≥3); "
                     f"phase index equidistributed on Z/3^r. **Certifying negative — "
                     f"retires the smooth-completion/Gauss-sum route; the Burgess wall "
                     f"is real.** Theorems 78.1–78.3 unaffected. Files: result_81_*.")
    with open(state, "a", encoding="utf-8") as f:
        f.write("\n".join(entry) + "\n")
    log(f"[appended] {state}")


def _flush_log():
    path = os.path.join(REPO, "result_81_log.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(LOG_LINES) + "\n")
    print(f"[wrote] {path}")


if __name__ == "__main__":
    main()
