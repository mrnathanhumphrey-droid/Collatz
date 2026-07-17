"""
R81b — Mahler/3-adic-analytic certification of the arg F-hat phase.

Standalone extension of R81 (does NOT import or touch result_81_* files).

Goal: R81 found arg F-hat(3a) has phase whose polynomial degree in b GROWS with r
(deg = floor(r/2)+2) and whose leading finite-difference has growing 3-adic depth
(v3 = 3,4,4,6 at r=3,4,5,6). The claim under test: this is a 3-adic ANALYTIC
function in the Mahler (binomial) basis, i.e. v3(c_k) grows with k (Amice) -> bounded
analytic norm, NOT an unbounded-degree pathology.

Method:
  F-hat(3a) = 3 * ghat(a),  ghat(a) = sum_{j=0}^{d-1} e_q(c*4^j) * e_d(-a j),
  q = 3^(r+1), d = 3^r.  ghat = length-d DFT of the chirp e_q(c*4^j)  ->  numpy FFT.
  |ghat| = sqrt(q) on support {a == a0 mod 3} (Thm 78.3 smoke check).
  Phase index: ghat(a)^2 / q is a root of unity (Thm cert); identify its index J2 in
  Z/2q exactly from float64 (residual << 1/2q). Work with s2(a)=J2(a) mod q in Z/q.
  Finite-difference s2 over b (a=a0+3b) -> Mahler coeffs c_k = Delta^k s2(0) mod q.
  Report: degree (highest k with Delta^k not ==0 mod q), and the v3(c_k) profile.

Self-validation: r=3..6 must reproduce agent's (deg, v3_leading) = (3,3),(4,4),(4,4),(5,6).
"""
import numpy as np

def v3(n, q):
    """3-adic valuation of integer n taken mod q (q=3^(r+1)); returns r+1 if n==0 mod q."""
    n %= q
    if n == 0:
        # count from q itself
        k = 0
        t = q
        while t % 3 == 0:
            t //= 3; k += 1
        return k  # = r+1  (n is 0 to full 3-adic precision available)
    k = 0
    while n % 3 == 0:
        n //= 3; k += 1
    return k

def analyze(r, ell, eps, verbose=False):
    q = 3**(r+1)
    d = 3**r
    a0 = 1 if eps == 0 else 2
    c = (pow(2, eps, q) * pow((1 + 3**r) % q, ell, q)) % q

    # chirp[j] = e_q(c * 4^j)
    pow4 = np.empty(d, dtype=np.int64)
    acc = 1
    for j in range(d):
        pow4[j] = acc
        acc = (acc * 4) % q
    m = (c * pow4) % q                      # integer residues in [0,q)
    chirp = np.exp(2j * np.pi * (m / q))

    ghat = np.fft.fft(chirp)                # ghat[a] = sum_j chirp[j] e^{-2pi i a j/d}

    # smoke: |ghat| = sqrt(q) on a==a0 mod 3, ~0 elsewhere
    mag = np.abs(ghat)
    supp_mask = (np.arange(d) % 3 == a0)
    on = mag[supp_mask]; off = mag[~supp_mask]
    smoke_ok = (abs(on.mean() - np.sqrt(q)) < 1e-6 * np.sqrt(q)
                and on.std() < 1e-6 * np.sqrt(q) and off.max() < 1e-6 * np.sqrt(q))

    # phase index of ghat^2/q in Z/2q  (root of unity)
    a_supp = np.arange(a0, d, 3)            # a = a0 + 3b, b = 0..3^(r-1)-1
    z = ghat[a_supp]**2 / q                 # unit modulus, = +/- zeta_q^s
    ang = np.angle(z) % (2*np.pi)
    J2f = ang * (2*q) / (2*np.pi)
    J2 = np.rint(J2f).astype(np.int64) % (2*q)
    resid = np.max(np.abs(z - np.exp(2j*np.pi*J2/(2*q))))   # certify integer index

    s2 = J2 % q                             # in Z/q  (drops the +/-1 sign sigma)

    # finite differences over b, mod q, representative in (-q/2, q/2]
    diffs = [s2.copy() % q]
    while len(diffs[-1]) > 1:
        dd = (np.diff(diffs[-1])) % q
        diffs.append(dd)
    # degree = highest k with Delta^k not identically 0 mod q
    degree = 0
    for k in range(len(diffs)):
        if np.any(diffs[k] % q != 0):
            degree = k
    # Mahler coeffs c_k = Delta^k s2(0), and their v3
    ck = [int(diffs[k][0]) % q for k in range(min(len(diffs), degree + 3))]
    v3_profile = [v3(ck_i, q) for ck_i in ck]
    v3_lead = v3(int(diffs[degree][0]), q) if degree < len(diffs) else None
    # is it a genuine degree-`degree` polynomial over Z/q? (Delta^{degree} constant)
    is_poly = bool(degree + 1 >= len(diffs) or np.all(diffs[degree] % q == diffs[degree][0] % q)) \
              and bool(degree + 1 >= len(diffs) or np.all(diffs[degree + 1] % q == 0))
    return dict(r=r, ell=ell, eps=eps, q=q, c=c, smoke_ok=smoke_ok, resid=float(resid),
                degree=degree, v3_lead=v3_lead, v3_profile=v3_profile, is_poly=is_poly)


if __name__ == "__main__":
    print("R81b: Mahler / 3-adic-analytic certification of arg F-hat phase")
    print("="*76)
    agent = {3:(3,3), 4:(4,4), 5:(4,4), 6:(5,6)}  # (degree, v3_leading) from R81 disposition
    print(f"{'r':>2} {'q=3^(r+1)':>10} {'deg':>4} {'floor(r/2)+2':>12} {'v3_lead':>7}"
          f" {'max_resid':>10} {'smoke':>6} {'agent(deg,v3)':>14} {'match':>6}")
    summary = {}
    for r in range(3, 10):
        # aggregate across the 6 (ell,eps) families; report the max degree / the leading v3
        res = [analyze(r, ell, eps) for ell in (0,1,2) for eps in (0,1)]
        deg = max(x['degree'] for x in res)
        # leading v3 among families achieving the max degree
        v3l = max(x['v3_lead'] for x in res if x['degree'] == deg)
        maxres = max(x['resid'] for x in res)
        smoke = all(x['smoke_ok'] for x in res)
        exp_deg = r//2 + 2
        ag = agent.get(r)
        match = "" if ag is None else ("YES" if (deg, v3l) == ag else "NO")
        summary[r] = (deg, v3l)
        print(f"{r:>2} {3**(r+1):>10} {deg:>4} {exp_deg:>12} {v3l:>7} {maxres:>10.2e}"
              f" {str(smoke):>6} {str(ag):>14} {match:>6}")

    print("\nPer-family v3(c_k) profiles (Amice signature — should grow with k):")
    for r in range(3, 10):
        x = analyze(r, 1, 0)  # representative family ell=1,eps=0
        print(f"  r={r}: deg={x['degree']:>2}  is_poly={x['is_poly']}  "
              f"v3(c_k) k=0..: {x['v3_profile']}")

    print("\nVERDICT INPUTS:")
    degs = [summary[r][0] for r in range(3,10)]
    v3s  = [summary[r][1] for r in range(3,10)]
    print(f"  degree r=3..9      : {degs}   (predicted floor(r/2)+2: {[r//2+2 for r in range(3,10)]})")
    print(f"  v3_leading r=3..9  : {v3s}")
    print(f"  v3_leading monotone non-decreasing after r=3? "
          f"{all(v3s[i] <= v3s[i+1] for i in range(len(v3s)-1))}")
