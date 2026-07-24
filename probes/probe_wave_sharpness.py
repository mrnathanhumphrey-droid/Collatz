"""
Recompute the waveform self-similarity rho-scan and report the SHARPNESS of the peak
(Wilson's one-line question: is the corr(rho) maximum sharp at rho=3.00, or flat over [2.7,3.3]?).

Read-only: reconstructs the exact/validated eps chain (r=1..17) and, for a grid of rho,
overlaps lobe-3 (eps>0, r>=9.01) onto the sign-flipped, r-rescaled lobe-2 (eps<0, r in [2.65,9.01]),
then reports Pearson corr of the shapes vs rho. No state mutation, no build_nu.
"""
import numpy as np

# eps r=1..12 (exact r<=8, float r=9..12) from probe_crossing EPS_F; extend to 17 via eps_{r+1}=eps_r+2 Lam_r
EPS = {1: 0.2, 2: 9.523809523809525e-3, 3: -5.091986325893010e-3, 4: -2.452258248318762e-3,
       5: -1.151746915130986e-3, 6: -4.979056652200001e-4, 7: -1.175236830400000e-3,
       8: -7.455463672900000e-4, 9: -7.520257156400000e-6, 10: 7.207509171100000e-4,
       11: 1.501967012082273e-3, 12: 2.274713720558208e-3}
LAM = {12: 3.3677e-4, 13: 3.1971e-4, 14: 2.8672e-4, 15: 2.6193e-4, 16: 2.3426e-4}  # nu-validated
for r in range(12, 17):
    EPS[r + 1] = EPS[r] + 2 * LAM[r]

rs = np.array(sorted(EPS))
ev = np.array([EPS[r] for r in rs])

# crossings (linear interp): +->- near 2.65, -->+ near 9.01
def cross(a, b):
    return a - EPS[a] * (b - a) / (EPS[b] - EPS[a])
c1 = cross(2, 3)     # ~2.65
c2 = cross(9, 10)    # ~9.01
print(f"# crossings: c1(+->-)={c1:.3f}  c2(-->+)={c2:.3f}")

# lobe2 = negative excursion on [c1, c2]; lobe3 = positive excursion on [c2, 17]
def seg(lo, hi):
    m = (rs >= lo - 1e-9) & (rs <= hi + 1e-9)
    return rs[m], ev[m]

r2, e2 = seg(c1, c2)   # lobe 2 (negative)
r3, e3 = seg(c2, 17)   # lobe 3 (positive, still rising)

# For a rho: map lobe3's r-support back by /rho into lobe2's frame, interpolate lobe2 there,
# compare -A*eps2(r3/rho) to eps3(r3). Use phase u in [0,1] across each lobe for a clean overlap.
def lobe_interp(rr, ee, u):
    # u in [0,1] across the lobe support
    x = rr[0] + u * (rr[-1] - rr[0])
    return np.interp(x, rr, ee)

U = np.linspace(0.05, 0.95, 40)   # avoid the transient-contaminated very-early lobe2 edge
e3u = lobe_interp(r3, e3, U)
print(f"\n# rho-scan: corr( eps3(u) , -eps2(u) ) is phase-shape only; the rho enters as the")
print(f"# ratio of lobe LENGTHS (len3/len2). Report corr and best-fit amplitude A vs rho.")
print(f"{'rho':>6} {'corr':>8} {'A(amp)':>9} {'resid%':>7}")
best = (None, -2)
lens = (r2[-1] - r2[0], r3[-1] - r3[0])
for rho in np.arange(2.4, 3.7 + 1e-9, 0.05):
    # a candidate rho predicts len3 = rho*len2 (lobe lengthening by rho); score how the
    # sign-flipped lobe2 shape matches lobe3 under that phase alignment, weighted toward
    # the region [c2, c2+rho*len2] we actually have data for.
    frac = min(1.0, lens[1] / (rho * lens[0]))   # fraction of predicted lobe3 we've observed
    Uobs = np.linspace(0.05, 0.05 + 0.90 * frac, 40)
    a = lobe_interp(r3, e3, Uobs)
    b = -lobe_interp(r2, e2, Uobs)               # sign-flipped lobe2 over same phase
    A = np.dot(a, b) / np.dot(b, b)
    resid = np.linalg.norm(a - A * b) / np.linalg.norm(a)
    cc = np.corrcoef(a, b)[0, 1]
    print(f"{rho:>6.2f} {cc:>8.4f} {A:>9.4f} {resid*100:>6.1f}%")
    if cc > best[1]:
        best = (rho, cc)
print(f"\n# best corr at rho={best[0]:.2f} (corr={best[1]:.4f})")
print(f"# SHARPNESS: compare corr at rho=3.00 vs rho in [2.7,3.3] above -- flat => rho imposed; peaked => measured")
