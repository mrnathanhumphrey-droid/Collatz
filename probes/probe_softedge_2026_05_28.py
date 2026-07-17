"""
Probe the soft-edge kernel g(d)=sum_{a>=1} 2^{-a} exp(-2pi i 2^{-(d+a)}) and its integral,
and EMPIRICALLY decompose log|S(kappa)| (early vs late phases) to find what sets the peak kappa*.
"""
import numpy as np
from scipy.integrate import quad
from probe_saddle_extract_2026_05_28 import build_VU
L=np.log2(3)

# ---- (1) the soft-edge kernel g(d) ----
def g(d, amax=200):
    a=np.arange(1,amax+1); return np.sum(2.0**-a * np.exp(-2j*np.pi*2.0**-(d+a)))
print("soft-edge kernel g(d)=sum 2^-a exp(-2pi i 2^-(d+a)):")
print(f"  {'d':>5} {'|g(d)|':>9} {'-log|g|':>9}")
for d in [-1,0,0.5,1,2,3,4,6,8,10]:
    gd=abs(g(d)); print(f"  {d:>5} {gd:>9.5f} {(-np.log(gd) if gd>0 else 9.9):>9.5f}")
I_full=quad(lambda d:-np.log(abs(g(d))), -1, 30, limit=200)[0]
I_pos =quad(lambda d:-np.log(abs(g(d))), 0, 30, limit=200)[0]
print(f"  integral_-1^inf -log|g| dd = {I_full:.5f}")
print(f"  integral_0^inf  -log|g| dd = {I_pos:.5f}")
print(f"  1/sqrt3={1/np.sqrt(3):.5f}  |g(-inf-ish, d=-3)|={abs(g(-3)):.4f} (should be O(1) oscillatory)")

# ---- (2) empirical peak decomposition: early vs late phases ----
n=120
print(f"\nEMPIRICAL peak decomposition, n={n}: keep ALL / EARLY only (zero j>n/2) / LATE only (zero j<=n/2)")
print(f"  {'kappa':>6} {'k':>4} {'|S_all|':>11} {'|S_early|':>11} {'|S_late|':>11}")
half=n//2
peak=None;best=-1
for dk in range(-14,5):
    k=round(L*n)+dk
    S_all =abs(build_VU(n,k)[0])
    S_early=abs(build_VU(n,k,zero_mask=set(range(half+1,n+1)))[0])  # keep early, zero late
    S_late =abs(build_VU(n,k,zero_mask=set(range(1,half+1)))[0])    # keep late, zero early
    kap=k-L*n
    mark=""
    if S_all>best: best=S_all; peak=kap; 
    print(f"  {kap:>6.2f} {k:>4} {S_all:>11.4e} {S_early:>11.4e} {S_late:>11.4e}")
print(f"  -> peak at kappa={peak:.2f} (target c_inf~-6.86)")
