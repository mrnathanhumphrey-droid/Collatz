"""
Genuine LD calc for the top-edge amplitude.
Depth walk d_j: increments (a-log2 3), a~Geom(2), drift +0.415, start h=log2(3)-kappa.
Top-edge decoherence = ruin (d_j hits <=0). Cramer root gamma=ln2 (x=1/2 exact).
Compute ruin prob psi(h) ~ C 2^{-h}; compare to -log|S_early|(kappa) ~ A_top 2^kappa;
assemble peak kappa* = log2( r / (A_top ln2) ), r = (ln3/2)/v* (tail rate from 1/sqrt3 + v*).
"""
import numpy as np
from probe_saddle_extract_2026_05_28 import build_VU
L=np.log2(3); ln2=np.log(2); ln3=np.log(3)

# --- verify Cramer root ---
def cramer(g): 
    x=np.exp(-g); return np.exp(g*L)*x/(2-x)  # E[e^{-g(a-log2 3)}]
print(f"Cramer check: E[e^(-gamma(a-log2 3))] at gamma=ln2 = {cramer(ln2):.6f} (should=1). x=1/2 exact root.")

# --- ruin probability psi(h)=P(walk ever <=0 from h), increments a-log2 3 ---
H=45.0; dx=0.005; G=int(H/dx)+1; hs=np.arange(G)*dx
amax=80; wa=np.array([2.0**-a for a in range(1,amax+1)]); shifts=np.array([a-L for a in range(1,amax+1)])
psi=np.zeros(G)  # psi(h)=0 init for h>0
for it in range(4000):
    new=np.zeros(G)
    for w,sh in zip(wa,shifts):
        tgt=hs+sh                      # new height
        idx=np.clip(np.round(tgt/dx).astype(int),-1,G-1)
        val=np.where(tgt<=0,1.0,psi[idx])   # ruin if lands <=0
        new+=w*val
    if np.max(np.abs(new-psi))<1e-12: psi=new; break
    psi=new
# fit psi(h) ~ C 2^{-h} on moderate h
msk=(hs>=4)&(hs<=14); slope,b=np.polyfit(hs[msk],np.log(psi[msk]),1)
C=np.exp(b); rate=-slope
print(f"ruin psi(h) ~ C*exp(-rate*h): rate={rate:.4f} (ln2={ln2:.4f}), C={C:.4f}")
print(f"  => -log|S_early| predicted (if full decoherence) = psi(log2 3 - kappa) ~ (C/3) 2^kappa, C/3={C/3:.4f}")

# --- measure A_top from operator: -log|S_early| = A_top 2^kappa at deep kappa ---
n=160; half=n//2
print(f"\nmeasure -log|S_early| at deep kappa (n={n}):")
amps=[]
for dk in [-18,-16,-14,-12,-10]:
    k=round(L*n)+dk; Se=abs(build_VU(n,k,zero_mask=set(range(half+1,n+1)))[0])
    kap=k-L*n; mlog=-np.log(Se); A=mlog/2**kap
    amps.append(A); print(f"  kappa={kap:.2f}: -log|S_early|={mlog:.3e}  A_top=-log/2^kappa={A:.4f}")
A_top=np.median(amps)
print(f"  A_top (median) = {A_top:.4f}   [ruin pred C/3={C/3:.4f}; ratio A_top/(C/3)={A_top/(C/3):.3f} = decoherence/overshoot factor]")

# --- tail rate r and assemble kappa* ---
vstar=1.78; r=(ln3/2)/vstar
kstar=np.log2(r/(A_top*ln2))
print(f"\ntail rate r=(ln3/2)/v* = {r:.4f} (v*={vstar})")
print(f"PREDICTED c_inf = log2( r/(A_top ln2) ) = {kstar:.3f}    [measured c_inf ~ -6.86]")
# also if we used full-decoherence C/3 instead of measured A_top:
print(f"  (using ruin C/3 instead of A_top: c_inf_pred={np.log2(r/((C/3)*ln2)):.3f})")
