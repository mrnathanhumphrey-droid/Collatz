"""
NULL TEST for the 'slow xi = +-2^k (small dlog)' claim.
Falsifier: is dlog-small concentration real, or is it just generic smoothness
(which would instead concentrate at ADDITIVE-small xi)? And is it beyond chance (shuffle)?
Also: honest k*(n) growth-law trend (is k*/n -> log2(3)=1.585 or stuck ~1.2?).
"""
import math, json, numpy as np
from scipy.stats import spearmanr
from probe_syrac_charfn_decay_2026_05_28 import syrac_offset_distribution
PI={13:"probe_self_similarity/pi_13_truncated.npz",14:"probe_self_similarity/pi_14_truncated.npz",
    15:"probe_self_similarity/pi_15_truncated.npz",16:"probe_self_similarity/pi_16_truncated.npz"}
def order2(N):
    m,v=1,2%N
    while v!=1: v=(v*2)%N; m+=1
    return m
def absmu_full(n):
    N=3**n
    if n<=12: P=syrac_offset_distribution(n)
    else: d=np.load(PI[n]); P=np.zeros(N); P[d['coprime']]=d['pi']
    return np.abs(np.fft.fft(P)),N
def full_dlog(N,ordr):
    dl=np.full(N,-1,dtype=np.int64); g=1
    for ell in range(ordr): dl[g]=ell; g=(g*2)%N
    return dl

rng=np.random.default_rng(0)
print("="*92)
print("NULL: does |muhat| rank with DLOG-distance (claim) or ADDITIVE-distance (generic smooth)?")
print("="*92)
print(f"{'n':>3} {'spearman vs -d_dlog':>20} {'spearman vs -d_add':>20} {'|mu|@xi=1':>10} {'max|mu|':>9} {'ratio':>7}")
results=[]
for n in range(8,14):                       # full dlog feasible n<=13
    am,N=absmu_full(n); ordr=order2(N)
    dl=full_dlog(N,ordr)
    xi=np.arange(N); nt=(xi%3!=0)&(dl>=0)
    vals=am[nt]; ells=dl[nt]; xis=xi[nt]
    d_dlog=np.minimum.reduce([ells,ordr-ells,np.abs(ells-ordr//2)])
    d_add=np.minimum(xis,N-xis)
    rho_dlog=spearmanr(vals,-d_dlog).statistic
    rho_add =spearmanr(vals,-d_add ).statistic
    mu_xi1=am[1]; mx=vals.max()
    print(f"{n:>3} {rho_dlog:>20.4f} {rho_add:>20.4f} {mu_xi1:>10.5f} {mx:>9.5f} {mu_xi1/mx:>7.3f}")
    # shuffle null on top-40 mean d_dlog
    order=np.argsort(vals)[::-1]; topd=d_dlog[order[:40]].mean()
    B=2000; nulls=np.empty(B)
    for b in range(B):
        idx=rng.choice(len(vals),40,replace=False); nulls[b]=d_dlog[idx].mean()
    z=(topd-nulls.mean())/nulls.std()
    pct=(nulls<topd).mean()
    results.append({"n":n,"spearman_dlog":float(rho_dlog),"spearman_add":float(rho_add),
                    "top40_mean_d_dlog":float(topd),"shuffle_mean":float(nulls.mean()),
                    "shuffle_z":float(z),"shuffle_pctile":float(pct),
                    "mu_xi1":float(mu_xi1),"max":float(mx)})
print("\nshuffle null (top-40 mean dlog-distance vs 2000 random-40 draws):")
for r in results:
    print(f"  n={r['n']:>2}: real top40 mean d_dlog={r['top40_mean_d_dlog']:.0f}  "
          f"shuffle mean={r['shuffle_mean']:.0f}  z={r['shuffle_z']:+.1f}  pctile={r['shuffle_pctile']:.3f}")

# also report a few additive-small xi magnitudes vs max (generic-smooth would peak here)
print("\n|muhat| at additively-small xi (generic-smooth bottleneck) vs max, n=12:")
am,N=absmu_full(12)
for xi in [1,2,4,5,7,8,10,11]:
    print(f"  xi={xi}: |mu|={am[xi]:.5f}  (max={am[ (np.arange(N)%3!=0) ].max():.5f})")

# k*(n) trend: argmax small-side dlog, n=6..16, via BSGS
def make_bsgs(N,ordr):
    m=int(math.isqrt(ordr))+1; baby={}; g=1
    for b in range(m): baby.setdefault(g,b); g=(g*2)%N
    return baby,pow(pow(2,-1,N),m,N),m
def dlog2(t,N,baby,factor,m,ordr):
    g=t%N
    for q in range(m+2):
        if g in baby: return (q*m+baby[g])%ordr
        g=(g*factor)%N
print("\nk*(n) GROWTH LAW (peak dlog; clean pow2-vs-3 mechanism predicts k*/n -> log2(3)=1.585):")
print(f"{'n':>3} {'argmax xi':>10} {'k*=small-side dlog':>18} {'k*/n':>7}")
kstar=[]
for n in range(6,17):
    am,N=absmu_full(n); ordr=order2(N); xi=np.arange(N); nt=xi%3!=0
    a=am.copy(); a[~nt]=-1; xistar=int(np.argmax(a))
    baby,factor,m=make_bsgs(N,ordr); ell=dlog2(xistar,N,baby,factor,m,ordr)
    ks=min(ell,ordr-ell,abs(ell-ordr//2))
    kstar.append((n,ks)); print(f"{n:>3} {xistar:>10} {ks:>18} {ks/n:>7.3f}")
ns=np.array([k[0] for k in kstar]); ks=np.array([k[1] for k in kstar])
slope=np.polyfit(ns,ks,1)[0]
print(f"\n  linear fit k* ~ {slope:.3f}*n  (vs log2(3)=1.585).  k*/n at n=16: {ks[-1]/ns[-1]:.3f}")
print(f"  trend of k*/n across n: {[round(k/n,2) for n,k in kstar]}")
json.dump({"null":results,"kstar":[{'n':int(n),'kstar':int(k)} for n,k in kstar],
           "kstar_slope":float(slope)},
          open("experiments_output/probe_resonant_null_2026_05_28.json","w"),indent=2)
print("\nsaved experiments_output/probe_resonant_null_2026_05_28.json")
