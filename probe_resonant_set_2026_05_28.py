"""
Characterize the slow-decaying (resonant) frequencies of the Syracuse char fn.
Q1: are the slowest xi literal small powers of 2?
Q2: in discrete-log coordinates (xi = 2^ell mod 3^n, 2 is a primitive root),
    where do the large |mu_hat| sit?
Q3: how does the COUNT of slow frequencies grow with n? (sparsity of exceptional set)
"""
import math, json, numpy as np
from probe_syrac_charfn_decay_2026_05_28 import syrac_offset_distribution

PI={13:"probe_self_similarity/pi_13_truncated.npz",14:"probe_self_similarity/pi_14_truncated.npz",
    15:"probe_self_similarity/pi_15_truncated.npz",16:"probe_self_similarity/pi_16_truncated.npz"}

def get_absmu(n):
    N=3**n
    if n<=12: P=syrac_offset_distribution(n)
    else:
        d=np.load(PI[n]); P=np.zeros(N); P[d['coprime']]=d['pi']
    mu=np.fft.fft(P); xi=np.arange(N); nt=xi%3!=0
    return np.abs(mu), nt, N

def order2(N):
    m,v=1,2%N
    while v!=1: v=(v*2)%N; m+=1
    return m

def make_bsgs(N,ordr):
    m=int(math.isqrt(ordr))+1
    baby={}; g=1
    for b in range(m):
        baby.setdefault(g,b); g=(g*2)%N
    factor=pow(pow(2,-1,N),m,N)  # 2^{-m}
    return baby,factor,m
def dlog2(target,N,baby,factor,m,ordr):
    gamma=target%N
    for q in range(m+2):
        if gamma in baby: 
            x=q*m+baby[gamma]
            return x%ordr
        gamma=(gamma*factor)%N
    return None

def is_pow2_int(x):  # literal integer power of 2
    return x>0 and (x&(x-1))==0

rows=[]
print("="*100)
for n in range(8,17):
    absmu,nt,N=get_absmu(n)
    ordr=order2(N)
    a=absmu.copy(); a[~nt]=-1  # mask trivial
    mx=a.max()
    order=np.argsort(a)[::-1]
    topK=order[:40]
    baby,factor,m=make_bsgs(N,ordr)
    # classify top-K
    n_litpow2=0; n_neglitpow2=0; dlogs=[]; small_dlog=0
    rec=[]
    for xi in topK[:12]:
        xi=int(xi); ell=dlog2(xi,N,baby,factor,m,ordr)
        dlogs.append(ell)
        d_small=min(ell, ordr-ell)                 # distance to nearest +-2^small
        d_neg=min(abs(ell-ordr//2), ordr-abs(ell-ordr//2))  # distance to -1*2^small
        litp=is_pow2_int(xi); neglit=is_pow2_int(N-xi)
        if litp: n_litpow2+=1
        if neglit: n_neglitpow2+=1
        rec.append((xi,float(absmu[xi]),ell,d_small,d_neg,litp,neglit))
    # full-window dlog distance for ALL top-40
    for xi in topK:
        xi=int(xi); ell=dlog2(xi,N,baby,factor,m,ordr)
        if min(ell,ordr-ell,abs(ell-ordr//2))<=64: small_dlog+=1
        if is_pow2_int(xi): n_litpow2+=0
    # sparsity: count xi with |mu|>=tau*max for several tau
    cnt={}
    for tau in (0.9,0.7,0.5):
        cnt[tau]=int(np.sum(a>=tau*mx))
    print(f"n={n}: max={mx:.5f}  ord=2*3^{n-1}={ordr}")
    print(f"  top-12 (xi, |mu|, dlog ell, dist-to-+2^k, dist-to- -2^k, litPow2, -litPow2):")
    for xi,am,ell,ds,dn,lp,nl in rec:
        tag = "<2^k" if lp else ("<-2^k" if nl else "")
        print(f"    xi={xi:>9}  |mu|={am:.5f}  ell={ell:>8}  d+={ds:>7}  d-={dn:>7}  {tag}")
    print(f"  of top-40: {small_dlog} have min(ell, ord-ell, |ell-ord/2|) <= 64  (near +-2^small)")
    print(f"  count |mu|>=tau*max:  tau=0.9:{cnt[0.9]}  0.7:{cnt[0.7]}  0.5:{cnt[0.5]}")
    rows.append({"n":n,"max":float(mx),"ord":ordr,"top12_dlog":dlogs,
                 "near_pow2_of_top40":small_dlog,"count_ge":{str(k):v for k,v in cnt.items()}})
    print()

# plot |mu| vs dlog for n=12
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    n=12; N=3**n; ordr=order2(N)
    P=syrac_offset_distribution(n); absmu=np.abs(np.fft.fft(P))
    dl=np.full(N,-1,dtype=np.int64); g=1
    for ell in range(ordr):
        dl[g]=ell; g=(g*2)%N
    xi=np.arange(N); nt=(xi%3!=0)&(dl>=0)
    ells=dl[nt]; vals=absmu[nt]
    # signed distance to nearest of {small +2^k, small -2^k}
    dpos=np.minimum(ells,ordr-ells); dneg=np.abs(ells-ordr//2)
    dist=np.minimum(dpos,dneg)
    fig,ax=plt.subplots(1,2,figsize=(15,5.5))
    ax[0].scatter(ells,vals,s=2,alpha=.3)
    ax[0].set_xlabel("dlog_2(xi)  (ell, with -1 at ord/2)"); ax[0].set_ylabel("|mu_hat(xi)|")
    ax[0].set_title(f"n=12: char fn vs discrete log (ord={ordr})"); ax[0].grid(alpha=.3)
    ax[1].loglog(dist+1,vals,".",ms=2,alpha=.3)
    ax[1].set_xlabel("distance of dlog to nearest +-2^small  (+1)")
    ax[1].set_ylabel("|mu_hat|"); ax[1].set_title("n=12: large |mu| concentrate at small dlog-distance")
    ax[1].grid(alpha=.3,which="both")
    fig.tight_layout(); fig.savefig("probe_resonant_set_2026_05_28.png",dpi=110)
    print("plot saved: probe_resonant_set_2026_05_28.png")
except Exception as e:
    print("plot skip:",e)

json.dump(rows,open("experiments_output/probe_resonant_set_2026_05_28.json","w"),indent=2)
print("json saved.")
