"""
Tighten k*(n) slope -> is it exactly log2(3)? Push n=16..24 (uint64 overflows past
n=20, so exact mulmod via hi/lo split, valid through n~26). Validate vs exact argmax
at n=16 (=20). Fit slope + FIXED-log2(3) residual-trend test (the real discriminator).
"""
import numpy as np, json, time
A_MAX=80

def make_mulmod(N):
    shift=(int(N).bit_length()//2)+1
    sh=np.uint64(shift); mask=np.uint64((1<<shift)-1); two_sh=np.uint64(1<<shift); Nu=np.uint64(N)
    def mm(a,b):  # (a*b)%N exact for a,b<N<2^(64-shift)
        b1=b>>sh; b0=b&mask
        hi=((a*b1)%Nu); hi=(hi*two_sh)%Nu
        lo=(a*b0)%Nu
        return (hi+lo)%Nu
    return mm,Nu

def mc_kscan(n,M,klist,seed=0):
    N=3**n; inv2=pow(2,-1,N)
    pw=np.array([pow(inv2,a,N) for a in range(A_MAX+1)],dtype=np.uint64)
    mm,Nu=make_mulmod(N); three=np.uint64(3); one=np.uint64(1)
    rng=np.random.default_rng(seed); S=np.zeros(M,dtype=np.uint64)
    for _ in range(n):
        a=rng.geometric(0.5,size=M); np.clip(a,1,A_MAX,out=a)
        S=mm(pw[a], (one+three*S)%Nu)
    out={}
    for k in klist:
        c=np.uint64(pow(2,k,N)); t=mm(c,S).astype(np.float64)
        ang=-2.0*np.pi*t/N
        out[k]=float(np.hypot(np.cos(ang).mean(),np.sin(ang).mean()))
    return out
def peak(kv):
    ks=np.array(sorted(kv)); vs=np.array([kv[k] for k in ks]); i=int(np.argmax(vs))
    if 0<i<len(ks)-1:
        y0,y1,y2=vs[i-1],vs[i],vs[i+1]; d=(y0-2*y1+y2)
        return ks[i]+(0.5*(y0-y2)/d if d!=0 else 0.0), int(ks[i])
    return float(ks[i]),int(ks[i])

M=int(3e7)
print(f"exact-mulmod MC, M={M:,}\n{'n':>3} {'argmax(int)':>11} {'k*(parab)':>10} {'|val|':>9} {'k*/n':>7} {'t(s)':>6}")
rows=[]
for n in range(16,25):
    t0=time.time(); klist=list(range(max(6,int(1.1*n)-10),int(1.7*n)+9))
    kv=mc_kscan(n,M,klist); ks,ki=peak(kv); dt=time.time()-t0
    rows.append({"n":n,"kstar":float(ks),"kint":ki,"kon":float(ks/n)})
    print(f"{n:>3} {ki:>11} {ks:>10.2f} {kv[ki]:>9.5f} {ks/n:>7.3f} {dt:>6.1f}",flush=True)

ns=np.array([r["n"] for r in rows]); ks=np.array([r["kstar"] for r in rows])
s,b=np.polyfit(ns,ks,1)
# slope SE
yhat=s*ns+b; resid=ks-yhat; sigma=resid.std(ddof=2); se=sigma/(np.sqrt(len(ns))*ns.std())
print(f"\nFREE fit n=16..24: slope={s:.3f} +- {se:.3f}  intercept={b:.2f}  (log2(3)=1.585)")
# fixed-slope log2(3): residual trend?
bf=np.mean(ks-1.585*ns); rf=ks-(1.585*ns+bf)
# trend of fixed-slope residual vs n (should be flat if true slope=1.585)
tr=np.polyfit(ns,rf,1)[0]
print(f"FIXED slope=1.585 fit: offset={bf:.2f}, residuals={[round(x,2) for x in rf]}")
print(f"   residual-vs-n trend slope={tr:+.4f}  (flat≈0 supports 1.585; strong neg => true slope<1.585)")
# fixed-slope 1.5 for contrast
b5=np.mean(ks-1.5*ns); r5=ks-(1.5*ns+b5); tr5=np.polyfit(ns,r5,1)[0]
print(f"   for contrast, FIXED slope=1.50 residual-trend slope={tr5:+.4f}")
print(f"\nk*/n: {[round(r['kon'],3) for r in rows]}")
json.dump({"rows":rows,"free_slope":float(s),"free_slope_se":float(se),
           "fixed1585_offset":float(bf),"fixed1585_resid_trend":float(tr)},
          open("experiments_output/probe_kstar_push2_2026_05_28.json","w"),indent=2)
print("saved.")
