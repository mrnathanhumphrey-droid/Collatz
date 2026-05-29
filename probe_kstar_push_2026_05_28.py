"""
Pin k*(n): does k*/n climb toward log2(3)=1.585 or stall ~1.2?
Monte-Carlo |S_chi(n)(2^k)| = |E exp(-2pi i 2^k X_n/3^n)| via per-sample recursion
X = U(1+3X) (n steps, U=2^{-a}, a~Geom(2)). uint64 modmul exact for n<=20.
Reuse X samples across all k. Parabola-fit peak -> continuous k*.
Validate vs EXACT argmax at n=15,16 (known: 18, 20).
"""
import numpy as np, json, time
A_MAX=64
EXACT_KSTAR={15:18,16:20}  # small-side dlog of exact argmax (from exact FFT)

def mc_charfn_over_k(n, M, klist, seed=0):
    N=3**n; inv2=pow(2,-1,N)
    pw=np.array([pow(inv2,a,N) for a in range(A_MAX+1)],dtype=np.uint64)
    rng=np.random.default_rng(seed)
    S=np.zeros(M,dtype=np.uint64)
    three=np.uint64(3); Nu=np.uint64(N); one=np.uint64(1)
    for _ in range(n):
        a=rng.geometric(0.5,size=M); np.clip(a,1,A_MAX,out=a)
        u=pw[a]
        S=(u*(((one+three*S)%Nu)))%Nu
    # S = X_n samples. scan k
    Sf=S.astype(np.float64)
    out={}
    for k in klist:
        c=np.uint64(pow(2,k,N))
        t=((c*S)%Nu).astype(np.float64)
        ang=-2.0*np.pi*t/N
        re=np.cos(ang).mean(); im=np.sin(ang).mean()
        out[k]=float(np.hypot(re,im))
    return out,N

def peak_parabola(kv):
    ks=np.array(sorted(kv)); vs=np.array([kv[k] for k in ks])
    i=int(np.argmax(vs))
    if 0<i<len(ks)-1:
        y0,y1,y2=vs[i-1],vs[i],vs[i+1]
        denom=(y0-2*y1+y2)
        delta=0.5*(y0-y2)/denom if denom!=0 else 0.0
        return ks[i]+delta, vs[i], int(ks[i])
    return ks[i], vs[i], int(ks[i])

M=int(4e7)
print(f"Monte-Carlo k* push, M={M:,} samples/n\n")
print(f"{'n':>3} {'argmax k (int)':>14} {'k* (parabola)':>14} {'|Schi(2^k*)|':>13} {'k*/n':>7} {'t(s)':>6} {'check':>10}")
rows=[]
for n in range(15,21):
    t0=time.time()
    kc=max(5,int(1.0*n)); klist=list(range(max(4,kc-12),int(1.7*n)+8))
    kv,N=mc_charfn_over_k(n,M,klist)
    kstar,vstar,kint=peak_parabola(kv)
    dt=time.time()-t0
    chk=""
    if n in EXACT_KSTAR:
        chk=f"exact={EXACT_KSTAR[n]}"
    rows.append({"n":n,"kstar_int":kint,"kstar_parab":float(kstar),"val":float(vstar),
                 "kstar_over_n":float(kstar/n),"kv":{str(k):v for k,v in kv.items()}})
    print(f"{n:>3} {kint:>14} {kstar:>14.2f} {vstar:>13.5f} {kstar/n:>7.3f} {dt:>6.1f} {chk:>10}")

ks=np.array([r["kstar_parab"] for r in rows]); ns=np.array([r["n"] for r in rows])
print(f"\nk*/n trend (n=15..20): {[round(r['kstar_over_n'],3) for r in rows]}")
print(f"log2(3) = 1.585  | linear slope over n=15..20: {np.polyfit(ns,ks,1)[0]:.3f}")
# combine with exact lower-n for full picture
print("\nFULL k*/n picture (exact n=6..16 + MC n=17..20):")
exact_lown={6:6,7:8,8:9,9:10,10:12,11:13,12:14,13:16,14:17,15:18,16:20}
allk=[]
for n in range(6,17): allk.append((n,exact_lown[n],n>=15 and "exact"))
for r in rows:
    if r["n"]>=17: allk.append((r["n"],r["kstar_parab"],"MC"))
for n,k,_ in allk: print(f"  n={n:>2}  k*={k if isinstance(k,int) else round(k,2):>6}  k*/n={k/n:.3f}")
json.dump(rows,open("experiments_output/probe_kstar_push_2026_05_28.json","w"),indent=2)
print("\nsaved.")
