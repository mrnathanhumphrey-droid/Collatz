"""
EXACT 1-D transfer operator for |S_chi(n)(2^k)|.
State = m = k - s_j (integer). Phase factor depends only on (m, j):
  theta_j(m) = (2^m mod 3^{n-j+1}) / 3^{n-j+1},  factor e(-theta_j(m)).
Backward recursion V_j(m) = sum_{a>=1} 2^{-a} e(-theta_j(m-a)) V_{j+1}(m-a), V_{n+1}=1.
Then |S_chi(2^k)| = |V_1(k)|  -- ONE pass gives all k.
Exact (complex128, no MC noise). O(n*width). Scales to n~100.
"""
import numpy as np, json, time

def Schi_all_k(n, a_max=70, buffer=90):
    k_hi = int(1.9*n)+5
    m_lo = -(2*n + buffer)
    xs = np.arange(m_lo, k_hi+1)
    W = len(xs)
    wa = np.array([2.0**-a for a in range(1, a_max+1)])
    V = np.ones(W, dtype=np.complex128)
    for j in range(n, 0, -1):
        Mj = 3**(n-j+1); ordj = 2*3**(n-j)
        ph = np.empty(W, dtype=np.complex128)
        for idx in range(W):
            x = int(xs[idx])
            val = pow(2, x % ordj, Mj)
            ph[idx] = np.exp(-2j*np.pi*(val/Mj))
        G = ph*V
        Vn = np.zeros(W, dtype=np.complex128)
        for a in range(1, a_max+1):
            Vn[a:] += wa[a-1]*G[:W-a]
        V = Vn
    base = -m_lo
    return {k: abs(V[k+base]) for k in range(max(1,int(1.15*n)-6), k_hi-3)}

def peak(kv):
    ks=np.array(sorted(kv)); vs=np.array([kv[k] for k in ks]); i=int(np.argmax(vs))
    if 0<i<len(ks)-1:
        lo,hi=max(0,i-2),min(len(ks),i+3); x=ks[lo:hi].astype(float); y=vs[lo:hi]
        c=np.polyfit(x,y,2); 
        return (-c[1]/(2*c[0]) if c[0]<0 else float(ks[i])), float(vs[i]), int(ks[i])
    return float(ks[i]),float(vs[i]),int(ks[i])

# ---- VALIDATION vs FFT at n=12 ----
print("VALIDATION (1-D transfer op vs exact FFT, n=12):")
from probe_syrac_charfn_decay_2026_05_28 import syrac_offset_distribution
P=syrac_offset_distribution(12); mu=np.abs(np.fft.fft(P)); N=3**12
kv12=Schi_all_k(12)
maxerr=0.0
for k in list(kv12)[:8]:
    xi=pow(2,k,N); err=abs(kv12[k]-mu[xi]); maxerr=max(maxerr,err)
    print(f"  k={k:>2}: transfer={kv12[k]:.8f}  fft={mu[xi]:.8f}  |diff|={err:.2e}")
print(f"  max|diff| over checked k = {maxerr:.2e}  -> {'EXACT MATCH' if maxerr<1e-7 else 'MISMATCH!'}")

# ---- PUSH k*(n) to large n ----
L=np.log2(3)
ns=[16,20,24,28,32,40,48,56,64,72,80,90]
print(f"\nk*(n) via exact transfer op (no MC noise):\n n   k*        k*/n    |S|peak   t(s)")
rows=[]
for n in ns:
    t0=time.time(); kv=Schi_all_k(n); kp,vp,ki=peak(kv); dt=time.time()-t0
    rows.append({"n":n,"kstar":kp,"peak":vp}); print(f"{n:>3} {kp:>9.4f} {kp/n:>6.3f} {vp:>9.5f} {dt:>6.1f}",flush=True)

n=np.array([r["n"] for r in rows]); k=np.array([r["kstar"] for r in rows])
# local slopes at the high end -- does dk*/dn -> log2(3)?
print("\nlocal slope dk*/dn (consecutive):")
for i in range(1,len(n)):
    print(f"  n={n[i-1]}->{n[i]}: {(k[i]-k[i-1])/(n[i]-n[i-1]):.4f}")
# free fits on growing tails -- does slope -> 1.585?
print("\nfree constant-slope fit on tails:")
for lo in [16,32,48,64]:
    m=n>=lo; s,b=np.polyfit(n[m],k[m],1); print(f"  n>={lo}: slope={s:.4f}  intercept={b:.3f}")
# fixed log2(3) + c + A/n
X=np.vstack([np.ones_like(n,float),1.0/n]).T; co,*_=np.linalg.lstsq(X,k-L*n,rcond=None)
print(f"\nfixed-slope log2(3): c_inf={co[0]:.3f}  A={co[1]:.2f}  max|resid|={np.abs((k-L*n)-X@co).max():.3f}")
json.dump({"rows":rows,"validation_maxerr":maxerr},open("experiments_output/probe_transfer_op_2026_05_28.json","w"),indent=2)
print("saved.")

# ===== EXTENSION: push to large n, nail slope + offset (exact, noise-free) =====
if __name__=="__main__":
    print("\n"+"="*70+"\nLARGE-n PUSH (exact)\n"+"="*70)
    L=np.log2(3)
    ns2=[40,60,80,100,120,140,160,180,200,220,240]
    rows2=[]
    for n in ns2:
        kv=Schi_all_k(n); kp,vp,ki=peak(kv); rows2.append({"n":n,"kstar":kp,"peak":vp})
        print(f"  n={n:>3}: k*={kp:.4f}  k*/n={kp/n:.4f}  |S|peak={vp:.3e}",flush=True)
    n=np.array([r["n"] for r in rows2]); k=np.array([r["kstar"] for r in rows2])
    print("\nlocal slope dk*/dn (does it -> log2(3)=1.58496?):")
    for i in range(1,len(n)):
        print(f"  n={n[i-1]:>3}->{n[i]:>3}: {(k[i]-k[i-1])/(n[i]-n[i-1]):.5f}")
    print("\nFREE 3-param fit k* = s n + c + A/n on n>=lo:")
    for lo in [40,80,120,160]:
        m=n>=lo; X=np.vstack([n[m],np.ones(m.sum()),1.0/n[m]]).T
        co,*_=np.linalg.lstsq(X,k[m],rcond=None)
        print(f"  n>={lo}: slope={co[0]:.5f}  c={co[1]:.3f}  A={co[2]:.2f}   (log2 3={L:.5f})")
    print("\nFIXED slope=log2(3), fit c_inf + A/n on n>=lo:")
    for lo in [40,120,200]:
        m=n>=lo; X=np.vstack([np.ones(m.sum()),1.0/n[m]]).T
        co,*_=np.linalg.lstsq(X,(k-L*n)[m],rcond=None)
        res=np.abs((k-L*n)[m]-X@co).max()
        print(f"  n>={lo}: c_inf={co[0]:.4f}  A={co[1]:.3f}  max|resid|={res:.4f}")
    # test closed forms for c_inf, A on the cleanest (n>=120) fit
    m=n>=120; X=np.vstack([np.ones(m.sum()),1.0/n[m]]).T; co,*_=np.linalg.lstsq(X,(k-L*n)[m],rcond=None)
    cinf,A=co
    print(f"\nCLEAN (n>=120) under slope=log2(3): c_inf={cinf:.4f}, A={A:.3f}")
    print("  c_inf candidates:")
    for nm,v in [("-(4+2log2 3)",-(4+2*L)),("-log2(108)",-np.log2(108)),("-(2+3log2 3)",-(2+3*L)),
                 ("-log2(112)",-np.log2(112)),("-(6+log2 ...)",None),("-7",-7.0),("-(4+log2 24)",-(4+np.log2(24)))]:
        if v is not None: print(f"     {nm:>16} = {v:.4f}  (diff {v-cinf:+.4f})")
    print("  A candidates:")
    for nm,v in [("16",16),("4 E[b]=16",16),("(2log2 3)^2=10.05",(2*L)**2),("8",8),("2/(2-L)^2",2/(2-L)**2),
                 ("4/(2-L)^2",4/(2-L)**2),("24",24),("E[b]Var=16",16)]:
        print(f"     {nm:>16} = {v:.3f}  (diff {v-A:+.3f})")
    json.dump({"rows2":rows2},open("experiments_output/probe_transfer_op_large_2026_05_28.json","w"),indent=2)
    print("saved.")
