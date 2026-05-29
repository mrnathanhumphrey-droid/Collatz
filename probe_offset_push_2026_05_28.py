"""Pin c_inf and A in k*(n) = log2(3) n + c_inf + A/n. Extend MC to n=26, M=5e7,
5-point parabola peak. Then fit and test closed-form candidates from the b-distribution
(b_j = sum of 2 Geom(2): mean 4, var 4; cumulative b_[1,j] mean 4j var 4j)."""
import numpy as np, json, time
A_MAX=90
def make_mm(N):
    sh=(int(N).bit_length()//2)+1
    s=np.uint64(sh); mask=np.uint64((1<<sh)-1); t=np.uint64(1<<sh); Nu=np.uint64(N)
    def mm(a,b):
        b1=b>>s; b0=b&mask
        return (((((a*b1)%Nu)*t)%Nu)+((a*b0)%Nu))%Nu
    return mm,Nu
def kstar(n,M,seed=0):
    N=3**n; inv2=pow(2,-1,N)
    pw=np.array([pow(inv2,a,N) for a in range(A_MAX+1)],dtype=np.uint64)
    mm,Nu=make_mm(N); three=np.uint64(3); one=np.uint64(1)
    rng=np.random.default_rng(seed); S=np.zeros(M,dtype=np.uint64)
    for _ in range(n):
        a=rng.geometric(0.5,size=M); np.clip(a,1,A_MAX,out=a); S=mm(pw[a],(one+three*S)%Nu)
    klist=list(range(max(6,int(1.3*n)-10),int(1.7*n)+6)); vv={}
    for k in klist:
        c=np.uint64(pow(2,k,N)); t=mm(c,S).astype(np.float64); ang=-2*np.pi*t/N
        vv[k]=float(np.hypot(np.cos(ang).mean(),np.sin(ang).mean()))
    ks=np.array(sorted(vv)); vs=np.array([vv[k] for k in ks]); i=int(np.argmax(vs))
    lo,hi=max(0,i-2),min(len(ks),i+3); x=ks[lo:hi].astype(float); y=vs[lo:hi]
    co=np.polyfit(x,y,2); kp=-co[1]/(2*co[0]) if co[0]<0 else float(ks[i])
    return float(kp), int(ks[i]), vv[ks[i]]

M=int(5e7); L=np.log2(3)
print(f"M={M:,}\n n   k*(parab)  k*/n    t(s)")
rows=[]
for n in range(16,27):
    t0=time.time(); kp,ki,v=kstar(n,M); dt=time.time()-t0
    rows.append({"n":n,"kstar":kp}); print(f"{n:>2} {kp:>9.3f} {kp/n:>6.3f} {dt:>6.1f}",flush=True)
n=np.array([r["n"] for r in rows]); k=np.array([r["kstar"] for r in rows])
# fit k = L n + c + A/n  (slope FIXED at log2 3) and free
X1=np.vstack([np.ones_like(n,float),1.0/n]).T
co1,*_=np.linalg.lstsq(X1,k-L*n,rcond=None); cinf,A=co1
res1=(k-L*n)-X1@co1
Xf=np.vstack([n,np.ones_like(n,float),1.0/n]).T
cof,*_=np.linalg.lstsq(Xf,k,rcond=None)
print(f"\nFIXED-slope log2(3): k* = log2(3) n + ({cinf:.3f}) + ({A:.2f})/n   max|resid|={np.abs(res1).max():.3f}")
print(f"FREE 3-term fit: slope={cof[0]:.4f}  c={cof[1]:.3f}  A={cof[2]:.2f}")
print(f"\nclosed-form candidates for c_inf={cinf:.3f}:")
for name,val in [("-(4+2log2 3)",-(4+2*L)),("-(4+log2 3)",-(4+L)),("-2(1+log2 3)",-2*(1+L)),
                 ("-(E[b]+2 log2 3)=-(4+3.17)",-(4+2*L)),("-(6+...)",None),
                 ("-log2(108)=-log2(4*27)",-np.log2(108)),("-log2(96)",-np.log2(96)),
                 ("-(2+3 log2 3)",-(2+3*L)),("-4 log2 3 + ...",None)]:
    if val is not None: print(f"   {name:>28} = {val:.3f}  (diff {val-cinf:+.3f})")
print(f"\nclosed-form candidates for A={A:.2f}: var(b_[1,j]) coeff related?")
for name,val in [("4*log2(3)^2*?",4*L*L),("2/0.83^2",2/0.83**2),("(2 log2 3)^2",(2*L)**2),
                 ("var-step 4 * ...",4),("16",16),("8 log2 3",8*L),("4/(2-log2 3)^2",4/(2-L)**2)]:
    print(f"   {name:>20} = {val:.3f}  (diff {val-A:+.3f})")
json.dump({"rows":rows,"cinf":float(cinf),"A":float(A),"free":cof.tolist()},
          open("experiments_output/probe_offset_push_2026_05_28.json","w"),indent=2)
print("saved.")
