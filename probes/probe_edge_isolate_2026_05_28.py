"""Clean edge isolation (no half-split). mode 'top': theta=0 for m<0 -> pure TOP decay;
'tail': theta=0 for m>(n-j+1)log2 3 -> pure TAIL decay. Check TOP+TAIL ~ -log|S_full|."""
import numpy as np
L=np.log2(3); ln2=np.log(2)
def Sval(n,k,mode='full',a_max=70,buffer=90):
    m_lo=-(2*n+buffer); k_hi=int(1.9*n)+5; xs=np.arange(m_lo,k_hi+1); W=len(xs); base=-m_lo
    wa=np.array([2.0**-a for a in range(1,a_max+1)])
    V=np.ones(W,dtype=np.complex128)
    for j in range(n,0,-1):
        Mj=3**(n-j+1); ordj=2*3**(n-j); top=(n-j+1)*L
        ph=np.empty(W,dtype=np.complex128)
        for idx in range(W):
            x=int(xs[idx])
            if mode=='top' and x<0: ph[idx]=1.0; continue
            if mode=='tail' and x>top: ph[idx]=1.0; continue
            ph[idx]=np.exp(-2j*np.pi*(pow(2,x%ordj,Mj)/Mj))
        G=ph*V; Vn=np.zeros(W,dtype=np.complex128)
        for a in range(1,a_max+1): Vn[a:]+=wa[a-1]*G[:W-a]
        V=Vn
    return abs(V[k+base])
for n in [120,160]:
    print(f"\n=== n={n} ===")
    print(f"  {'kappa':>6} {'-logSfull':>10} {'TOP':>9} {'TAIL':>9} {'TOP+TAIL':>9} {'A_top':>9}")
    rows=[]
    for dk in range(-14,4):
        k=round(L*n)+dk; kap=k-L*n
        Sf=Sval(n,k,'full'); Tp=Sval(n,k,'top'); Tl=Sval(n,k,'tail')
        TOP=-np.log(Tp); TAIL=-np.log(Tl); FULL=-np.log(Sf)
        rows.append((kap,FULL,TOP,TAIL))
        Atop=TOP/2**kap if kap<-3 else float('nan')
        print(f"  {kap:>6.2f} {FULL:>10.4f} {TOP:>9.4f} {TAIL:>9.4f} {TOP+TAIL:>9.4f} {Atop:>9.4f}")
    R=np.array(rows); kap=R[:,0]
    deep=kap<=-8; A_top=np.median(R[deep,2]/2**kap[deep]) if deep.sum() else float('nan')
    r_slope=-np.polyfit(kap,R[:,3],1)[0]
    kstar_pred=np.log2(r_slope/(A_top*ln2)) if (A_top>0 and r_slope>0) else float('nan')
    kstar_meas=kap[np.argmin(R[:,1])]
    print(f"  additivity check: max|FULL-(TOP+TAIL)| = {np.max(np.abs(R[:,1]-R[:,2]-R[:,3])):.4f}")
    print(f"  A_top(deep)={A_top:.4f}  tail rate r={r_slope:.4f}  -> c_inf_pred={kstar_pred:.3f}  (measured peak {kstar_meas:.2f})")
