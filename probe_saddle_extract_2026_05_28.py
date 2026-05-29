"""
Extract the saddle structure of |S_chi(n)(2^k)| from the EXACT 1-D transfer operator.
Goal: locate WHERE (which j, which m) the phase decoherence that fixes the peak lives,
so the c_inf boundary-layer integral can be built for the right region.

Diagnostics:
 (A) saddle path m_j* via forward*backward occupation |U_j(m) V_{j+1}(m)|.
 (B) phase magnitude |theta_j| along the saddle path.
 (C) decay localization: zero theta_j on a window, measure |S| recovery -> which j carry decay.
 (D) how the decay-carrying region moves as k varies around k* (what balances to set the peak).
"""
import numpy as np

def build_VU(n, k, a_max=70, buffer=90, zero_mask=None):
    """Return |S(k)| and (optionally) full V_j, U_j arrays. zero_mask: set of j with theta_j:=0."""
    m_lo=-(2*n+buffer); k_hi=int(1.9*n)+5; xs=np.arange(m_lo,k_hi+1); W=len(xs); base=-m_lo
    wa=np.array([2.0**-a for a in range(1,a_max+1)])
    if zero_mask is None: zero_mask=set()
    # precompute phases ph_j over m-grid
    def phase_row(j):
        if j in zero_mask: return np.ones(W,dtype=np.complex128)
        Mj=3**(n-j+1); ordj=2*3**(n-j); ph=np.empty(W,dtype=np.complex128)
        for idx in range(W):
            x=int(xs[idx]); ph[idx]=np.exp(-2j*np.pi*(pow(2,x%ordj,Mj)/Mj))
        return ph
    # backward V
    Vs=[None]*(n+2); V=np.ones(W,dtype=np.complex128); Vs[n+1]=V.copy()
    phs={}
    for j in range(n,0,-1):
        ph=phase_row(j); phs[j]=ph; G=ph*V; Vn=np.zeros(W,dtype=np.complex128)
        for a in range(1,a_max+1): Vn[a:]+=wa[a-1]*G[:W-a]
        V=Vn; Vs[j]=V.copy()
    S=V[k+base]
    # forward U: U_0=delta_k; U_j(m)=ph_j(m) sum_a wa U_{j-1}(m+a)
    U=np.zeros(W,dtype=np.complex128); U[k+base]=1.0; Us=[None]*(n+1); Us[0]=U.copy()
    for j in range(1,n+1):
        Un=np.zeros(W,dtype=np.complex128)
        for a in range(1,a_max+1): Un[:W-a]+=wa[a-1]*U[a:]   # U_{j-1}(m+a)
        Un*=phs[j]; U=Un; Us[j]=U.copy()
    return S, xs, Vs, Us, base

def saddle_path(n,k):
    S,xs,Vs,Us,base=build_VU(n,k)
    print(f"  n={n}, k={k}: |S|={abs(S):.6e}")
    print(f"  {'j':>4} {'m_j*(occ peak)':>14} {'<m>_occ':>9} {'spread':>7} {'|theta_j*|':>11}")
    path=[]
    for j in range(1,n+1):
        occ=Us[j]*Vs[j+1]   # complex occupation; total sums to S
        a=np.abs(occ); 
        if a.sum()==0: continue
        i=int(np.argmax(a)); mstar=int(xs[i])
        mean=float((xs*a).sum()/a.sum()); 
        spread=float(np.sqrt((a*(xs-mean)**2).sum()/a.sum()))
        Mj=3**(n-j+1); ordj=2*3**(n-j); th=pow(2,mstar%ordj,Mj)/Mj
        path.append((j,mstar,mean,spread,th))
        if j<=6 or j%max(1,n//12)==0 or j>n-3:
            print(f"  {j:>4} {mstar:>14} {mean:>9.2f} {spread:>7.2f} {th:>11.5f}")
    return path,abs(S)

def decay_localization(n,k,nbins=12):
    S0=abs(build_VU(n,k)[0])
    print(f"\n  decay localization (zero theta_j on window -> |S| recovery), |S_full|={S0:.3e}:")
    print(f"  {'window j':>14} {'|S_masked|':>12} {'recovery x':>11}")
    edges=np.linspace(1,n+1,nbins+1).astype(int)
    recs=[]
    for b in range(nbins):
        j1,j2=edges[b],edges[b+1]; mask=set(range(j1,j2))
        Sm=abs(build_VU(n,k,zero_mask=mask)[0]); rec=Sm/S0
        recs.append((j1,j2,rec)); print(f"  [{j1:>4},{j2:>4}) {Sm:>12.3e} {rec:>11.3f}")
    return recs

if __name__=="__main__":
    L=np.log2(3)
    n=72; k=round(L*n-6.86)  # near k*
    print("="*70); print(f"SADDLE PATH  (k* ~ {k})"); print("="*70)
    saddle_path(n,k)
    decay_localization(n,k)
    # (D) move k around k* +-3, see where the decay region & |S| go
    print("\n  k-scan around peak (which window dominates the decay?):")
    for dk in (-4,-2,0,2,4):
        kk=k+dk; S=abs(build_VU(n,kk)[0]); print(f"    k={kk} (kappa={kk-L*n:+.2f}): |S|={S:.4e}")
