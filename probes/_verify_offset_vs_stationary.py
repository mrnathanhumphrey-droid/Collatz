import numpy as np, scipy.sparse as sp
from probe_syrac_charfn_decay_2026_05_28 import syrac_offset_distribution

EPS_KNOWN = {1:+2.0e-01,2:+9.5238095238e-03,3:-5.0919863259e-03,4:-2.4522582483e-03,
 5:-1.1517469151e-03,6:-4.9790566522e-04,7:-1.1752368304e-03,8:-7.4554636729e-04,
 9:-7.5202571564e-06,10:+7.2075091711e-04,11:+1.5019670121e-03}

def order_of_two(N):
    m,v=1,2%N
    while v!=1: v=(v*2)%N; m+=1
    return m
def build_K(k,v_max=100):
    N=3**k; M=order_of_two(N); ve=min(M,v_max); inv2=pow(2,-1,N)
    pw=[]; p=inv2
    for v in range(ve): pw.append(p); p=(p*inv2)%N
    mask=np.ones(N,bool); mask[::3]=False; idxN=np.where(mask)[0]; n=len(idxN)
    sidx=-np.ones(N,np.int64); sidx[idxN]=np.arange(n)
    w=np.array([2.0**-(v+1) for v in range(ve)]); w/=w.sum()
    rows=[];cols=[];vals=[]; base=(3*idxN+1)%N
    for v in range(ve):
        t=(base*pw[v])%N; rows.append(np.arange(n)); cols.append(sidx[t]); vals.append(np.full(n,w[v]))
    K=sp.csr_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n))
    K.sum_duplicates(); return K,idxN

print("n | eps(offset L2-7/15)  eps_KNOWN(stationary)   diff   | ||K^T P_X - P_X||_1 (stationarity test)")
for n in range(1,9):
    N=3**n
    P=syrac_offset_distribution(n)
    absmu=np.abs(np.fft.fft(P)); xi=np.arange(N); nt=xi%3!=0
    eps_off=float(np.sum(absmu[nt]**2))-7.0/15.0
    K,idxN=build_K(n)
    Pu=P[idxN]               # offset dist restricted to units
    resid=float(np.linalg.norm(K.T@Pu - Pu,1))
    # also compare offset dist to true stationary of K
    pi=np.full(len(idxN),1/len(idxN)); KT=K.T.tocsr()
    for _ in range(2000):
        pn=KT@pi; pn/=pn.sum()
        if np.linalg.norm(pn-pi,1)<1e-15: pi=pn; break
        pi=pn
    dist_diff=float(np.linalg.norm(Pu-pi,1))
    ek=EPS_KNOWN[n]
    print(f"{n} | {eps_off:+.8e}  {ek:+.8e}  {eps_off-ek:+.2e} | stat_resid={resid:.2e}  ||P_X - pi_stationary||_1={dist_diff:.2e}")
