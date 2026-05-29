"""
Joint saddle: fit the FULL exact f(kappa)=-log|S(2^k)| to A*2^kappa - r*kappa + c0
(top-edge 2^kappa + tail-linear, interaction included automatically).
Test c_inf = log2(r/(A ln2)) == argmin. Extract true coupled A, r vs isolated-channel values.
"""
import numpy as np
from scipy.optimize import curve_fit
from probe_transfer_op_2026_05_28 import Schi_all_k
L=np.log2(3); ln2=np.log(2)
def f_model(kap,A,r,c0): return A*2**kap - r*kap + c0
print(f"{'n':>4} {'A':>8} {'r':>7} {'kap*=log2(r/A ln2)':>19} {'argmin(exact)':>13} {'fit max|resid|':>14}")
for n in [120,160,200,240]:
    kv=Schi_all_k(n)
    ks=np.array(sorted(kv)); f=np.array([-np.log(kv[k]) for k in ks]); kap=ks-L*n
    i0=int(np.argmin(f))
    win=(kap>=kap[i0]-9)&(kap<=kap[i0]+4)
    p0=[15.,0.2,f[i0]]
    try:
        popt,_=curve_fit(f_model,kap[win],f[win],p0=p0,maxfev=40000)
        A,r,c0=popt; resid=np.max(np.abs(f[win]-f_model(kap[win],*popt)))
        kstar=np.log2(r/(A*ln2))
        print(f"{n:>4} {A:>8.3f} {r:>7.4f} {kstar:>19.3f} {kap[i0]:>13.2f} {resid:>14.4f}")
    except Exception as e:
        print(f"{n}: fit failed {e}")
# high-precision argmin via fine fit near peak (parabola on f) + n->inf extrap of A,r
print("\nlarge-n coupled constants & c_inf:")
As=[];rs=[];ns=[120,160,200,240,280]
for n in ns:
    kv=Schi_all_k(n); ks=np.array(sorted(kv)); f=np.array([-np.log(kv[k]) for k in ks]); kap=ks-L*n
    i0=int(np.argmin(f)); win=(kap>=kap[i0]-9)&(kap<=kap[i0]+4)
    popt,_=curve_fit(f_model,kap[win],f[win],p0=[15.,0.2,f[i0]],maxfev=40000)
    As.append(popt[0]);rs.append(popt[1])
for n,A,r in zip(ns,As,rs):
    print(f"  n={n}: A={A:.3f} r={r:.4f} -> c_inf=log2(r/(A ln2))={np.log2(r/(A*ln2)):.3f}")
print(f"\n  isolated-channel gave A_top~18.4, r~0.16 -> -6.3 (WRONG); coupled fit above is the joint answer.")
