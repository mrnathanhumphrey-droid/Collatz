"""
Extend max_xi|mu_hat(xi)| decay to n=16 and decide its functional form.
n<=12: direct offset computation. n=13..16: load on-disk pi_k (= offset dist
by the verified identity X_n == stationary(K_n)), FFT, take max/mean over xi!=0 mod 3.
Fit three models on the tail and report which wins + extrapolation.
"""
import json, os, time
import numpy as np
from scipy.optimize import curve_fit
from probe_syrac_charfn_decay_2026_05_28 import syrac_offset_distribution, v3

def stats_from_full(P, N):
    mu = np.fft.fft(P)
    xi = np.arange(N); nt = xi % 3 != 0
    a = np.abs(mu[nt])
    xis = xi[nt]; im = int(np.argmax(a)); xistar = int(xis[im])
    return {"max": float(a.max()), "mean": float(a.mean()),
            "median": float(np.median(a)),
            "L2": float(np.sum(a**2)), "argmax_xi": xistar, "argmax_v3": v3(xistar)}

PI = {13:"probe_self_similarity/pi_13_truncated.npz",
      14:"probe_self_similarity/pi_14_truncated.npz",
      15:"probe_self_similarity/pi_15_truncated.npz",
      16:"probe_self_similarity/pi_16_truncated.npz"}

rows=[]
print(f"{'n':>3} {'max|mu|':>11} {'mean|mu|':>11} {'argmax xi':>10} {'src':>8} {'t(s)':>7}")
for n in range(1,17):
    t0=time.time(); N=3**n
    if n<=12:
        P=syrac_offset_distribution(n); src="direct"
    else:
        d=np.load(PI[n]); P=np.zeros(N); P[d['coprime']]=d['pi']; src="pi_npz"
    s=stats_from_full(P,N); s["n"]=n; s["t"]=time.time()-t0
    rows.append(s)
    print(f"{n:>3} {s['max']:>11.6f} {s['mean']:>11.6f} {s['argmax_xi']:>10} {src:>8} {s['t']:>7.2f}")

ns=np.array([r['n'] for r in rows],float)
mx=np.array([r['max'] for r in rows],float)
ln=np.log(mx)

def fit_window(lo):
    m=ns>=lo; x=ns[m]; y=ln[m]; K=len(y)
    # power: y = a - A ln x
    pw=np.polyfit(np.log(x),y,1); A=-pw[0]; sse_pw=np.sum((y-np.polyval(pw,np.log(x)))**2)
    # geometric: y = a - r x
    gm=np.polyfit(x,y,1); rho=np.exp(gm[0]); sse_gm=np.sum((y-np.polyval(gm,x))**2)
    # stretched: y = a - c x^b
    def f(x,a,c,b): return a-c*np.power(x,b)
    try:
        p,_=curve_fit(f,x,y,p0=[0.0,0.5,0.8],bounds=([-5,1e-6,0.05],[5,50,2.0]),maxfev=20000)
        sse_st=np.sum((y-f(x,*p))**2); st=(float(p[0]),float(p[1]),float(p[2]))
    except Exception as e:
        sse_st=np.nan; st=None
    def aic(sse,p_): return K*np.log(sse/K)+2*p_
    return {"lo":lo,"K":K,
            "power":{"A":float(A),"sse":float(sse_pw),"aic":float(aic(sse_pw,2))},
            "geom":{"rho":float(rho),"sse":float(sse_gm),"aic":float(aic(sse_gm,2))},
            "stretched":{"a_c_beta":st,"sse":float(sse_st),"aic":float(aic(sse_st,3)) if sse_st==sse_st else None}}

print("\n=== 3-way fits (lower AIC = better) ===")
fits={}
for lo in (2,4,8):
    fw=fit_window(lo); fits[lo]=fw
    st=fw['stretched']
    print(f"\nwindow n>={lo} (K={fw['K']}):")
    print(f"  power  n^-A : A={fw['power']['A']:.3f}   AIC={fw['power']['aic']:.2f}")
    print(f"  geom rho^n  : rho={fw['geom']['rho']:.4f}  AIC={fw['geom']['aic']:.2f}")
    if st['a_c_beta']:
        a,c,b=st['a_c_beta']
        print(f"  stretched   : exp(-{c:.3f}*n^{b:.3f})  AIC={st['aic']:.2f}   <-- beta")

print("\n=== local geometric ratio (plateau? -> geometric; ->1 -> subgeometric) ===")
for i in range(1,len(ns)):
    print(f"  n={int(ns[i]):>2}: ratio={mx[i]/mx[i-1]:.4f}   local A={-np.log(mx[i]/mx[i-1])/np.log(ns[i]/ns[i-1]):.3f}")

# extrapolation under stretched + geom to n=20,27
a,c,b=fits[8]['stretched']['a_c_beta']
rho=fits[8]['geom']['rho']; ga=np.polyfit(ns[ns>=8],ln[ns>=8],1)
print("\n=== extrapolation of max|mu| ===")
for nn in (18,20,27):
    print(f"  n={nn}: stretched={np.exp(a-c*nn**b):.3e}   geometric={np.exp(np.polyval(ga,nn)):.3e}")

json.dump({"rows":rows,"fits":fits},open("experiments_output/probe_syrac_decay_fit_2026_05_28.json","w"),indent=2)
print("\nsaved experiments_output/probe_syrac_decay_fit_2026_05_28.json")
