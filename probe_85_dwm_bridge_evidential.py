"""
Probe 85 Phase 1 — the evidential DWM bridge at j=2 (design A).

Extends probe_82's step-3 moment construction to n=5,6; float pi (Phase-0 finding).
Evidential test: is the DWM j=2 operator phase, at the evidential level r=n-3, the
SAME character as R81's independently-certified F-hat chirp at that r? j=1 is kept
as native scaffolding (NOT tested; pre-reg A). n=5 -> r=2 (divide out Probe-84's
derived omega_3^l offset); n=6 -> r=3 (clean, certified Mahler regime).

VALIDATION GATE: at n=3 the NATIVE construction must reproduce the archived Syracuse
targets G1=0.10783, G2=0.6089 (transcription anchor) before any n>=5 result is trusted.

Guards: raw 2^{-v-v'} weights at integration (NO b_prior pre-averaging); no magnitude
filters; do not tune. THEOREM_C_745 + Thms 78.1-78.3 not at stake.
"""
import sys, os, cmath, math, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"C:\Collatz")
from bilinear_pair_operator import build_markov_rational

LOG=[]
def log(m=""):
    print(m); LOG.append(str(m))

def float_pi(K):
    Kf=np.array([[float(x) for x in row] for row in K]); n=Kf.shape[0]
    pi=np.full(n,1.0/n)
    for _ in range(80):
        nx=pi@Kf; nx/=nx.sum()
        if np.abs(nx-pi).sum()<1e-15: return nx
        pi=nx
    return pi

def fhat_chirp_phase(r, ell, eps):
    """R81 F-hat chirp phase e_{3^{r+1}}(c*4^u) as an array over the <4> orbit u=0..3^r-1.
    Returns dict {4^u mod 3^{r+1} : unit-modulus phase}. This is the INDEPENDENTLY
    certified (R81/R81b) side of the bridge, computed from the flat-measure chirp."""
    q=3**(r+1); c=(pow(2,eps,q)*pow((1+3**r)%q,ell,q))%q
    orbit={}
    x=1  # 4^0
    for u in range(3**r):
        orbit[x]=cmath.exp(2j*math.pi*((c*x)%q)/q)   # e_q(c * 4^u)
        x=(x*4)%q
    return q,c,orbit

def dwm_j2_chirp_phase(n):
    """DWM j=2 operator phase, reduced to the evidential modulus 3^{n-2}=3^{r+1}, r=n-3.
    j=2: x_j = 3^{2}*2^{-b}; e_{3^n}(9 * z) = e_{3^{n-2}}(z). The chirp multiplier on the
    2-adic orbit element w=2^{-b}*(2^{-v}-2^{-v'}) is the unit 2^{-b}; we read the pure
    chirp e_{3^{r+1}}(unit * <orbit element>) as a function of the orbit element."""
    N=3**n; qr=3**(n-2)  # = 3^{r+1}, r=n-3
    # pure chirp on the reduced modulus: e_{qr}(w) as function of w in (Z/qr)
    # the DWM operator multiplies by e_{qr}(-xi * 2^{-b} * (2^{-v}-2^{-v'})); the chirp
    # kernel (independent of the state index xi and the -sign) is e_{qr}( orbit_elt ).
    return qr, {w: cmath.exp(2j*math.pi*(w%qr)/qr) for w in range(qr)}

# ---- moment construction (native), parametrized by n; validated at n=3 ----
def moments_native(n, V_MAX):
    N=3**n; TPI=2j*math.pi/N
    K,coprime=build_markov_rational(n); pi_f=float_pi(K)
    dim=len(coprime); idx={r:i for i,r in enumerate(coprime)}
    inv2=pow(2,-1,N); powinv=[pow(inv2,v,N) for v in range(0,4*V_MAX+2)]
    idx1=idx[1]
    def Mt(v,vp,j,b):
        M=np.zeros((dim,dim),complex)
        if v==vp: return M
        xj=(pow(3,2*j-2,N)*pow(inv2,b,N))%N; pd=(powinv[v]-powinv[vp])%N
        for i,xi in enumerate(coprime):
            t=idx.get((xi*powinv[v+vp])%N,-1)
            if t>=0: M[i,t]+=cmath.exp(-TPI*xi*xj*pd)
        return M
    def Off(j,b):
        M=np.zeros((dim,dim),complex)
        for v in range(1,V_MAX+1):
            for vp in range(1,V_MAX+1):
                if v!=vp: M+=2.0**(-v-vp)*Mt(v,vp,j,b)
        return M
    Off1=Off(1,0)
    g1=0j; g2s=0j; g2t=0j; g2d=0j; g2v=0j
    for v1 in range(1,V_MAX+1):
        for vp1 in range(1,V_MAX+1):
            if v1==vp1: continue
            w1=2.0**(-v1-vp1); b1=v1+vp1; X1=Mt(v1,vp1,1,0)-Off1; Off2=Off(2,b1)
            for v2 in range(1,V_MAX+1):
                for vp2 in range(1,V_MAX+1):
                    if v2==vp2: continue
                    ww=w1*2.0**(-v2-vp2); X2=Mt(v2,vp2,2,b1)-Off2
                    P121=X1@X2@X1; g1+=ww*P121.sum()
                    P1212=(X1@X2)@(X1@X2)
                    g2s+=ww*P1212.sum(); g2t+=ww*np.einsum('i,ii->',pi_f,P1212)
                    g2d+=ww*P1212[idx1,idx1]; g2v+=ww*(pi_f@P1212@pi_f)
    return dict(G1=g1.real,G2_sum=g2s.real,G2_tr=g2t.real,G2_d1=g2d.real,G2_vac=g2v.real)

if __name__=="__main__":
    log("# PROBE 85 Phase 1 — evidential DWM bridge (design A: j=2 F-hat substitution)\n")
    V=16
    # ---- VALIDATION GATE ----
    log("## VALIDATION GATE — n=3 native vs archived Syracuse targets")
    m3=moments_native(3,V)
    tg={'G1':0.10783,'G2_sum':0.6089,'G2_tr':0.05357,'G2_d1':0.05742,'G2_vac':0.004775}
    ok=True
    for k,t in tg.items():
        r=m3[k]/t; good=abs(r-1)<0.02; ok=ok and good
        log(f"   {k:7s}: got {m3[k]:+.5e}  target {t:.4e}  ratio {r:.4f}  {'OK' if good else 'FAIL'}")
    log(f"   GATE: {'PASS' if ok else 'FAIL — construction wrong, aborting'}\n")
    if not ok: sys.exit("validation gate failed")

    # ---- EVIDENTIAL PHASE AGREEMENT: DWM j=2 chirp vs R81 F-hat chirp ----
    log("## EVIDENTIAL — does the DWM j=2 chirp == R81 F-hat chirp at r=n-3?")
    for n in (5,6):
        r=n-3
        qr_d,dwm=dwm_j2_chirp_phase(n)
        log(f"\n n={n} (r={r}, modulus 3^{r+1}={qr_d}):")
        for ell in (0,1,2):
            for eps in (0,1):
                qf,c,fh=fhat_chirp_phase(r,ell,eps)
                # compare the two chirps on the shared orbit <4> mod 3^{r+1}
                # F-hat chirp key = 4^u; DWM chirp evaluated at same orbit element * c-multiplier?
                # bridge claim: F-hat's e_q(c*4^u) == DWM's e_q(4^u') under the step-1 unit/sign map.
                # test: is {arg fh(4^u)} a relabeling of {arg dwm} that is CONSTANT-offset (global phase)?
                fu=np.array([cmath.phase(fh[k]) for k in sorted(fh)])
                # DWM pure chirp sampled on same <4> orbit points (keys of fh are 4^u mod qr)
                du=np.array([cmath.phase(dwm[k%qr_d]) for k in sorted(fh)])
                diff=(fu-du)
                # a genuine bridge => diff is constant (global phase; = omega_3^l offset at r=2)
                spread=float(np.std(np.angle(np.exp(1j*(diff-diff.mean())))))
                off=float(np.mean(diff))
                log(f"   (l={ell},e={eps}) c={c}: phase-diff spread={spread:.2e} rad, "
                    f"global-offset={off:+.4f} ({'CONSTANT->same chirp' if spread<1e-9 else 'STRUCTURED->differs'})")
    log("\n(Interpretation: spread~0 across a family => DWM j=2 chirp = F-hat chirp up to a")
    log(" global phase [the Probe-84 omega_3^l at r=2, ~0 at r=3]. Nonzero spread => not the")
    log(" same chirp -> the identification is species-level only, not this-object-level.)")
    with open(r"C:\Collatz\result_85_log.txt","w",encoding="utf-8") as f: f.write("\n".join(LOG)+"\n")
    print("[wrote] result_85_log.txt")
