"""
PROBE P1LVL (Wilson) -- the PER-LEVEL amplitude-vs-phase gate (2026-07-26). Fixes P1: apply the bridge at each shell.

Tower telescoping: gamma_r(k) = 1 + Sum_{j<=r} A_j(k), A_j(k) = Sum_{3-nmid a}|rho-hat_j(a)|^2 e(ak/3^j) = the PRIMITIVE
shell at level j (rho-hat_r(3^{r-j}b)=rho-hat_j(b) by the tower). The bridge reaches EVERY shell AT ITS OWN LEVEL, where
3-nmid a and |tau_a|=sqrt(q) (gate-verified 7.5e-13). class mean = A_1 (level-1); cascade c_k = Sum_{j>=2}A_j(k), and
A_j(1) = d1^{(j)} S_j, so "cascade positive" == "d1 positive" (the arc's question, new coordinates).

P1 per-level: scramble arg(pi-hat_j)=arg(W-hat_j) at each level j, bridge at level j -> A'_j(k), accumulate
gamma' = 1 + Sum_j A'_j, cascade' = Sum_{j>=2} A'_j. Does the cascade survive (amplitude) or collapse (phase)?
GATE FIRST: (1) per-level bridge A_j == direct A_j; (2) 1 + Sum_{j=1}^R A_j == certified gamma_R (telescoping holds).

Reuses probe_p1.build_level/bridge (BRIDGE2 exact, factored+verified). Cheap.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level, bridge, Mclass


def scramble(What, Nn, mode, seed=0):
    mag = np.abs(What)
    if mode == 'zero':
        return mag.astype(complex)
    rng = np.random.default_rng(seed)
    th = rng.uniform(-np.pi, np.pi, Nn); th[0] = 0.0
    for s in range(1, Nn // 2 + 1):
        th[Nn - s] = -th[s]
    return mag * np.exp(1j * th)


def shell_A(L, What_vec, ks):
    """A_j(k) = sum_{3-nmid a} |rho-hat_j(a)|^2 e(ak/Nn), rho-hat from the bridge at THIS level."""
    rhat = bridge(L, What_vec)
    Nn = L['Nn']; a = np.arange(Nn); prim = a % 3 != 0
    P = np.abs(rhat[prim]) ** 2; ap = a[prim]
    return {k: float(np.sum(P * np.cos(2 * np.pi * ap * k / Nn))) for k in ks}


def shell_A_direct(L, ks):
    rd = np.conj(np.fft.fft(L['rho']))
    Nn = L['Nn']; a = np.arange(Nn); prim = a % 3 != 0
    P = np.abs(rd[prim]) ** 2; ap = a[prim]
    return {k: float(np.sum(P * np.cos(2 * np.pi * ap * k / Nn))) for k in ks}


def main():
    t0 = time.time()
    R = 7
    KS = [1, 2, 3, 4, 5, 6]
    print(f"# PROBE P1LVL -- per-level bridge shells A_j(k); cascade = Sum_{{j>=2}} A_j (n up to {R})\n")
    Ls = {j: build_level(j) for j in range(1, R + 1)}

    # --- GATE 1: per-level bridge == direct ---
    gate1 = 0.0
    Areal = {}
    for j in range(1, R + 1):
        ab = shell_A(Ls[j], Ls[j]['What'], KS)
        ad = shell_A_direct(Ls[j], KS)
        gate1 = max(gate1, max(abs(ab[k] - ad[k]) for k in KS))
        Areal[j] = ad
    print(f"## GATE1 per-level bridge A_j == direct A_j: max abs diff = {gate1:.2e} "
          f"[{'OK' if gate1 < 1e-9 else 'FAIL'}]")

    # --- GATE 2: telescoping 1 + Sum_j A_j == certified gamma_R ---
    gamma_tel = {k: 1.0 + sum(Areal[j][k] for j in range(1, R + 1)) for k in KS}
    rhoR = Ls[R]['rho']
    gamma_cert = {k: 3.0 ** R * float(np.dot(rhoR, np.roll(rhoR, -k))) for k in KS}
    gate2 = max(abs(gamma_tel[k] - gamma_cert[k]) for k in KS)
    print(f"## GATE2 telescoping 1+Sum_j A_j == 3^R<rho,shift rho>: max abs = {gate2:.2e} "
          f"[{'OK' if gate2 < 1e-9 else 'FAIL'}]\n")

    # --- per-shell A_j(1): the d1^{(j)} S_j sequence ---
    print("## per-shell A_j(1) (= d1^{(j)} S_j). class mean = 1 + A_1; cascade = Sum_{j>=2} A_j")
    print("   " + "  ".join(f"A_{j}(1)={Areal[j][1]:+.5f}" for j in range(1, R + 1)))
    casc_real = {k: sum(Areal[j][k] for j in range(2, R + 1)) for k in KS}
    print(f"   1+A_1(1) = {1+Areal[1][1]:.4f} (class mean 2/3={2/3:.4f})  cascade c_1 = {casc_real[1]:+.5f}\n")

    # --- SCRAMBLE per level ---
    for mode, lab in (('zero', 'ZERO-phase'), ('random', 'RANDOM-phase')):
        Ascr = {}
        neg = False
        for j in range(1, R + 1):
            Wp = scramble(Ls[j]['What'], Ls[j]['Nn'], mode, seed=j)
            Ascr[j] = shell_A(Ls[j], Wp, KS)
            rhop = np.fft.ifft(np.conj(bridge(Ls[j], Wp))).real
            neg = neg or (rhop.min() < -1e-9)
        casc_scr = {k: sum(Ascr[j][k] for j in range(2, R + 1)) for k in KS}
        print(f"## SCRAMBLE = {lab} (per level)  [inverse spectrum {'SIGNED' if neg else 'nonneg'}]")
        print(f"   per-shell A_j(1) scr: " + "  ".join(f"{Ascr[j][1]:+.4f}" for j in range(1, R + 1)))
        print(f"   {'k':>2} {'3|k':>4} {'c_k real':>10} {'c_k scr':>10} {'ratio':>7}")
        for k in KS:
            cr = casc_real[k]; cs = casc_scr[k]
            print(f"   {k:>2} {'yes' if k%3==0 else 'no':>4} {cr:>+10.5f} {cs:>+10.5f} "
                  f"{cs/cr if abs(cr)>1e-6 else float('nan'):>7.2f}")
        print()

    print("   [cascade survives ~measured => AMPLITUDE (phase claim dead); collapses to ~0 => PHASE-carried.")
    print("    dichotomy = sign(c_k) 3|k vs 3-nmid k under scramble. c_1 real should ~ +0.06.]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
