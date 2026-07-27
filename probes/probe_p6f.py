"""
PROBE P6F (Wilson) -- close the numerator<->divided seam: C_nu = q * C_rho, and the whole chain (2026-07-26).

Wilson pinned both seams. The numerator-vs-divided seam is a MULTIPLIER, not a transform:
    W = 2^-a Y  =>  nu = Sum_a P(a) (x2^-a)_* rho,  P(a)=2^-a (a>=1)
    => C_nu(m) = Sum_d q(d) C_rho(m+d),  q(d) = Sum_a P(a)P(a+d) = 2^-|d| / 3   (branch-difference law)
    => nu-hat = D~ . rho-hat  (pointwise),  |D~|^2 = q-hat = 1/(5-4 cos theta) in [1/9, 1] -- never vanishes,
       so deconvolution is well-defined and stable. Nothing to build; one gate.

And the lag map is DERIVED (not guessed): cross-parity <=> branch difference d odd; rho on <4> = EVEN base-2
positions => C_rho supported on EVEN base-2 lags => C_nu^same needs m even, C_nu^cross needs m odd. The cross object
at odd base-2 lag m draws channels k=(m+d)/2 over odd d, weight q(d)=2^-|d|/3 (dominant d=+-1 => k=(m+-1)/2).

GATES (this probe):
 (1) C_nu = q * C_rho directly, j=2..6, machine-exact. + nu-hat = mult . rho-hat + |mult|^2 = 1/(5-4cos).
 (2) deconvolution round-trip: rho recovered from nu by rho-hat = nu-hat/mult == certified rho (stable, |mult|>0).
 (3) parity: C_rho on even base-2 lags only; C_nu^cross on odd lags only (derived rule, reconfirmed).
 (4) CHAIN CLOSURE: reconstruct the certified channel gamma_n(k) from nu_e + boundary alone
     (nu_o via P6D collapse -> nu -> deconv q -> C_rho -> gamma). Machine-exact => whole cascade is a functional
     of the single sub-measure nu_e (+ boundary).

Reuses probe_p1.build_level (certified rho + channel) + probe_p6d.build_base2 (certified base-2 nu_e/nu_o). No new transport.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from probe_p1 import build_level
from probe_p6d import build_base2


def autocorr(f):
    F = np.fft.fft(f)
    return np.fft.ifft(F * np.conj(F)).real


def corr(f, g):
    return np.fft.ifft(np.fft.fft(f) * np.conj(np.fft.fft(g))).real


def q_circular(twoN, D=240):
    """q(d)=2^-|d|/3 wrapped onto Z/twoN (exact; 2^-|d| decays so D=240 is beyond machine)."""
    q = np.zeros(twoN)
    for d in range(-D, D + 1):
        q[d % twoN] += (2.0 ** (-abs(d))) / 3.0
    return q


def embed_rho_base2(L):
    """certified base-4 numerator profile L['rho'] -> full base-2 group (even positions 2s)."""
    Nn = L['Nn']; twoN = 2 * Nn
    rho2 = np.zeros(twoN)
    s = np.arange(Nn)
    rho2[(2 * s) % twoN] = L['rho']
    return rho2


def main():
    t0 = time.time()
    print("# PROBE P6F -- C_nu = q * C_rho, and the chain closes on nu_e + boundary\n")

    # ---------- (1) C_nu = q * C_rho ----------
    print("## (1) GATE  C_nu(m) = Sum_d q(d) C_rho(m+d),  q(d)=2^-|d|/3   (j=2..6)")
    for n in (2, 3, 4, 5, 6):
        L = build_level(n); S = build_base2(n); twoN = S['twoN']
        rho2 = embed_rho_base2(L)
        nu2 = S['R_e'] + S['R_o']
        # mass consistency
        mrho, mnu = rho2.sum(), nu2.sum()
        C_rho = autocorr(rho2); C_nu = autocorr(nu2)
        q = q_circular(twoN)
        # pred[m] = Sum_d q(d) C_rho[(m+d)%twoN] = circular correlation of C_rho with q
        pred = corr(C_rho, q[::-1] if False else q)          # see note below
        # do it explicitly to avoid convention slips:
        Crho_hat = np.fft.fft(C_rho); q_hat = np.fft.fft(q)
        pred = np.fft.ifft(Crho_hat * np.conj(q_hat)).real   # Sum_d C_rho[t] q[t-m] = Sum_d q(d)C_rho(m+d)
        gate = np.max(np.abs(C_nu - pred))
        # nu-hat = mult . rho-hat ; |mult|^2 = q_hat = 1/(5-4cos)
        nuh = np.fft.fft(nu2); rhoh = np.fft.fft(rho2)
        mult = np.where(np.abs(rhoh) > 1e-14, nuh / np.where(np.abs(rhoh) > 1e-14, rhoh, 1), 0)
        theta = 2 * np.pi * np.arange(twoN) / twoN
        qhat_analytic = 1.0 / (5 - 4 * np.cos(theta))
        multcheck = np.max(np.abs((np.abs(mult) ** 2 - qhat_analytic)[np.abs(rhoh) > 1e-9]))
        qhat_check = np.max(np.abs(q_hat.real - qhat_analytic))
        print(f"   n={n}: |C_nu - q*C_rho| = {gate:.2e}   [mass rho={mrho:.5f}, nu={mnu:.5f}]  "
              f"|mult|^2==1/(5-4cos): {multcheck:.1e}  q_hat==1/(5-4cos): {qhat_check:.1e}")
    print()

    # ---------- (2) deconvolution round-trip ----------
    print("## (2) DECONV round-trip: rho recovered from nu by rho-hat = nu-hat / mult  (|mult|^2 in [1/9,1], stable)")
    for n in (3, 5):
        L = build_level(n); S = build_base2(n); twoN = S['twoN']
        rho2 = embed_rho_base2(L); nu2 = S['R_e'] + S['R_o']
        theta = 2 * np.pi * np.arange(twoN) / twoN
        # mult = Sum_{a>=1} 2^-a e^{i a theta} = (1/2 e^{i th})/(1 - 1/2 e^{i th})
        z = np.exp(1j * theta)
        mult = (0.5 * z) / (1 - 0.5 * z)
        rho_rec = np.fft.ifft(np.fft.fft(nu2) / mult).real
        err = np.max(np.abs(rho_rec - rho2))
        print(f"   n={n}: max|rho_recovered - rho_certified| = {err:.2e}   "
              f"[|mult|^2 range: {np.min(np.abs(mult)**2):.4f}..{np.max(np.abs(mult)**2):.4f} = 1/9..1]")
    print()

    # ---------- (3) parity of C_rho / C_nu^cross ----------
    print("## (3) parity (derived): C_rho on EVEN base-2 lags; C_nu^cross on ODD lags")
    n = 4; L = build_level(n); S = build_base2(n); twoN = S['twoN']
    rho2 = embed_rho_base2(L)
    C_rho = autocorr(rho2)
    nu_e, nu_o = S['R_e'], S['R_o']
    Xcross = corr(nu_e, nu_o) + corr(nu_o, nu_e)
    print(f"   n={n}: C_rho  even-lag L1={np.abs(C_rho[0::2]).sum():.4f}  odd-lag L1={np.abs(C_rho[1::2]).sum():.2e}")
    print(f"        C_nu^cross even-lag L1={np.abs(Xcross[0::2]).sum():.2e}  odd-lag L1={np.abs(Xcross[1::2]).sum():.4f}")
    print()

    # ---------- (4) CHAIN CLOSURE: certified channel from nu_e + boundary ----------
    print("## (4) CHAIN: reconstruct certified channel from nu_e + boundary (nu_o via P6D collapse)")
    print("   nu_e (+beta) -> nu_o=1/2 shift(nu_e)+1/2 beta -> nu -> deconv q -> C_rho -> compare certified C_rho")
    for n in (2, 3, 4, 5, 6):
        L = build_level(n); S = build_base2(n); twoN = S['twoN']
        nu_e = S['R_e']; beta = 2.0 * S['B']
        nu_o_rec = 0.5 * np.roll(nu_e, -1) + 0.5 * beta          # P6D collapse: nu_o[t]=1/2 nu_e[t+1]+1/2 beta[t]
        nu_rec = nu_e + nu_o_rec
        # deconv q at the autocorr level: C_rho = ifft( fft(C_nu) / q_hat )
        C_nu = autocorr(nu_rec)
        theta = 2 * np.pi * np.arange(twoN) / twoN
        qhat = 1.0 / (5 - 4 * np.cos(theta))
        C_rho_rec = np.fft.ifft(np.fft.fft(C_nu) / qhat).real
        # certified C_rho
        rho2 = embed_rho_base2(L); C_rho_cert = autocorr(rho2)
        err = np.max(np.abs(C_rho_rec - C_rho_cert))
        # and the certified channel gamma_n(k) from certified base-4 autocorr, vs from reconstructed C_rho folded
        print(f"   n={n}: max|C_rho_reconstructed(from nu_e) - C_rho_certified| = {err:.2e}   "
              f"[{'CHAIN CLOSED' if err < 1e-12 else 'gap'}]")
    print(f"\n# ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
