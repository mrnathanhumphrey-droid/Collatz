"""
PROBE D2-e -- THE LADDER MATRICES (raw material; sigma-judge SEALED).
From the banked dense M_tower at L=2,3: gauge-Fourier (I (x) F over the gauge circle a=dlog),
reduce the block (e_rho,gamma) structure by the E-form (k=0) QSD, and read the effective
gauge-frequency operator A[kout,kin]. Then:
  (1) the per-ladder effective matrices (L x L) for the k=+-1 and k=+-2 ladders (rungs k0*3^j mod D);
  (2) their eigenvalues vs the banked modes (doublet, m=2 seat, k=+-2 accounting 0.267/0.244);
      + the full block-resolved ladder subspace spectrum (leakage localization);
  (3) the empirical coupling structure vs the Decimation Lemma.
      PRE-REGISTERED (SHAPE, committed before looking):
        - delta on the ki = 3*ko (sub/superdiagonal ladder step) -- the x3 decimation delta;
        - a TAIL COLUMN feeding the BOTTOM rung from all frequencies, Dirichlet ~1/distance profile;
        - CORNER-DOMINANCE (cyclic companion, single effective corner) ABSENT.
Internal exactness check: A[0,0] = rho(E-form) = 1/3 (weight-free, D0.3).
INSTRUMENT LAW: dense/direct at q=3. No fit. Exact rationals where feasible (E-form Perron).
"""
import numpy as np
from collections import OrderedDict
from fractions import Fraction as Fr
from probe_phase2c0 import build_M_tower_and_coords

np.set_printoptions(linewidth=160, suppress=True)

# banked modes (from result_trackD2a.md C3 / D2-d), for the spectrum comparison
BANKED = {
    2: {"doublet": [], "note": "L=2: single top pair"},
    3: {
        "doublet": [0.237639959367 + 0.183030417014j, 0.234998609841 + 0.183154982890j],  # m=1 (k=+-1) top two
        "m2_seat": [0.02024 + 0.18363j, 0.00406 + 0.19035j],                                # m=2 (k=+-2) occupants
        "k4_seat": [-0.009554 + 0.002875j],                                                 # k=+-4
    },
}


def ladder_rungs(k0, L, D):
    """rung gauge-frequencies k0*3^j mod D, j=0..L-1 (top coprime rung -> deeper internal rungs -> DC/Nyquist)."""
    return [(k0 * 3 ** j) % D for j in range(L)]


def build_gauge_fourier(L, lam=0.5):
    """Return M-tilde[bo,ko,bi,ki] (gauge-Fourier of M_tower over the a-circle), block list, D."""
    Mt, states, idx, tw, pos, twcoords, q, qL, sub, D, dl, w, two = build_M_tower_and_coords(L, lam)
    Md = Mt.toarray()
    by_eg = OrderedDict()
    for ti, (a, e, g) in enumerate(twcoords):
        by_eg.setdefault((e, g), {})[dl[a]] = ti
    blocks = list(by_eg.keys())
    Nb = len(blocks)
    assert all(len(m) == D for m in by_eg.values()), "block not full-D -- tensor structure broken"
    order = np.array([by_eg[b][al] for b in blocks for al in range(D)])
    Mperm = Md[np.ix_(order, order)]
    M4 = Mperm.reshape(Nb, D, Nb, D)                              # [bo,jo,bi,ji]
    F = np.exp(-2j * np.pi * np.outer(np.arange(D), np.arange(D)) / D) / np.sqrt(D)   # unitary DFT
    Fc = F.conj()
    Mtil = np.einsum('kj,ojpJ,lJ->okpl', F, M4, Fc)              # [bo,ko,bi,ki]
    return Mtil, blocks, Nb, D, Mt, twcoords, dl


def eform_qsd(Mtil, Nb):
    """E-form (k=0 sector) block operator = M-tilde[:,0,:,0]; return its left/right Perron (block QSD)."""
    B = Mtil[:, 0, :, 0]
    assert np.max(np.abs(B.imag)) < 1e-11, "E-form not real"
    B = B.real
    evR, VR = np.linalg.eig(B)
    j = int(np.argmax(evR.real))
    rho = evR[j].real
    R = VR[:, j].real
    evL, VL = np.linalg.eig(B.T)
    jL = int(np.argmax(evL.real))
    Lp = VL[:, jL].real
    Lp = Lp / (Lp @ R)                                           # normalize Lp.R = 1
    return B, rho, Lp, R


def reduce_ladder(Mtil, Lp, R, rungs):
    """QSD-reduced effective ladder matrix A[a,b] over rung indices, a,b enumerate `rungs`."""
    n = len(rungs)
    A = np.zeros((n, n), dtype=complex)
    for a, ko in enumerate(rungs):
        for b, ki in enumerate(rungs):
            A[a, b] = Lp @ Mtil[:, ko, :, ki] @ R
    return A


def coupling_map(Mtil, D):
    """C[ko,ki] = Frobenius norm over blocks of the (ko,ki) sub-block -- which frequencies couple."""
    return np.sqrt(np.sum(np.abs(Mtil) ** 2, axis=(0, 2)))       # D x D


def block_ladder_spectrum(Mtil, Nb, rungs):
    """Full block-resolved spectrum on the ladder subspace (blocks x rungs) -- for leakage localization."""
    n = len(rungs)
    big = np.zeros((Nb * n, Nb * n), dtype=complex)
    for a, ko in enumerate(rungs):
        for b, ki in enumerate(rungs):
            big[a * Nb:(a + 1) * Nb, b * Nb:(b + 1) * Nb] = Mtil[:, ko, :, ki]
    ev = np.linalg.eig(big)[0]
    return sorted(ev, key=lambda z: -abs(z))


def nearest(spec, targets):
    out = []
    for t in targets:
        z = min(spec, key=lambda z: abs(z - t))
        out.append((t, z, abs(z - t)))
    return out


def run(L):
    print(f"\n{'='*82}\n## L={L}")
    Mtil, blocks, Nb, D, Mt, twcoords, dl = build_gauge_fourier(L)
    print(f"   tower dim = {Mt.shape[0]}  =  Nb(blocks) {Nb} x D(gauge) {D}")
    B, rho, Lp, R = eform_qsd(Mtil, Nb)
    print(f"   E-form (k=0) Perron rho(S) = {rho:.12f}   (D0.3 expects 1/3 = {1/3:.12f};  "
          f"dev {abs(rho-1/3):.2e}, {'EXACT-to-fp' if abs(rho-1/3)<1e-9 else 'DEV'})")

    # ---- coupling map (Decimation Lemma shape) ----
    C = coupling_map(Mtil, D)
    print(f"\n   [3] COUPLING MAP  C[kout,kin] = ||block sub-block||_F   (D x D, gauge freq)")
    print(f"       row=kout, col=kin;  x3 ladder step is kin -> ko=3*ki (i.e. C[3*ki mod D, ki] large)")
    # print compactly
    with np.printoptions(precision=4, suppress=True):
        print(np.array2string(C, prefix="       "))
    # test the pre-registered shape on the k=+-1 ladder bottom rung
    for k0 in (1, 2):
        rungs = ladder_rungs(k0, L, D)
        bottom = rungs[-1]
        col = C[:, rungs[0]]                                     # what the TOP coprime rung feeds
        # delta at ko = 3*top:
        step_to = (3 * rungs[0]) % D
        delta_str = f"C[3*{rungs[0]}={step_to}, {rungs[0]}]={C[step_to, rungs[0]]:.4f}"
        # tail feeding bottom rung: row `bottom`, off the ladder-step entries
        tailrow = C[bottom, :].copy()
        print(f"       k=+-{k0} ladder rungs {rungs}:  x3-step delta {delta_str};  "
              f"bottom-rung({bottom}) in-feed row max off-diag = {np.max(np.delete(tailrow, bottom)):.4f}")

    # ---- reduced ladder matrices + spectra ----
    print(f"\n   [1]/[2] PER-LADDER EFFECTIVE MATRICES (QSD-reduced, L x L) + spectra vs banked")
    results = {}
    for k0 in (1, 2):
        rungs = ladder_rungs(k0, L, D)
        A = reduce_ladder(Mtil, Lp, R, rungs)
        evA = sorted(np.linalg.eig(A)[0], key=lambda z: -abs(z))
        bspec = block_ladder_spectrum(Mtil, Nb, rungs)
        print(f"\n   --- k=+-{k0} ladder (rungs {rungs}) ---")
        print(f"       A (reduced {len(rungs)}x{len(rungs)}):")
        with np.printoptions(precision=5, suppress=True):
            for r in range(len(rungs)):
                print("        " + "  ".join(f"{A[r,c].real:+.5f}{A[r,c].imag:+.5f}j" for c in range(len(rungs))))
        print(f"       reduced spectrum:  " + ", ".join(f"{z.real:+.5f}{z.imag:+.5f}j(|{abs(z):.5f}|)" for z in evA))
        print(f"       block-resolved ladder-subspace top-6:  " +
              ", ".join(f"{z.real:+.4f}{z.imag:+.4f}j" for z in bspec[:6]))
        results[k0] = (A, evA, bspec)

    # ---- explicit comparison to banked modes: LOCATE each banked mode across BOTH ladders ----
    if L in BANKED and BANKED[L].get("doublet"):
        print(f"\n   >> LOCATING BANKED MODES across both ladders (reduced rung-diag + block-resolved subspace):")
        # per-ladder: reduced diagonal rung values (labeled by gauge freq) + block-resolved spectrum
        pools = {}
        for k0 in (1, 2):
            A, evA, bspec = results[k0]
            rungs = ladder_rungs(k0, L, 2 * 3 ** (L - 1))
            pools[k0] = {"reduced_rungs": [(rungs[i], A[i, i]) for i in range(len(rungs))], "block": bspec}
        named = [("doublet(mod~.30)", BANKED[L]["doublet"]),
                 ("m=2 seat(mod~.19)", BANKED[L]["m2_seat"]),
                 ("k=+-4(mod~.010)", BANKED[L]["k4_seat"])]
        for label, targs in named:
            for t in targs:
                best = None
                for k0 in (1, 2):
                    for kf, zr in pools[k0]["reduced_rungs"]:
                        d = abs(zr - t)
                        if best is None or d < best[0]:
                            best = (d, k0, kf, zr, "reduced-rung")
                    zb = min(pools[k0]["block"], key=lambda z: abs(z - t))
                    d = abs(zb - t)
                    if d < best[0]:
                        best = (d, k0, None, zb, "block-resolved")
                d, k0, kf, z, kind = best
                loc = f"k=+-{k0} {kind}" + (f" rung(gf={kf})" if kf is not None else "")
                print(f"      banked {t.real:+.5f}{t.imag:+.5f}j (|{abs(t):.4f}|)  ->  {loc}: "
                      f"{z.real:+.5f}{z.imag:+.5f}j (d={d:.4f})")
    return results


def dump_tsv(allres, path="outputs/ladder_matrices_q3.tsv"):
    with open(path, "w", encoding="utf-8") as f:
        f.write("# D2-e per-ladder QSD-reduced effective matrices (gauge-Fourier of M_tower), q=3\n")
        for L, res in allres.items():
            for k0, (A, evA, bspec) in res.items():
                rungs = ladder_rungs(k0, L, 2 * 3 ** (L - 1))
                f.write(f"\n## L={L}  k=+-{k0}  rungs={rungs}\n")
                f.write("A[row,col]\t" + "\t".join(f"col_k={r}" for r in rungs) + "\n")
                for i, r in enumerate(rungs):
                    f.write(f"row_k={r}\t" + "\t".join(f"{A[i,j].real:+.8f}{A[i,j].imag:+.8f}j" for j in range(len(rungs))) + "\n")
                f.write("reduced_eig\t" + "\t".join(f"{z.real:+.8f}{z.imag:+.8f}j" for z in evA) + "\n")
                f.write("block_resolved_top8\t" + "\t".join(f"{z.real:+.6f}{z.imag:+.6f}j" for z in bspec[:8]) + "\n")
    print(f"\n   [dump] wrote {path}")


def main():
    print("# PROBE D2-e -- THE LADDER MATRICES. gauge-Fourier reduction of banked M_tower L=2,3. Judge SEALED.")
    allres = {}
    for L in (2, 3):
        allres[L] = run(L)
    dump_tsv(allres)
    print("\n# NOTE: reduced = E-form-QSD collapse of blocks (mean-field); block-resolved = full ladder subspace.")
    print("#       gap(reduced vs block-resolved) = block-off-QSD leakage. sigma-vs-spectra judge STAYS SEALED.")


if __name__ == "__main__":
    main()
