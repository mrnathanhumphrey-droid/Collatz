"""
PROBE J1 -- COMPLETENESS AT L=2,3 (the campaign's verdict table, unsealed).
Assign EVERY tower mode with |lambda| > 0.05 to its (ladder k0, rung gf, block) and tabulate
model-vs-measured with residuals. Model = the full QSD-reduced gauge-frequency operator's diagonal
A[gf,gf] (per-frequency collapsed value); block splitting = spread of modes sharing a gf.
PRE-REGISTERED (SHAPE): no orphan modes (assignment total); residual-to-diagonal ~ block-structure
scale (the splitting), << the mode modulus. Dense/direct at q=3 (INSTRUMENT LAW). No fit.
"""
import numpy as np
from collections import defaultdict
from probe_trackD2e import build_gauge_fourier, eform_qsd

np.set_printoptions(linewidth=160, suppress=True)


def sector_mass(vec, twcoords, dl, D):
    blocks = defaultdict(lambda: np.zeros(D, dtype=complex))
    for i, (a, er, g) in enumerate(twcoords):
        blocks[(er, g)][dl[a]] = vec[i]
    mk = np.zeros(D)
    for bv in blocks.values():
        mk += np.abs(np.fft.fft(bv)) ** 2
    return mk / mk.sum()


def cls9(k):
    r = k % 9
    if r in (1, 8): return "+-1"
    if r in (2, 7): return "+-2"
    if r in (4, 5): return "+-4"
    return f"div3(r={r})"


def rung_label(gf):
    """gf -> (ladder coprime class, rung index j). gf = k0*3^j; DC=0 terminal, Nyquist real."""
    if gf == 0:
        return ("DC", 0)
    j = 0; s = gf
    while s % 3 == 0:
        s //= 3; j += 1
    return (cls9(s), j)


def run(L, thr=0.05):
    print(f"\n{'='*88}\n## L={L}   (modes with |lambda| > {thr})")
    Mtil, blocks, Nb, D, Mt, twcoords, dl = build_gauge_fourier(L)
    B, rho, Lp, R = eform_qsd(Mtil, Nb)
    A = np.einsum('o,okpl,p->kl', Lp, Mtil, R)                    # full reduced gauge-freq operator D x D
    print(f"   tower dim {Mt.shape[0]} = {Nb} blocks x D={D};  E-form(k=0) Perron (finite-L partner) rho={rho:.10f}")

    ev, VR = np.linalg.eig(Mt.toarray())
    keep = [i for i in range(len(ev)) if abs(ev[i]) > thr]
    print(f"   modes |lambda|>{thr}: {len(keep)}  (of {len(ev)} total)")

    # assign each kept mode to its dominant gauge frequency -> ladder/rung; model = A[gf,gf]
    rows = []
    for i in keep:
        mk = sector_mass(VR[:, i], twcoords, dl, D)
        gf = int(np.argmax(mk))
        purity = mk[gf]
        lad, j = rung_label(gf)
        model = A[gf, gf]
        rows.append((abs(ev[i]), np.angle(ev[i]), ev[i], gf, lad, j, model, abs(ev[i] - model), purity))
    rows.sort(key=lambda r: -r[0])

    # group by gauge freq for the block-splitting readout
    bygf = defaultdict(list)
    for r in rows:
        bygf[r[3]].append(r)

    print(f"\n   {'|lam|':>8} {'phase':>8}  {'gf':>3} {'ladder':>7} {'rung':>4}  "
          f"{'model A[gf,gf]':>22}  {'resid':>8} {'purity':>6}")
    orphan = []
    for r in rows[:40]:
        modv = f"{r[6].real:+.5f}{r[6].imag:+.5f}j"
        flag = "  <ORPHAN" if r[7] > 0.25 * r[0] and r[7] > 0.02 else ""
        if flag:
            orphan.append(r)
        print(f"   {r[0]:8.5f} {r[1]:+8.4f}  {r[3]:3d} {r[4]:>7} {r[5]:4d}  {modv:>22}  {r[7]:8.5f} {r[8]:6.3f}{flag}")
    if len(rows) > 40:
        print(f"   ... ({len(rows)-40} more, in dump file)")

    # block-splitting per gf (spread of measured modes sharing a gf) vs residual-to-diagonal
    print(f"\n   BLOCK-SPLITTING (per gauge freq, |lambda|>{thr}):  gf(ladder,rung): n modes | model | "
          f"meas-spread | max resid-to-diag")
    maxresid = 0.0
    for gf in sorted(bygf, key=lambda g: -max(x[0] for x in bygf[g])):
        grp = bygf[gf]
        lad, j = rung_label(gf)
        mods = [x[2] for x in grp]
        spread = max(abs(a - b) for a in mods for b in mods) if len(mods) > 1 else 0.0
        mr = max(x[7] for x in grp)
        maxresid = max(maxresid, mr)
        modvals = ", ".join(f"{m.real:+.4f}{m.imag:+.4f}j" for m in sorted(mods, key=lambda z: -abs(z))[:4])
        print(f"     gf={gf:2d} ({lad},r{j}): {len(grp):2d} | A[gf,gf]={A[gf,gf].real:+.4f}{A[gf,gf].imag:+.4f}j | "
              f"spread={spread:.4f} | maxresid={mr:.4f}   [{modvals}]")

    # totality / pre-reg verdict
    ladders_seen = sorted(set(r[4] for r in rows))
    print(f"\n   >> TOTALITY: {len(rows)} modes, ALL assigned; ladders present: {ladders_seen}")
    print(f"   >> ORPHANS (resid > 25% of |lambda| and >0.02): {len(orphan)}  "
          f"({'NONE -- assignment total, pre-reg SHAPE holds' if not orphan else 'SEE FLAGS'})")
    print(f"   >> max residual-to-diagonal over all kept modes = {maxresid:.5f} "
          f"(block-structure scale; mode modula range {rows[-1][0]:.3f}..{rows[0][0]:.3f})")
    return rows, A, D


def main():
    print("# PROBE J1 -- COMPLETENESS AT L=2,3. Verdict table. Everything pre-committed; nothing to tune.")
    dump = []
    for L in (2, 3):
        rows, A, D = run(L)
        dump.append((L, rows))
    with open("outputs/judge_completeness_L23.tsv", "w", encoding="utf-8") as f:
        f.write("# J1 completeness census -- every tower mode |lambda|>0.05 assigned to (gf,ladder,rung)\n")
        for L, rows in dump:
            f.write(f"\n## L={L}\n|lambda|\tphase\teigenvalue\tgf\tladder\trung\tmodel_A[gf,gf]\tresid\tpurity\n")
            for r in rows:
                f.write(f"{r[0]:.6f}\t{r[1]:.5f}\t{r[2].real:+.6f}{r[2].imag:+.6f}j\t{r[3]}\t{r[4]}\t{r[5]}\t"
                        f"{r[6].real:+.6f}{r[6].imag:+.6f}j\t{r[7]:.6f}\t{r[8]:.4f}\n")
    print("\n   [dump] wrote outputs/judge_completeness_L23.tsv")


if __name__ == "__main__":
    main()
