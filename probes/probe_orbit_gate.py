"""
PROBE ORBIT-GATE -- does the qx+1 transfer operator block-diagonalize by <2>-character-orbits? (2026-07-28)

Wilson's YOLO conjecture (<2>-graded fine structure): the 2-adic renewal r -> (qr+1) 2^-v on G=(Z/q^k)*
commutes with x2, so its transfer operator block-diagonalizes over the character group Ghat, graded by
the <2>-orbit of each character; a block/orbit of size L contributes a denominator 2^L - 1; summing
rational-in-q coefficients over orbits reproduces the fine-structure invariant (c~_q, den theorem).

GATE (falsifier #3, done empirically -- do NOT assume the grading): build the transfer operator K in the
CHARACTER basis and measure which partition of Ghat, if any, makes K block-diagonal.
  candidate A: chi -> chi^2 orbits  (a -> 2a mod phi)         [Wilson's stated mechanism]
  candidate B: x2-translation eigenclasses (equal chi(2))     [dual of mult-by-element-2, the diagonal one]
Report off-block mass for each vs a random-partition null. Whichever kills the off-block mass IS the grading.
If neither -> the character-squaring mechanism is a corpse (a different structure carries the 2^L-1).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from c_tilde_q17_probe import order_of_two


def generator(N, q, phi):
    def is_gen(g):
        x = 1; seen = 0; period = 0
        for _ in range(phi):
            x = (x * g) % N; period += 1
            if x == 1:
                break
        return period == phi
    for g in range(2, N):
        if g % q != 0 and is_gen(g):
            return g
    raise RuntimeError("no generator (G not cyclic?)")


def build_K_char(q, k, vmax=60):
    N = q ** k
    phi = q ** (k - 1) * (q - 1)
    M = order_of_two(N)
    vm = min(M, vmax)
    inv2 = pow(2, -1, N)
    powinv2 = [pow(inv2, v, N) for v in range(1, vm + 1)]
    w = np.array([0.5 ** v for v in range(1, vm + 1)]); w /= w.sum()
    units = [r for r in range(N) if r % q != 0]
    idx = {r: i for i, r in enumerate(units)}
    n = len(units)
    K = np.zeros((n, n))
    for r in units:
        for vi in range(vm):
            tgt = ((q * r + 1) * powinv2[vi]) % N
            K[idx[r], idx[tgt]] += w[vi]
    # discrete log wrt a generator
    g = generator(N, q, phi)
    dlog = {}; x = 1
    for i in range(phi):
        dlog[x] = i; x = (x * g) % N
    s = dlog[2 % N]                                  # ind(2); x2 = +s in exponent
    # unitary character matrix  F[unit_j, a] = exp(2pi i a dlog(u_j)/phi)/sqrt(phi)
    dl = np.array([dlog[r] for r in units])
    a = np.arange(phi)
    F = np.exp(2j * np.pi * np.outer(dl, a) / phi) / np.sqrt(phi)
    Kc = F.conj().T @ K @ F                          # transfer operator in character basis
    return Kc, phi, s, M


def orbits_squaring(phi):
    """chi -> chi^2 orbits: a -> 2a mod phi."""
    seen = [False] * phi; orbs = []
    for a in range(phi):
        if seen[a]:
            continue
        o = []; b = a
        while not seen[b]:
            seen[b] = True; o.append(b); b = (2 * b) % phi
        orbs.append(o)
    return orbs


def classes_translation(phi, s):
    """x2-translation eigenclasses: group a by value (a*s mod phi) [equal chi(2) eigenvalue]."""
    d = {}
    for a in range(phi):
        d.setdefault((a * s) % phi, []).append(a)
    return list(d.values())


def offblock_mass(Kc, partition):
    """fraction of |Kc|^2 mass lying BETWEEN different blocks."""
    lab = np.empty(Kc.shape[0], dtype=int)
    for bi, blk in enumerate(partition):
        for a in blk:
            lab[a] = bi
    P = np.abs(Kc) ** 2
    tot = P.sum()
    same = 0.0
    for a in range(Kc.shape[0]):
        same += P[a, lab == lab[a]].sum()
    return 1.0 - same / tot


def main():
    print("# PROBE ORBIT-GATE -- does K block-diagonalize by <2>-character-orbits?\n")
    rng = np.random.default_rng(0)
    for q, k in ((7, 2), (13, 2)):
        Kc, phi, s, M = build_K_char(q, k)
        orbsA = orbits_squaring(phi)
        clsB = classes_translation(phi, s)
        offA = offblock_mass(Kc, orbsA)
        offB = offblock_mass(Kc, clsB)
        # random-partition null with the same block-size multiset as A
        sizesA = sorted(len(o) for o in orbsA)
        perm = rng.permutation(phi); rp = []; i = 0
        for sz in sizesA:
            rp.append(list(perm[i:i + sz])); i += sz
        offRnd = offblock_mass(Kc, rp)
        print(f"## q={q}, k={k}: |G|=phi={phi}, ord(2 mod q^k)=M={M}, ind(2)=s={s}")
        print(f"   candidate A  chi->chi^2 orbits ({len(orbsA)} blocks, sizes {sorted(set(len(o) for o in orbsA))}):"
              f"  off-block mass = {offA:.4e}")
        print(f"   candidate B  x2-translation eigenclasses ({len(clsB)} blocks of {len(clsB[0])}):"
              f"  off-block mass = {offB:.4e}")
        print(f"   random-partition null (same sizes as A):                       off-block mass = {offRnd:.4e}")
        winner = "A (chi->chi^2)" if offA < 1e-9 else "B (x2-eigenclass)" if offB < 1e-9 else "NEITHER"
        print(f"   => block-diagonal under: {winner}")
        # DIAGNOSIS: does K commute with x2 (the mechanism's premise)? and is the +1 the culprit?
        c_with, c_without = commute_check(q, k)
        print(f"   commute check ||K P2 - P2 K||_max: WITH +1 = {c_with:.3e}  |  WITHOUT +1 (qr only, on 1-units) = {c_without:.3e}")
        print(f"   => the +1 {'BREAKS' if c_with > 1e-9 else 'preserves'} x2-equivariance"
              f" (mechanism premise {'FALSE' if c_with > 1e-9 else 'holds'}); "
              f"the multiplicative-only map {'DOES' if c_without < 1e-9 else 'does not'} commute.\n")


def commute_check(q, k, vmax=60):
    """||K P2 - P2 K|| for the real map (qr+1) and for the multiplicative-only surrogate r->r*2^-v
       restricted to the 1-units U_1 = {1 mod q} (a genuine <2>-invariant subgroup), to isolate the +1."""
    N = q ** k; M = order_of_two(N); vm = min(M, vmax); inv2 = pow(2, -1, N)
    powinv2 = [pow(inv2, v, N) for v in range(1, vm + 1)]
    w = np.array([0.5 ** v for v in range(1, vm + 1)]); w /= w.sum()
    units = [r for r in range(N) if r % q != 0]; idx = {r: i for i, r in enumerate(units)}
    n = len(units)
    P2 = np.zeros((n, n))
    for r in units:
        P2[idx[(2 * r) % N], idx[r]] = 1.0
    K = np.zeros((n, n))
    for r in units:
        for vi in range(vm):
            K[idx[((q * r + 1) * powinv2[vi]) % N], idx[r]] += w[vi]
    c_with = np.abs(K @ P2 - P2 @ K).max()
    # pure-halving surrogate on ALL units: drop the whole affine (qr+1), map r -> r*2^-v (stays a unit).
    # This is the <2>-equivariant map; if it commutes and the real one doesn't, the affine qr+1 is the culprit.
    Kb = np.zeros((n, n))
    for r in units:
        for vi in range(vm):
            Kb[idx[(r * powinv2[vi]) % N], idx[r]] += w[vi]
    c_without = np.abs(Kb @ P2 - P2 @ Kb).max()
    return c_with, c_without


if __name__ == "__main__":
    main()
