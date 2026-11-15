"""Commutation with the hard-square transfer matrix on the alphabet A_4, near the presentation.

Paper 1 (arXiv:2611.06780), Section 5.2 of the J. Phys. A version, and Remark 2.11.  On A_4 a
face weight is a pair (W^0, W^1) of binary face weights on the seven admissible faces, W^tau
acting on the faces whose top-left height has parity tau.  On a cylinder of even circumference
M, a row of heights is a hard-square row x together with the parity p of its first height (site
k has parity p + k mod 2).  Consecutive rows have opposite parities, and in the strip from a
bottom row (x, 1 - p) to a top row (y, p) the face at column k has type p + k mod 2.  The entry
(bottom, top) of T_M(W) is the product of the M faces of the strip.

Question.  Let Q be a presentation of hard squares on A_4, that is a staggered pure hard-square
face of activity z != 0 (Proposition 2.9(2)).  Is there an analytic family W(u) of A_4 face
weights with W(u0) = Q and [T_M(Q), T_M(W(u))] = 0 for all u, other than the families with
T_M(W(u)) proportional to T_M(Q)?  No local relation is assumed.  This is the transfer-matrix
counterpart, on A_4, of Theorem B.  Partners far from Q are not addressed: with fourteen
weights, the global elimination of project.commuting is out of reach.

Method.
  * Reduction.  Corner factors and a diagonal similarity of the rows, both realized on faces,
    carry every presentation of activity z to Q0 (+) Q0, with Q0 = (1,1,1,1,1,z,z) on both
    types (check_reduction).
  * Trivial symmetries, which multiply T_M(W) by a scalar: the overall scale, W^0 -> mu W^0
    with W^1 -> W^1 / mu, and four corner gauges, two on the top row and two on the bottom
    row.  Their orbit through Q lies in the solution set and has dimension 6.
  * Tangent space.  If the Zariski tangent space of the solution set at Q has dimension 6, the
    solution set near Q is that orbit, so T_M(W(u)) is proportional to T_M(Q) near u0, and by
    analyticity for all u.  In log coordinates W' = Q h the trivial directions are transverse
    to {h^tau_0 = h^tau_2 = h^tau_3 = 0}, so the tangent space is trivial exactly when the
    eight columns h^tau_1, h^tau_4, h^tau_5, h^tau_6 have rank 8.  Exact ranks at fixed z;
    over Q[z], the gcd of their 8 x 8 minors is a power of z, so the rank is 8 for every z != 0.
  * Control.  With hard hexagons (1,1,1,1,1,z,0) on both types the tangent space has one more
    direction, and it is unstaggered: the direction of Baxter's family on the binary alphabet.

Run:  uv run python -m project.commuting_a4 [--long]   (--long adds M = 10)
"""
from __future__ import annotations

import random
import sys
import time
from fractions import Fraction

import sympy as sp

z = sp.symbols("z")
# faces (a b ; c d): a b on top, c d at the bottom; the seven admissible ones, labelled w0..w6
FACES = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (0, 1, 1, 0), (1, 0, 0, 0), (1, 0, 0, 1)]
IDX = {f: i for i, f in enumerate(FACES)}
# coordinate c = 7 tau + j is the weight of W^tau on face j
NONTRIVIAL = [7 * t + j for t in (0, 1) for j in (1, 4, 5, 6)]
TANGENT_ACTIVITIES = [sp.Integer(2), sp.Integer(-1), sp.Rational(1, 3), sp.Rational(-1, 4), sp.Rational(7, 5)]


def q0(zval):
    return [1, 1, 1, 1, 1, zval, zval]


# ---------------------------------------------------------------------------
# Rows, strips and transfer matrices
# ---------------------------------------------------------------------------
def hs_rows(M):
    out = []
    for bits in range(2 ** M):
        r = tuple((bits >> k) & 1 for k in range(M))
        if not any(r[k] and r[(k + 1) % M] for k in range(M)):
            out.append(r)
    return out


def rows(M):
    """Rows of heights of A_4 on a cylinder of even circumference M: (hard-square row, parity)."""
    assert M % 2 == 0
    return [(x, p) for p in (0, 1) for x in hs_rows(M)]


def strips(M):
    """{(bottom index, top index): [(type, face label)] over the M columns} for admissible strips."""
    R = rows(M)
    S = {}
    for i, (x, pb) in enumerate(R):
        for j, (y, pt) in enumerate(R):
            if pt == pb:
                continue
            faces = []
            for k in range(M):
                f = (y[k], y[(k + 1) % M], x[k], x[(k + 1) % M])
                if f not in IDX:
                    break
                faces.append(((pt + k) % 2, IDX[f]))
            else:
                S[(i, j)] = faces
    return R, S


def tm(W0, W1, M):
    """T_M(W) for W = (W^0, W^1), as a sympy matrix."""
    R, S = strips(M)
    T = sp.zeros(len(R), len(R))
    for (i, j), faces in S.items():
        val = sp.Integer(1)
        for t, f in faces:
            val *= (W0, W1)[t][f]
        T[i, j] = val
    return T


def corner_map(W, fa, fb, fc, fd):
    """W[a b; c d] fa^a fb^b fc^c fd^d on the seven faces."""
    return [W[i] * fa ** a * fb ** b * fc ** c * fd ** d for i, (a, b, c, d) in enumerate(FACES)]


def symmetry(W0, W1, lam, mu, al0, al1, ep0, ep1):
    """The trivial symmetries; T_M is multiplied by lam^M."""
    return (corner_map([lam * mu * w for w in W0], al0, 1 / al1, ep0, 1 / ep1),
            corner_map([lam / mu * w for w in W1], al1, 1 / al0, ep1, 1 / ep0))


def similarity(W0, W1, de0, de1):
    """Faces realizing T_M(W) -> D T_M(W) D^-1, with D(x, p) = prod over occupied sites k of de_{p+k}."""
    return corner_map(W0, 1 / de0, 1, de1, 1), corner_map(W1, 1 / de1, 1, de0, 1)


def row_factor(row, de0, de1):
    x, p = row
    out = sp.Integer(1)
    for k, occ in enumerate(x):
        if occ:
            out *= (de0, de1)[(p + k) % 2]
    return out


def check_reduction(M=4) -> bool:
    """Symmetries, the similarity, and the reduction of every presentation to Q0 (+) Q0."""
    a = sp.symbols("a0:7")
    b = sp.symbols("b0:7")
    lam, mu, al0, al1, ep0, ep1, de0, de1 = sp.symbols("lam mu al0 al1 ep0 ep1 de0 de1", nonzero=True)
    T = tm(a, b, M)
    ok = (tm(*symmetry(a, b, lam, mu, al0, al1, ep0, ep1), M) - lam ** M * T).applyfunc(sp.cancel) == sp.zeros(*T.shape)
    R = rows(M)
    D = sp.diag(*[row_factor(r, de0, de1) for r in R])
    ok &= (tm(*similarity(a, b, de0, de1), M) - D * T * D.inv()).applyfunc(sp.cancel) == sp.zeros(*T.shape)
    # a general presentation (s_tau A_tau^a B_tau^b C_tau^c D_tau^d) with A_t B_{1-t} C_{1-t} D_t = z
    s0, s1, B0, C0, D0, B1, C1, D1 = sp.symbols("s0 s1 B0 C0 D0 B1 C1 D1", nonzero=True)
    A0, A1 = z / (B1 * C1 * D0), z / (B0 * C0 * D1)
    one = [1] * 7
    P0 = corner_map([s0 * w for w in one], A0, B0, C0, D0)
    P1 = corner_map([s1 * w for w in one], A1, B1, C1, D1)
    # scale each type (an element of the torus generated by lam and mu), then the corner gauges
    # al_t = B_{1-t}, ep_t = D_{1-t}, then the similarity de_t = 1 / (C_{1-t} D_t)
    G0, G1 = [w / s0 for w in P0], [w / s1 for w in P1]
    G0, G1 = symmetry(G0, G1, sp.Integer(1), sp.Integer(1), B1, B0, D1, D0)
    G0, G1 = similarity(G0, G1, 1 / (C1 * D0), 1 / (C0 * D1))
    ok &= [sp.cancel(w) for w in G0] == q0(z) and [sp.cancel(w) for w in G1] == q0(z)
    return bool(ok)


# ---------------------------------------------------------------------------
# Sparse exact linear algebra
# ---------------------------------------------------------------------------
def by_row(A):
    out = {}
    for (i, j), v in A.items():
        out.setdefault(i, []).append((j, v))
    return out


def commutator(A, B):
    """AB - BA for sparse matrices given as {(i, j): value}."""
    Ar, Br = by_row(A), by_row(B)
    C = {}
    for X, Y, sign in ((Ar, Br, 1), (Br, Ar, -1)):
        for i, row in X.items():
            for j, x in row:
                for k, y in Y.get(j, ()):
                    C[(i, k)] = C.get((i, k), 0) + sign * x * y
    return {key: v for key, v in C.items() if v != 0}


def echelon(columns):
    """Rank and pivot rows of the matrix whose columns are sparse dicts over the same keys."""
    keys = sorted(set().union(*columns))
    pivots = []
    for key in keys:
        v = [Fraction(c.get(key, 0)) for c in columns]
        for p, b in pivots:
            if v[p]:
                f = v[p] / b[p]
                v = [vi - f * bi for vi, bi in zip(v, b)]
        nz = next((i for i, x in enumerate(v) if x), None)
        if nz is not None:
            pivots.append((nz, v))
            if len(pivots) == len(columns):
                break
    return len(pivots), [v for _, v in pivots]


def frac(x):
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


# ---------------------------------------------------------------------------
# Tangent space at a numeric face weight (plain coordinates)
# ---------------------------------------------------------------------------
def jacobian(W0, W1, M):
    """K = T_M(W) and the fourteen derivatives dT_M/dW^tau_j, sparse, at numeric weights."""
    W = ([frac(w) for w in W0], [frac(w) for w in W1])
    _, S = strips(M)
    K, D = {}, [dict() for _ in range(14)]
    for key, faces in S.items():
        vals = [W[t][f] for t, f in faces]
        pre, suf = [Fraction(1)], [Fraction(1)]
        for v in vals:
            pre.append(pre[-1] * v)
        for v in reversed(vals):
            suf.append(suf[-1] * v)
        suf.reverse()
        if pre[-1]:
            K[key] = pre[-1]
        for k, (t, f) in enumerate(faces):
            d = pre[k] * suf[k + 1]
            if d:
                c = 7 * t + f
                D[c][key] = D[c].get(key, 0) + d
    return K, D


def orbit_generators(W0, W1):
    """Derivatives at the identity of the six trivial symmetries, as 14-vectors."""
    W0, W1 = [frac(w) for w in W0], [frac(w) for w in W1]
    A = [f[0] for f in FACES]
    B = [f[1] for f in FACES]
    C = [f[2] for f in FACES]
    Dd = [f[3] for f in FACES]
    g = lambda u, v: [u[j] * W0[j] for j in range(7)] + [v[j] * W1[j] for j in range(7)]
    one, neg = [1] * 7, [-1] * 7
    return [g(one, one), g(one, neg), g(A, [-x for x in B]), g([-x for x in B], A),
            g(C, [-x for x in Dd]), g([-x for x in Dd], C)]


def vectors_rank(vectors):
    """Rank of a list of 14-vectors (each one a column)."""
    return echelon([{c: x for c, x in enumerate(v) if x} for v in vectors])[0]


def tangent(W0, W1, M):
    """(dimension of the tangent space of the solution set at W, dimension of the trivial orbit, kernel)."""
    K, D = jacobian(W0, W1, M)
    cols = [commutator(K, Dc) for Dc in D]
    rank, piv = echelon(cols)
    kernel = sp.Matrix(piv).nullspace() if piv else [sp.eye(14)[:, c] for c in range(14)]
    return 14 - rank, vectors_rank(orbit_generators(W0, W1)), kernel


def unstaggered_extra(W0, W1, kernel) -> bool:
    """Every tangent vector is a trivial direction plus an unstaggered one (W'^0 = W'^1)."""
    base = orbit_generators(W0, W1) + [[int(c == j or c == 7 + j) for c in range(14)] for j in range(7)]
    with_kernel = base + [[frac(x) for x in k] for k in kernel]
    return vectors_rank(with_kernel) == vectors_rank(base)


# ---------------------------------------------------------------------------
# Tangent space at Q0 (+) Q0 in log coordinates, with z symbolic
# ---------------------------------------------------------------------------
def log_columns(M):
    """Columns [K, N_c o K] of the tangent equations at Q0 (+) Q0 in log coordinates.

    K(i, j) = z^|top row j| on admissible strips, and N_c(i, j) counts the faces of coordinate c
    in the strip.  Entries are integer polynomials in z, stored as {exponent: coefficient}.
    """
    R, S = strips(M)
    top = [sum(x) for x, _ in R]
    n = [dict() for _ in range(14)]
    for key, faces in S.items():
        for t, f in faces:
            n[7 * t + f][key] = n[7 * t + f].get(key, 0) + 1
    adj = by_row({key: 1 for key in S})
    cols = []
    for c in range(14):
        C = {}
        for i, row in adj.items():
            for j, _ in row:
                for k, _ in adj.get(j, ()):
                    coef = n[c].get((j, k), 0) - n[c].get((i, j), 0)
                    if coef:
                        poly = C.setdefault((i, k), {})
                        e = top[j] + top[k]
                        poly[e] = poly.get(e, 0) + coef
        cols.append({key: {e: v for e, v in p.items() if v} for key, p in C.items()})
    return [{key: p for key, p in col.items() if p} for col in cols]


def evaluate(poly, zval):
    zf = frac(zval)
    return sum(Fraction(v) * zf ** e for e, v in poly.items())


def log_tangent_dimension(M, zval, cols=None) -> int:
    cols = cols or log_columns(M)
    num = [{key: evaluate(p, zval) for key, p in col.items()} for col in cols]
    return 14 - echelon(num)[0]


def nontrivial_minors_gcd(M, cols=None, seed=0) -> sp.Expr:
    """A multiple of the gcd of the 8 x 8 minors of the nontrivial columns, over Q[z].

    Rows are added until the product of the invariant factors is a power of z; since that
    product divides the gcd for the chosen rows only, the gcd for all rows divides it.
    """
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.matrices.normalforms import invariant_factors

    cols = cols or log_columns(M)
    sub = [cols[c] for c in NONTRIVIAL]
    keys = sorted(set().union(*sub))
    rng = random.Random(seed)
    rng.shuffle(keys)
    K = sp.QQ[z]
    to_poly = lambda p: K.from_sympy(sum(v * z ** e for e, v in p.items())) if p else K.zero
    chosen, seen = [], set()
    for start in range(0, len(keys), 40):
        for key in keys[start:start + 40]:
            r = tuple(tuple(sorted(col.get(key, {}).items())) for col in sub)
            if any(r) and r not in seen:
                seen.add(r)
                chosen.append([to_poly(col.get(key, {})) for col in sub])
        if len(chosen) < 8:
            continue
        inv = invariant_factors(DomainMatrix(chosen, (len(chosen), 8), K))
        if len(inv) == 8:
            prod = sp.factor(sp.prod([K.to_sympy(f) for f in inv]))
            if prod != 0 and sp.Poly(prod, z).is_monomial:
                return prod
    return sp.nan


def staggered_presentation(zval):
    """A presentation with type-dependent corner factors and s_1 / s_0 = 3/2 (A_t B_{1-t} C_{1-t} D_t = z)."""
    B0, C0, D0 = sp.Rational(1, 2), sp.Integer(2), sp.Integer(5)
    A1, B1, C1 = sp.Rational(2, 3), sp.Rational(7, 3), sp.Rational(1, 4)
    A0, D1 = zval / (B1 * C1 * D0), zval / (A1 * B0 * C0)
    return corner_map([1] * 7, A0, B0, C0, D0), corner_map([sp.Rational(3, 2)] * 7, A1, B1, C1, D1)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    long = "--long" in argv
    sizes = [6, 8, 10] if long else [6, 8]
    fails = []
    t0 = time.time()

    def check(label, cond):
        print(f"[{'PASS' if cond else 'FAIL'}] {label}", flush=True)
        if not cond:
            fails.append(label)

    print("== Reduction and symmetries ==")
    check("M=4: the six trivial symmetries multiply T_M by lam^M, the face similarity conjugates T_M by "
          "a diagonal matrix, and every presentation of activity z is carried to Q0 (+) Q0", check_reduction(4))

    print("\n== Tangent space at Q0 (+) Q0 ==")
    for M in [4] + sizes:
        cols = log_columns(M)
        dims = {str(zv): log_tangent_dimension(M, zv, cols) for zv in TANGENT_ACTIVITIES}
        want = 7 if M == 4 else 6
        check(f"M={M} ({len(rows(M))} rows): tangent dimension {dims} "
              f"({'one direction beyond the six trivial ones' if M == 4 else 'only the six trivial directions'})",
              set(dims.values()) == {want})
    tan, orb, _ = tangent(q0(2), q0(2), 6)
    check(f"M=6, z=2: plain coordinates agree (tangent {tan}), and the trivial orbit has dimension {orb}",
          tan == 6 and orb == 6)

    for M in (6,):
        g = nontrivial_minors_gcd(M)
        check(f"M={M}: over Q[z], the gcd of the 8 x 8 minors of the nontrivial columns divides {g}, "
              f"so the tangent space is trivial for every z != 0", g is not sp.nan)

    print("\n== A staggered presentation (no reduction used) ==")
    P0, P1 = staggered_presentation(sp.Integer(2))
    tan, orb, _ = tangent(P0, P1, 6)
    check(f"M=6, z=2, staggered corner factors and s_1/s_0 = 3/2: tangent {tan}, trivial orbit {orb}",
          tan == 6 and orb == 6)

    print("\n== Positive control: hard hexagons on both types ==")
    H = [1, 1, 1, 1, 1, 2, 0]
    for M in (6, 8):
        tan, orb, ker = tangent(H, H, M)
        check(f"M={M}, z=2: tangent {tan}, trivial orbit {orb}: one more direction, and it is unstaggered",
              tan == 7 and orb == 6 and unstaggered_extra(H, H, ker))

    print(("\nAll checks passed." if not fails else f"\n{len(fails)} check(s) FAILED.") + f" [{time.time() - t0:.0f} s]")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
