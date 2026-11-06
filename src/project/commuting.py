"""Row transfer matrices that commute with the hard-square one, for small sizes.

Question (Section 5.2 of the J. Phys. A version of the paper): which binary faces W,
with weights on the seven admissible faces, have a row transfer matrix commuting
with that of a pure hard-square face Q,

    [T_N(Q), T_N(W)] = 0 ?

No local relation (Yang--Baxter or otherwise) is assumed.  The protocol is
config/simulations/commuting-small-n.yaml and the results are discussed in
docs/simulations/commuting-small-n.md.  This module is not part of the arXiv
ancillary file; verify_certificates.py is unchanged.

Conventions are those of verify_certificates.py.  A face (a b ; c d) has top row
a b and bottom row c d; the admissible faces are labelled

    w0 = 0000, w1 = 0001, w2 = 0010, w3 = 0100, w4 = 0110, w5 = 1000, w6 = 1001,

and the pure hard-square face is Q = (1, D, C, B, BC, A, AD), with z = ABCD.  On
periodic hard-square rows of length N,

    T_N(W)(x, y) = prod_i W[y_i y_{i+1} ; x_i x_{i+1}]    (bottom row x, top row y).

Facts used, each checked by `check_reduction`:
  * T_N(Q) = D_CD C_N D_AB, where C_N is the compatibility matrix of hard-square rows
    and D_t = diag(t^|x|).  Conjugating by D_CD turns the problem for Q into the
    problem for Q0 = (1, 1, 1, 1, 1, z, z), with W replaced by W (CD)^(a - c).  So
    only z matters, and all computations below use Q0.
  * Every T_N is unchanged by W -> W t^(a - b) s^(c - d) (two "corner gauges"),
    and W -> lambda W multiplies it by lambda^N.  These three directions are the
    trivial symmetries of the problem.
  * delta_{b,c} = (1,1,0,0,1,1,1) and delta_{a,d} = (1,0,1,1,1,0,1) give the shift
    S and its inverse, which commute with T_N(Q).

Run:  uv run python -m project.commuting [--long]
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from fractions import Fraction

import sympy as sp

FACES = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (0, 1, 1, 0), (1, 0, 0, 0), (1, 0, 0, 1)]
LABEL = {f: k for k, f in enumerate(FACES)}
# characters of (scale, t, s) on the weights: W -> lambda t^(a-b) s^(c-d) W
CHI = [(1, a - b, c - d) for a, b, c, d in FACES]
W = sp.symbols("w0:7")
z = sp.symbols("z")
_t = sp.symbols("t_sat")

Q0 = [1, 1, 1, 1, 1, z, z]
SHIFT = [1, 1, 0, 0, 1, 1, 1]          # delta_{b,c}
SHIFT_INV = [1, 0, 1, 1, 1, 0, 1]      # delta_{a,d}
HARD_HEXAGONS = [1, 1, 1, 1, 1, z, 0]  # w6 = 0 excludes the a-d diagonal

FULL = (0, 1, 2, 3, 4, 5, 6)
# The components found for N = 5, 6, 7 (symbolic z, z != 0): support, fixed weights, Groebner basis.
EXPECTED = {
    FULL: ([0, 1, 3], [W[2] - 1, W[4] - 1, W[5] - z, W[6] - z]),        # W = Q0
    (0, 1, 4, 5, 6): ([0, 1, 4], [W[5] - 1, W[6] - 1]),                  # W = delta_{b,c}
    (0, 2, 3, 4, 6): ([0, 2, 3], [W[4] - 1, W[6] - 1]),                  # W = delta_{a,d}
}

FAILURES: list[str] = []


def check(label: str, condition: bool) -> None:
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if not condition:
        FAILURES.append(label)


# ---------------------------------------------------------------------------
# Rows, strips and transfer matrices
# ---------------------------------------------------------------------------
def hs_rows(N: int) -> list[tuple[int, ...]]:
    """Periodic hard-square rows of length N (no two adjacent 1s, cyclically)."""
    return [r for r in itertools.product((0, 1), repeat=N)
            if all(not (r[k] and r[(k + 1) % N]) for k in range(N))]


def strip_exponents(bottom, top) -> tuple[int, ...] | None:
    """Face-type counts of the strip with the given rows, or None if a face is not admissible."""
    N = len(bottom)
    c = [0] * 7
    for k in range(N):
        j = LABEL.get((top[k], top[(k + 1) % N], bottom[k], bottom[(k + 1) % N]))
        if j is None:
            return None
        c[j] += 1
    return tuple(c)


def strips(N: int) -> tuple[list, dict]:
    R = hs_rows(N)
    st = {}
    for a in R:
        for b in R:
            e = strip_exponents(a, b)
            if e is not None:
                st[(a, b)] = e
    return R, st


def monomial(weights, e) -> sp.Expr:
    return sp.Mul(*[sp.sympify(weights[j]) ** e[j] for j in range(7) if e[j]])


def transfer_matrix(weights, N: int) -> sp.Matrix:
    R, st = strips(N)
    return sp.Matrix(len(R), len(R), lambda i, j: monomial(weights, st[(R[i], R[j])])
                     if (R[i], R[j]) in st else 0)


def column_graph_has_closed_walk(support, N: int) -> bool:
    """T_N(W) is not identically zero on a support iff the column graph has a closed walk of length N.

    Column states are (top, bottom) = (0,0), (1,0), (0,1); the face (a b; c d) goes from (a, c) to (b, d).
    """
    state = {(0, 0): 0, (1, 0): 1, (0, 1): 2}
    M = sp.zeros(3, 3)
    for j in support:
        a, b, c, d = FACES[j]
        M[state[(a, c)], state[(b, d)]] = 1
    return (M ** N).trace() != 0


# ---------------------------------------------------------------------------
# Reduction checks
# ---------------------------------------------------------------------------
def check_reduction(N: int) -> bool:
    A, B, C, D = sp.symbols("A B C D")
    R = hs_rows(N)
    n = len(R)
    Q = [1, D, C, B, B * C, A, A * D]
    Cm = sp.Matrix(n, n, lambda i, j: int(all(not (R[i][k] and R[j][k]) for k in range(N))))
    DAB = sp.diag(*[(A * B) ** sum(r) for r in R])
    DCD = sp.diag(*[(C * D) ** sum(r) for r in R])
    TQ = transfer_matrix(Q, N)
    ok = (TQ - DCD * Cm * DAB).expand() == sp.zeros(n, n)
    # conjugation by D_CD: T(Q) -> T(Q0) with z = ABCD, and T(W) -> T(W (CD)^(a-c))
    ok &= (DCD.inv() * TQ * DCD - transfer_matrix([sp.sympify(q).subs(z, A * B * C * D) for q in Q0], N)
           ).expand() == sp.zeros(n, n)
    Wt = [W[j] * (C * D) ** (FACES[j][0] - FACES[j][2]) for j in range(7)]
    ok &= (DCD.inv() * transfer_matrix(W, N) * DCD - transfer_matrix(Wt, N)).expand() == sp.zeros(n, n)
    # corner gauges leave T_N unchanged
    t, s = sp.symbols("t s")
    Wg = [W[j] * t ** (FACES[j][0] - FACES[j][1]) * s ** (FACES[j][2] - FACES[j][3]) for j in range(7)]
    ok &= (transfer_matrix(Wg, N) - transfer_matrix(W, N)).expand() == sp.zeros(n, n)
    # delta_{b,c} and delta_{a,d} give mutually inverse permutations commuting with T_N(Q)
    S, Si = transfer_matrix(SHIFT, N), transfer_matrix(SHIFT_INV, N)
    ok &= all(sum(S.row(i)) == 1 for i in range(n)) and S * Si == sp.eye(n)
    ok &= (TQ * S - S * TQ).expand() == sp.zeros(n, n)
    return bool(ok)


# ---------------------------------------------------------------------------
# Commutator equations
# ---------------------------------------------------------------------------
def commutator_equations(N: int, target) -> list[sp.Expr]:
    """Distinct nonzero entries of [T_N(target), T_N(W)], one per pair of rows up to rotation."""
    R, st = strips(N)
    out, seen, done = [], set(), set()
    for x in R:
        for y in R:
            key = min((x[k:] + x[:k], y[k:] + y[:k]) for k in range(N))
            if key in done:
                continue
            done.add(key)
            poly = sp.Integer(0)
            for m in R:
                if (x, m) in st and (m, y) in st:
                    poly += (monomial(target, st[(x, m)]) * monomial(W, st[(m, y)])
                             - monomial(W, st[(x, m)]) * monomial(target, st[(m, y)]))
            poly = sp.expand(poly)
            if poly != 0 and poly not in seen:
                seen.add(poly)
                out.append(poly)
    return out


def gauge_slice(support) -> list[int]:
    """A maximal set of weights in the support with independent characters (fixed to 1)."""
    fixed: list[int] = []
    for j in support:
        if sp.Matrix([CHI[g] for g in fixed + [j]]).rank() > len(fixed):
            fixed.append(j)
    return fixed


def krull_dimension(gb, gens) -> int:
    lead = [sp.Poly(g, *gens).monoms(order=gb.order)[0] for g in gb.exprs]
    for k in range(len(gens), 0, -1):
        for U in itertools.combinations(range(len(gens)), k):
            if all(any(m[i] for i in range(len(gens)) if i not in U) for m in lead):
                return k
    return 0


def solve_support(eqs, support, symbolic_z: bool):
    """Solutions with exactly this support (nonzero weights), modulo scale and corner gauges.

    Returns (fixed, free, basis, dimension) with basis a lex Groebner basis of the saturated
    ideal in the free weights (and z), or None when there is no solution.
    """
    fixed = gauge_slice(support)
    free = [j for j in support if j not in fixed]
    sub = {W[j]: 0 for j in range(7) if j not in support}
    sub.update({W[g]: 1 for g in fixed})
    E = [e for e in (sp.expand(e.subs(sub)) for e in eqs) if e != 0]
    unknowns = [W[j] for j in free] + ([z] if symbolic_z else [])
    if not E:
        return fixed, free, [], len(unknowns)
    if not unknowns:
        return None
    nonzero = sp.Mul(*unknowns)
    gb = sp.groebner(E + [_t * nonzero - 1], _t, *unknowns, order="lex")
    if list(gb.exprs) == [1]:
        return None
    elim = [g for g in gb.exprs if not g.has(_t)]
    if not elim:
        return fixed, free, [], len(unknowns)
    gb2 = sp.groebner(elim, *unknowns, order="lex")
    return fixed, free, list(gb2.exprs), krull_dimension(gb2, unknowns)


def components(N: int, target, symbolic_z: bool) -> dict:
    """{support: (fixed, free, basis, dimension)} over all supports on which T_N(W) is not zero."""
    eqs = commutator_equations(N, target)
    out = {}
    for k in range(7, 0, -1):
        for S in itertools.combinations(range(7), k):
            if column_graph_has_closed_walk(S, N):
                res = solve_support(eqs, S, symbolic_z)
                if res is not None:
                    out[S] = res
    return out


# ---------------------------------------------------------------------------
# Tangent space at Q0
# ---------------------------------------------------------------------------
def tangent_rows(N: int, zval):
    """Coefficient vectors (in h) of [T_N(Q0), d/de T_N(Q0 (1 + e h))|_{e=0}], divided by z^|y|.

    At Q0 all seven weights are nonzero and the derivative is T_N(Q0) o sum_j h_j n_j, where
    n_j(x, y) counts faces of type j in the strip.  The entry (x, y) of the commutator, divided
    by z^|y|, is  sum_{m compatible with x and y} z^|m| (F(m, y) - F(x, m)),  F = sum_j h_j n_j.
    """
    R, st = strips(N)
    nb = {a: [b for b in R if (a, b) in st] for a in R}
    seen = set()
    for x in R:
        for y in R:
            key = min((x[k:] + x[:k], y[k:] + y[:k]) for k in range(N))
            if key in seen:
                continue
            seen.add(key)
            ys = set(nb[y])
            v = [0] * 7
            for m in nb[x]:
                if m in ys:
                    a, b = st[(m, y)], st[(x, m)]
                    zm = zval ** sum(m)
                    for j in range(7):
                        v[j] += zm * (a[j] - b[j])
            yield v


def tangent_dimension(N: int, zval: Fraction) -> int:
    """Exact dimension of the tangent space of {W : [T_N(Q0), T_N(W)] = 0} at Q0 (z = zval)."""
    basis: list[tuple[int, list[Fraction]]] = []
    for v in tangent_rows(N, zval):
        v = [Fraction(x) for x in v]
        for p, b in basis:
            if v[p]:
                f = v[p] / b[p]
                v = [x - f * y for x, y in zip(v, b)]
        piv = next((k for k in range(7) if v[k]), None)
        if piv is not None:
            basis.append((piv, v))
    return 7 - len(basis)


def tangent_invariant_factors(N: int) -> list[sp.Expr]:
    """Invariant factors over Q[z] of the tangent equations in the nontrivial directions.

    The trivial directions (scale, two corner gauges) are transverse to {h0 = h2 = h3 = 0},
    so the nontrivial tangent space is the kernel of the columns h1, h4, h5, h6.  Over the PID
    Q[z], the product of the invariant factors is the gcd of the 4x4 minors: the tangent space
    is trivial exactly at the z where it does not vanish.
    """
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.matrices.normalforms import invariant_factors

    cols = (1, 4, 5, 6)
    K = sp.QQ[z]
    rows_, seen = [], set()
    for v in tangent_rows(N, z):
        r = tuple(sp.expand(v[c]) for c in cols)
        if any(r) and r not in seen:
            seen.add(r)
            rows_.append([K.from_sympy(x) for x in r])
    inv = invariant_factors(DomainMatrix(rows_, (len(rows_), 4), K))
    return [sp.factor(K.to_sympy(f)) for f in inv]


# ---------------------------------------------------------------------------
# Positive control: hard hexagons
# ---------------------------------------------------------------------------
def baxter_F(weights) -> sp.Expr:
    """F(z_W, X, Y) = z XY(XY - X - Y) - (X - 1)(Y - 1) for a face with w0 w1 w2 w3 w5 != 0."""
    w0, w1, w2, w3, w4, w5, w6 = weights
    X, Y = w0 * w6 / (w1 * w5), w0 * w4 / (w2 * w3)
    zw = w1 * w2 * w3 * w5 / w0 ** 4
    return sp.together(zw * X * Y * (X * Y - X - Y) - (X - 1) * (Y - 1))


def on_baxter_surface(fixed, free, basis) -> bool:
    """The numerator of F reduces to 0 modulo the basis (component with w0..w6 in the slice)."""
    weights = [sp.Integer(1) if j in fixed else W[j] for j in range(7)]
    num = sp.numer(baxter_F(weights))
    gens = [W[j] for j in free]
    return sp.groebner(basis, *gens, order="lex").reduce(sp.expand(num))[1] == 0


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def describe(comp: dict) -> None:
    for S, (fixed, free, basis, dim) in comp.items():
        print(f"    support {S}: fixed {fixed}, free {free}, dimension {dim}, basis {basis}")


def matches_expected(comp: dict) -> bool:
    if set(comp) != set(EXPECTED):
        return False
    return all(comp[S][0] == EXPECTED[S][0]
               and [sp.expand(g) for g in comp[S][2]] == [sp.expand(g) for g in EXPECTED[S][1]]
               for S in EXPECTED)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--long", action="store_true", help="also N = 7 globally and the tangent up to N = 12")
    args = ap.parse_args(argv)
    t0 = time.time()

    print("== Reduction and evident solutions ==")
    for N in range(2, 7):
        check(f"N={N}: T_N(Q) = D_CD C_N D_AB, conjugation to Q0, corner gauges, shift commutes",
              check_reduction(N))

    print("\n== Small sizes: all solutions, symbolic z (z != 0) ==")
    for N in (2, 3):
        comp = components(N, Q0, symbolic_z=True)
        full = comp.get(FULL)
        print(f"  N={N}: {len(comp)} supports with solutions; full support: dimension "
              f"{full[3] - 1 if full else None} for fixed z (modulo scale and corner gauges)")
        check(f"N={N}: positive-dimensional solutions on the full support besides Q0",
              full is not None and full[3] - 1 >= 1)
    comp4 = components(4, Q0, symbolic_z=True)
    describe(comp4)
    check("N=4: solutions only on the supports of Q0, delta_{b,c} and delta_{a,d}", set(comp4) == set(EXPECTED))
    q4 = sp.groebner(comp4[FULL][2], W[2], W[4], W[5], W[6], z, order="lex")

    def in_q4(g):
        return q4.reduce(sp.expand(g))[1] == 0

    check("N=4: on the full support the only point is Q0, and it is not reduced "
          "((w4 - 1)^2 and (w6 - z)^2 lie in the ideal, w4 - 1 does not)",
          all(in_q4(g) for g in [W[2] - 1, W[5] - z, (W[4] - 1) ** 2, (W[6] - z) ** 2]) and not in_q4(W[4] - 1))
    sizes = (5, 6, 7) if args.long else (5, 6)
    for N in sizes:
        comp = components(N, Q0, symbolic_z=True)
        check(f"N={N}: exactly W ~ Q0, delta_{{b,c}}, delta_{{a,d}} for every z != 0, each reduced",
              matches_expected(comp))

    print("\n== Tangent space at Q0 (trivial part: dimension 3) ==")
    for N in (5, 6):
        inv = tangent_invariant_factors(N)
        check(f"N={N}: invariant factors over Q[z] are {inv}; the tangent space is trivial for every z != 0",
              inv == [1, 1, z ** 2, z ** 2])
    Nmax = 12 if args.long else 9
    for zval in (Fraction(2), Fraction(-1), Fraction(1, 3), Fraction(-1, 4), Fraction(7, 5)):
        dims = {N: tangent_dimension(N, zval) for N in range(2, Nmax + 1)}
        check(f"z={zval}: tangent dimensions {dims}; 3 (trivial) for 5 <= N <= {Nmax}",
              all(d == 3 for N, d in dims.items() if N >= 5) and dims[4] == 4)

    print("\n== Positive control: hard hexagons (z = 2) ==")
    hh = [x.subs(z, 2) if isinstance(x, sp.Basic) else x for x in HARD_HEXAGONS]
    bases = []
    for N in ((4, 5, 6) if args.long else (4, 5)):
        comp = components(N, hh, symbolic_z=False)
        fixed, free, basis, dim = comp[FULL]
        bases.append([sp.expand(g) for g in basis])
        check(f"N={N}: the full support carries a curve of commuting partners (dimension {dim}), "
              f"on Baxter's surface F = 0", dim == 1 and on_baxter_surface(fixed, free, basis))
    check("the curve is the same for every N computed", all(b == bases[0] for b in bases))

    print(f"\n{'All checks passed.' if not FAILURES else f'{len(FAILURES)} check(s) FAILED.'} "
          f"[{time.time() - t0:.0f} s]")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
