"""Independent cross-check of project.commuting (seven binary weights).

Nothing is imported from project.commuting.  What differs, and why it guards against
each possible failure of the first computation:

  * Equations: every entry of the explicit matrix product T_N(Q0) T_N(W) - T_N(W) T_N(Q0),
    built from a row enumeration written anew (no reduction modulo rotation).
  * Solving: no lex bases and no elimination.  On each support the slice fixes a
    different set of weights to 1 (the last independent ones instead of the first), the
    Groebner basis uses grevlex and the F5B algorithm, and the conclusion is read off by
    ideal membership: 1 in the ideal (no solution), or the coordinates of the expected
    point in the ideal (a single point, with its multiplicity tested separately).
  * Tangent space: the Jacobian of the explicit equations at Q0, by symbolic
    differentiation, instead of the log-derivative formula.

Run:  uv run python -m project.commuting_crosscheck
"""
from __future__ import annotations

import itertools
import sys

import sympy as sp

z = sp.symbols("z")
w = sp.symbols("w0:7")
sat = sp.symbols("u_sat")
# faces (a b ; c d): a b on top, c d at the bottom; the seven admissible ones
ADM = {(0, 0, 0, 0): 0, (0, 0, 0, 1): 1, (0, 0, 1, 0): 2, (0, 1, 0, 0): 3,
       (0, 1, 1, 0): 4, (1, 0, 0, 0): 5, (1, 0, 0, 1): 6}
CHAR = {k: (1, f[0] - f[1], f[2] - f[3]) for f, k in ADM.items()}


def rows(N):
    out = []
    for bits in range(2 ** N):
        r = tuple((bits >> k) & 1 for k in range(N))
        if not any(r[k] and r[(k + 1) % N] for k in range(N)):
            out.append(r)
    return out


def tm(weights, N):
    R = rows(N)

    def entry(x, y):  # x bottom, y top
        val = sp.Integer(1)
        for k in range(N):
            f = (y[k], y[(k + 1) % N], x[k], x[(k + 1) % N])
            if f not in ADM:
                return sp.Integer(0)
            val *= weights[ADM[f]]
        return val

    return sp.Matrix([[entry(x, y) for y in R] for x in R])


def equations(N, target):
    K, T = tm(target, N), tm(w, N)
    E = (K * T - T * K).expand()
    return sorted({e for e in E if e != 0}, key=sp.default_sort_key)


def slice_from_end(S):
    fixed = []
    for j in reversed(S):
        if sp.Matrix([CHAR[g] for g in fixed + [j]]).rank() > len(fixed):
            fixed.append(j)
    return sorted(fixed)


def closed_walk(S, N):
    st = {(0, 0): 0, (1, 0): 1, (0, 1): 2}
    M = sp.zeros(3, 3)
    for f, k in ADM.items():
        if k in S:
            M[st[(f[0], f[2])], st[(f[1], f[3])]] = 1
    return (M ** N).trace() != 0


def ideal_on(eqs, S, symbolic):
    fixed = slice_from_end(S)
    free = [j for j in S if j not in fixed]
    sub = {w[j]: 0 for j in range(7) if j not in S}
    sub.update({w[g]: 1 for g in fixed})
    E = [e for e in (sp.expand(e.subs(sub)) for e in eqs) if e != 0]
    gens = [w[j] for j in free] + ([z] if symbolic else [])
    return fixed, free, E, gens


def no_solution(eqs, S, symbolic):
    fixed, free, E, gens = ideal_on(eqs, S, symbolic)
    if not E:
        return False
    if not gens:
        return True
    gb = sp.groebner(E + [sat * sp.Mul(*gens) - 1], sat, *gens, order="grevlex", method="f5b")
    return list(gb.exprs) == [1]


def point_only(eqs, S, symbolic, point):
    """The saturated ideal contains w_j - point[j] for every free weight (the slice point is unique)."""
    fixed, free, E, gens = ideal_on(eqs, S, symbolic)
    gb = sp.groebner(E + [sat * sp.Mul(*gens) - 1], sat, *gens, order="grevlex", method="f5b")
    return all(gb.reduce(sp.expand(w[j] - point[j]))[1] == 0 for j in free), gb, free


def proportional(T, M):
    nz = [(i, j) for i in range(M.rows) for j in range(M.cols) if M[i, j] != 0]
    if not nz or any(T[i, j] == 0 for i, j in nz):
        return False
    c = T[nz[0]] / M[nz[0]]
    return (T - c * M).applyfunc(sp.simplify) == sp.zeros(M.rows, M.cols)


def quotient_dimension(gb, gens):
    """Number of standard monomials of a zero-dimensional ideal (grevlex basis)."""
    lead = [sp.Poly(g, *gens).monoms(order="grevlex")[0] for g in gb.exprs]
    bound = []
    for i in range(len(gens)):
        pure = [m[i] for m in lead if all(m[k] == 0 for k in range(len(gens)) if k != i)]
        if not pure:
            return None  # not zero-dimensional
        bound.append(min(pure))
    return sum(1 for e in itertools.product(*[range(b) for b in bound])
               if not any(all(e[k] >= m[k] for k in range(len(gens))) for m in lead))


def check_size(N, zval):
    """Supports with solutions, and for each solution point the class of T_N(W).

    Returns (supports, classes, reduced, number of equations): classes maps a support to the
    set of names among Q0, S, S^-1 to which T_N(W) is proportional at the points of the slice
    (None if a point matches none of them), and reduced tells whether the saturated ideal
    contains the linear equations of its points (a radical, zero-dimensional ideal).
    """
    target = [1, 1, 1, 1, 1, zval, zval]
    eqs = equations(N, target)
    ref = {"Q0": tm(target, N), "S": tm([1, 1, 0, 0, 1, 1, 1], N), "S^-1": tm([1, 0, 1, 1, 1, 0, 1], N)}
    supports, classes, reduced = [], {}, {}
    for k in range(7, 0, -1):
        for S in itertools.combinations(range(7), k):
            if not closed_walk(S, N) or no_solution(eqs, S, False):
                continue
            supports.append(S)
            fixed, free, E, gens = ideal_on(eqs, S, False)
            gb = sp.groebner(E + [sat * sp.Mul(*gens) - 1], sat, *gens, order="grevlex", method="f5b")
            pts = sp.solve(list(gb.exprs), [sat] + gens, dict=True)
            names = set()
            for p in pts:
                if set(p) != set(gens) | {sat}:
                    names.add(None)   # positive-dimensional: not expected here
                    continue
                wt = [p[w[j]] if j in free else (1 if j in fixed else 0) for j in range(7)]
                T = tm(wt, N)
                hit = [n for n, M in ref.items() if proportional(T, M)]
                names.add(hit[0] if hit else None)
            classes[S] = names
            # zero-dimensional ideal: it is radical iff the quotient has dimension = number of points
            reduced[S] = quotient_dimension(gb, [sat] + gens) == len(pts)
    return supports, classes, reduced, len(eqs)


def tangent_rank(N, zval):
    target = [1, 1, 1, 1, 1, zval, zval]
    eqs = equations(N, target)
    at = {w[j]: target[j] for j in range(7)}
    J = sp.Matrix([[sp.diff(e, w[j]).subs(at) for j in range(7)] for e in eqs])
    return J.rank()


def main() -> int:
    fails = []

    def check(label, cond):
        print(f"[{'PASS' if cond else 'FAIL'}] {label}", flush=True)
        if not cond:
            fails.append(label)

    evident = {(0, 1, 2, 3, 4, 5, 6): {"Q0"}, (0, 1, 4, 5, 6): {"S"}, (0, 2, 3, 4, 6): {"S^-1"}}
    for zval in (sp.Integer(3), sp.Rational(-2, 7), sp.Rational(5, 2)):
        for N in (4, 5, 6):
            supports, classes, reduced, neq = check_size(N, zval)
            want_reduced = {S: not (N == 4 and S == (0, 1, 2, 3, 4, 5, 6)) for S in supports}
            check(f"z={zval}, N={N} ({neq} explicit equations): solutions only on the supports of Q0, "
                  f"delta_{{b,c}}, delta_{{a,d}}; at every solution T_N(W) is a multiple of T_N(Q0), S or S^-1; "
                  f"{'Q0 is a double point' if N == 4 else 'all points reduced'}",
                  classes == evident and reduced == want_reduced)
    for zval in (sp.Integer(3), sp.Rational(-2, 7)):
        ranks = {N: tangent_rank(N, zval) for N in (4, 5, 6)}
        check(f"z={zval}: Jacobian ranks at Q0 {ranks}, i.e. tangent dimensions 4, 3, 3", ranks == {4: 3, 5: 4, 6: 4})
    print("All checks passed." if not fails else f"{len(fails)} check(s) FAILED.")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
