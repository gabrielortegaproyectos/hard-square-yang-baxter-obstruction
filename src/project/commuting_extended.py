"""Commuting partners of hard squares beyond the seven admissible weights.

Extension of project.commuting to the other binary faces (issue #1).

Nine weights.  The faces (10;10) and (01;01) keep both rows hard-square rows; they let a
site be occupied in two consecutive rows, and each such vertical pair uses one of each,
so they enter T_N only through the product w7 w8.  The transfer matrix still acts on
hard-square rows.  Besides multiples of K = T_N(Q0), S and S^-1, the evident partners
are now the identity I (W = 1 on 0000, 1010, 0101) and K S, K S^-1 (hard squares in which
the vertical exclusion is replaced by a diagonal one).

Sixteen faces.  With weights on faces with two horizontally adjacent occupied sites,
T_N(W) acts on all binary rows.  The hard-square face extended by zero, Q-hat, has
T_N(Q-hat) = T_N(Q) (+) 0 (hard-square rows first).  Since det C_N = +-1, T_N(Q) is
invertible for z != 0, so [T_N(Q-hat), T_N(W)] = 0 says exactly: the hard-square block of
T_N(W), which depends only on the nine weights above, commutes with T_N(Q), and the two
off-diagonal blocks vanish.  The block on the other rows is unconstrained.  For example
W_u = Q-hat + u [1111] has T_N(W_u) = T_N(Q) (+) u^N E, with E the projection on the
all-occupied row: a one-parameter family through hard squares whose transfer matrices
commute for every N, which by Proposition 4.3 of the paper is not a Yang--Baxter family.

Run:  uv run python -m project.commuting_extended [--long]
"""
from __future__ import annotations

import argparse
import itertools
import random
import sys
import time

import sympy as sp

from project.commuting import check, FAILURES, hs_rows, krull_dimension

FACES9 = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (0, 1, 1, 0), (1, 0, 0, 0), (1, 0, 0, 1),
          (1, 0, 1, 0), (0, 1, 0, 1)]
L9 = {f: k for k, f in enumerate(FACES9)}
CHI9 = [(1, a - b, c - d) for a, b, c, d in FACES9]
W9 = sp.symbols("w0:9")
_t = sp.symbols("t_sat")
VERTICAL = (7, 8)


def q9(z):
    return [1, 1, 1, 1, 1, z, z, 0, 0]


EVIDENT = {
    "I": [1, 0, 0, 0, 0, 0, 0, 1, 1],
    "S": [1, 1, 0, 0, 1, 1, 1, 0, 0],
    "S^-1": [1, 0, 1, 1, 1, 0, 1, 0, 0],
}


# ---------------------------------------------------------------------------
# Nine weights on hard-square rows
# ---------------------------------------------------------------------------
def strip9(bottom, top):
    N = len(bottom)
    c = [0] * 9
    for k in range(N):
        j = L9.get((top[k], top[(k + 1) % N], bottom[k], bottom[(k + 1) % N]))
        if j is None:
            return None
        c[j] += 1
    return tuple(c)


def strips9(N):
    R = hs_rows(N)
    st = {}
    for a in R:
        for b in R:
            e = strip9(a, b)
            if e is not None:
                st[(a, b)] = e
    return R, st


def mono9(weights, e):
    return sp.Mul(*[sp.sympify(weights[j]) ** e[j] for j in range(9) if e[j]])


def tm9(weights, N):
    R, st = strips9(N)
    return sp.Matrix(len(R), len(R), lambda i, j: mono9(weights, st[(R[i], R[j])]) if (R[i], R[j]) in st else 0)


def equations9(N, target):
    R, st = strips9(N)
    out, seen, done = [], set(), set()
    for x in R:
        for y in R:
            key = min((x[k:] + x[:k], y[k:] + y[:k]) for k in range(N))
            if key in done:
                continue
            done.add(key)
            p = sp.Integer(0)
            for m in R:
                if (x, m) in st and (m, y) in st:
                    p += mono9(target, st[(x, m)]) * mono9(W9, st[(m, y)]) - mono9(W9, st[(x, m)]) * mono9(target, st[(m, y)])
            p = sp.expand(p)
            if p != 0 and p not in seen:
                seen.add(p)
                out.append(p)
    return out


def walks9(S, N):
    st = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}
    M = sp.zeros(4, 4)
    for j in S:
        a, b, c, d = FACES9[j]
        M[st[(a, c)], st[(b, d)]] = 1
    return (M ** N).trace() != 0


def slice9(S):
    fixed = []
    for j in S:
        if sp.Matrix([CHI9[g] for g in fixed + [j]]).rank() > len(fixed):
            fixed.append(j)
    return fixed


def solve9(eqs, S):
    fixed = slice9(S)
    free = [j for j in S if j not in fixed]
    sub = {W9[j]: 0 for j in range(9) if j not in S}
    sub.update({W9[g]: 1 for g in fixed})
    E = [e for e in (sp.expand(e.subs(sub)) for e in eqs) if e != 0]
    unk = [W9[j] for j in free]
    if not E:
        return fixed, free, [], len(unk)
    if not unk:
        return None
    gb = sp.groebner(E + [_t * sp.Mul(*unk) - 1], _t, *unk, order="lex")
    if list(gb.exprs) == [1]:
        return None
    el = [g for g in gb.exprs if not g.has(_t)]
    if not el:
        return fixed, free, [], len(unk)
    gb2 = sp.groebner(el, *unk, order="lex")
    return fixed, free, list(gb2.exprs), krull_dimension(gb2, unk)


def parametrize(fixed, free, basis, S):
    """Points of the component as dictionaries, with the non-leading free weights kept symbolic."""
    base = {W9[g]: sp.Integer(1) for g in fixed}
    base.update({W9[j]: sp.Integer(0) for j in range(9) if j not in S})
    unk = [W9[j] for j in free]
    lead = []
    for g in basis:
        v = next(v for v in unk if g.has(v))
        if v not in lead:
            lead.append(v)
    sols = sp.solve(basis, lead, dict=True) if lead else [{}]
    out = []
    for s in sols:
        pt = dict(base)
        pt.update({v: v for v in unk})
        pt.update(s)
        out.append(pt)
    return out


def proportional(T, M):
    nzM = [(i, j) for i in range(M.rows) for j in range(M.cols) if M[i, j] != 0]
    if not nzM:
        return False
    c = sp.simplify(T[nzM[0]] / M[nzM[0]])
    if c == 0:
        return False
    return (T - c * M).applyfunc(sp.simplify) == sp.zeros(M.rows, M.cols)


def classify9(N, z):
    """{class name: [supports]} for all nine-weight W with T_N(W) != 0 commuting with T_N(Q0)."""
    target = q9(z)
    eqs = equations9(N, target)
    K = tm9(target, N)
    ref = {"K": K, **{n: tm9(f, N) for n, f in EVIDENT.items()}}
    ref["KS"] = K * ref["S"]
    ref["KS^-1"] = K * ref["S^-1"]
    out: dict[str, list] = {}
    for k in range(9, 0, -1):
        for S in itertools.combinations(range(9), k):
            if not walks9(S, N):
                continue
            res = solve9(eqs, S)
            if res is None:
                continue
            fixed, free, basis, dim = res
            for pt in parametrize(fixed, free, basis, S):
                T = tm9([pt[W9[j]] for j in range(9)], N)
                if T == sp.zeros(T.rows, T.cols):
                    continue
                name = next((n for n, M in ref.items() if proportional(T, M)), None)
                if name is None:
                    name = f"sporadic, rank {T.rank(simplify=True)}"
                out.setdefault(name, []).append(S)
    return out


# ---------------------------------------------------------------------------
# Sixteen faces: all binary rows
# ---------------------------------------------------------------------------
def all_rows(N):
    hs = hs_rows(N)
    rest = [r for r in itertools.product((0, 1), repeat=N) if r not in set(hs)]
    return hs, rest


def tm16(face, N):
    """face: dict (a,b,c,d) -> weight on all sixteen binary faces; rows: hard-square ones first."""
    hs, rest = all_rows(N)
    R = hs + rest

    def entry(x, y):
        v = sp.Integer(1)
        for k in range(N):
            v *= face.get((y[k], y[(k + 1) % N], x[k], x[(k + 1) % N]), 0)
            if v == 0:
                return sp.Integer(0)
        return v

    return sp.Matrix(len(R), len(R), lambda i, j: entry(R[i], R[j])), len(hs)


def compatibility_det(N):
    from sympy.polys.matrices import DomainMatrix

    R = hs_rows(N)
    n = len(R)
    M = DomainMatrix([[sp.ZZ(int(all(not (R[i][k] and R[j][k]) for k in range(N)))) for j in range(n)]
                      for i in range(n)], (n, n), sp.ZZ)
    return M.det()


def check_sixteen(N):
    z, u, v = sp.symbols("z u v")
    qhat = {f: q for f, q in zip(FACES9[:7], [1, 1, 1, 1, 1, z, z])}
    T, h = tm16(qhat, N)
    n = T.rows
    block = T[h:, :] == sp.zeros(n - h, n) and T[:, h:] == sp.zeros(n, n - h)
    full = (1, 1, 1, 1)
    Tu, _ = tm16({**qhat, full: u}, N)
    Tv, _ = tm16({**qhat, full: v}, N)
    E = sp.zeros(n, n)
    ones = n - 1  # the all-occupied row is the last one in itertools.product order
    E[ones, ones] = 1
    decoupled = (Tu - T - u ** N * E).expand() == sp.zeros(n, n)
    commute = ((T * Tu - Tu * T).expand() == sp.zeros(n, n)) and ((Tu * Tv - Tv * Tu).expand() == sp.zeros(n, n))
    return block, decoupled, commute


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
EXPECTED_CLASSES = {"K", "S", "S^-1", "I", "KS", "KS^-1"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--long", action="store_true", help="also z = -1/3 for N = 6")
    args = ap.parse_args(argv)
    t0 = time.time()

    print("== Nine weights on hard-square rows ==")
    for N in range(2, 8):
        z = sp.symbols("z")
        K = tm9(q9(z), N)
        ok = all(((K * M - M * K).expand() == sp.zeros(K.rows, K.rows)) for M in
                 [tm9(f, N) for f in EVIDENT.values()])
        ok &= tm9(EVIDENT["I"], N) == sp.eye(K.rows)
        ok &= tm9([1, 1, 1, 1, 0, z, z, z, 1], N) == K * tm9(EVIDENT["S"], N)
        ok &= tm9([1, 1, 1, 1, 1, z, 0, 1, z], N) == K * tm9(EVIDENT["S^-1"], N)
        check(f"N={N}: I, S, S^-1 commute with K; K S and K S^-1 are nine-weight transfer matrices", ok)
    zs = [sp.Integer(2)] + ([sp.Rational(-1, 3)] if args.long else [])
    for z in zs:
        for N in (4, 5, 6):
            cl = classify9(N, z)
            extra = sorted(c for c in cl if c not in EXPECTED_CLASSES)
            print(f"  z={z}, N={N}: classes {sorted(cl)}")
            if N < 6:
                check(f"z={z}, N={N}: the evident classes occur, and sporadic solutions exist at this size "
                      f"({len(extra)} classes outside the evident ones)", EXPECTED_CLASSES <= set(cl) and len(extra) > 0)
            else:
                check(f"z={z}, N=6: every nine-weight W with T_6(W) != 0 commuting with T_6(Q0) has T_6(W) a multiple "
                      f"of K, S, S^-1, I, K S or K S^-1", set(cl) == EXPECTED_CLASSES)
                check(f"z={z}, N=6: T_6(W) is a multiple of K only on supports without the vertical pair "
                      f"(w7 w8 = 0)", all(not set(VERTICAL) <= set(S) for S in cl["K"]))

    print("\n== Sixteen faces ==")
    dets = {N: compatibility_det(N) for N in range(2, 13)}
    check(f"det C_N = +-1 for 2 <= N <= 12 ({dets}), so T_N(Q) is invertible for z != 0",
          all(abs(d) == 1 for d in dets.values()))
    for N in range(2, 7):
        block, decoupled, commute = check_sixteen(N)
        check(f"N={N}: T_N(Q-hat) = T_N(Q) (+) 0; T_N(Q-hat + u[1111]) = T_N(Q-hat) + u^N E; "
              f"these commute with T_N(Q-hat) and with each other", block and decoupled and commute)

    print(f"\n{'All checks passed.' if not FAILURES else f'{len(FAILURES)} check(s) FAILED.'} "
          f"[{time.time() - t0:.0f} s]")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
