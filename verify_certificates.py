#!/usr/bin/env python3
"""Verification of the finite computations in
"A local Yang--Baxter obstruction for the hard-square model".

Every check is exact (sympy rational arithmetic or integer combinatorics), except the
calibration against the critical RSOS model, which is floating point and labelled as such.  Run with

    python verify_certificates.py            # all exact checks
    python verify_certificates.py --groebner # also the Groebner-basis cross-check
    python verify_certificates.py --classify # also the brute-force classification check

Requirements: Python >= 3.10, sympy >= 1.12, networkx >= 3.0.  Tested with Python 3.12,
sympy 1.14.0 and networkx 3.7.
"""
from __future__ import annotations

import argparse
import itertools
import sys

import networkx as nx
import sympy as sp

FAILURES: list[str] = []


def check(label: str, condition: bool) -> None:
    status = "PASS" if condition else "FAIL"
    print(f"[{status}] {label}")
    if not condition:
        FAILURES.append(label)


# ---------------------------------------------------------------------------
# Binary faces.  A face (a b ; c d) has top row a b and bottom row c d.
# ---------------------------------------------------------------------------
A, B, C, D, rho = sp.symbols("A B C D rho")
LABEL = {(0, 0, 0, 0): 0, (0, 0, 0, 1): 1, (0, 0, 1, 0): 2, (0, 1, 0, 0): 3,
         (0, 1, 1, 0): 4, (1, 0, 0, 0): 5, (1, 0, 0, 1): 6}
Q_BIN = [sp.Integer(1), D, C, B, B * C, A, A * D]


def binary_admissible(a, b, c, d):
    return a * b == 0 and c * d == 0 and a * c == 0 and b * d == 0


def w_bin(w, a, b, c, d):
    k = LABEL.get((a, b, c, d))
    return 0 if k is None else w[k]


def fyb_binary(P, Q, R, factor=1):
    """Return {boundary: LHS - factor*RHS} for the 64 binary boundaries (nonzero only)."""
    out = {}
    for bd in itertools.product((0, 1), repeat=6):
        a, b, c, d, e, f = bd
        lhs = sum(w_bin(P, f, g, a, b) * w_bin(Q, e, d, f, g) * w_bin(R, d, c, g, b) for g in (0, 1))
        rhs = sum(w_bin(R, e, g, f, a) * w_bin(Q, g, c, a, b) * w_bin(P, e, d, g, c) for g in (0, 1))
        ex = sp.expand(lhs - factor * rhs)
        if ex != 0:
            out[bd] = ex
    return out


def section_binary(run_groebner: bool) -> None:
    print("\n== Binary seven-weight faces ==")
    adm = sorted(f for f in itertools.product((0, 1), repeat=4) if binary_admissible(*f))
    check("exactly seven admissible binary faces, labelled as in (2.4)", adm == sorted(LABEL))
    p = sp.symbols("p0:7")
    r = sp.symbols("r0:7")
    E = fyb_binary(p, Q_BIN, r)
    check("18 nonzero boundary equations", len(E) == 18)
    listed = {
        (0, 0, 1, 0, 0, 0): p[0] * r[3] + D * p[3] * r[4] - B * p[1] * r[0],
        (0, 0, 1, 0, 1, 0): A * p[0] * r[3] + A * D * p[3] * r[4] - B * p[6] * r[5],
        (0, 0, 0, 1, 0, 0): B * p[0] * r[5] - p[3] * r[0] - A * p[4] * r[3],
        (0, 1, 0, 1, 0, 0): B * p[1] * r[6] - D * p[3] * r[0] - A * D * p[4] * r[3],
        (0, 0, 0, 1, 0, 1): B * C * p[5] * r[5] - p[3] * r[2] - A * p[4] * r[4],
        (0, 1, 0, 1, 0, 1): B * C * p[6] * r[6] - D * p[3] * r[2] - A * D * p[4] * r[4],
        (0, 0, 0, 0, 1, 0): A * D * p[3] * r[2] + A * p[0] * r[0] - p[5] * r[5],
        (1, 0, 0, 1, 0, 0): B * p[2] * r[5] - C * p[3] * r[1],
        (0, 0, 1, 0, 0, 1): C * p[5] * r[3] - B * p[1] * r[2],
        (0, 1, 0, 0, 0, 0): p[1] * r[1] - D * p[0] * r[0] - A * D * p[2] * r[3],
        (0, 0, 0, 0, 0, 1): C * p[5] * r[0] - A * p[2] * r[4] - p[0] * r[2],
        (0, 1, 0, 0, 0, 1): C * p[6] * r[1] - A * D * p[2] * r[4] - D * p[0] * r[2],
    }
    check("equations displayed in Section 3.1 and in the proof of the near-full support theorem agree with the enumeration",
          all(sp.expand(E[k] - v) == 0 for k, v in listed.items()))
    E1 = p[6] * r[5] - A * p[1] * r[0]
    E2 = p[1] * r[6] - D * p[0] * r[5]
    E3 = p[6] * r[6] - D * p[5] * r[5]
    E4 = p[5] * r[5] - A * p[0] * r[0] - A * D * p[3] * r[2]
    L = listed
    check("B*E1 = A*eq(0,0,1,0,0,0) - eq(0,0,1,0,1,0)", sp.expand(B * E1 - (A * L[(0, 0, 1, 0, 0, 0)] - L[(0, 0, 1, 0, 1, 0)])) == 0)
    check("B*E2 = eq(0,1,0,1,0,0) - D*eq(0,0,0,1,0,0)", sp.expand(B * E2 - (L[(0, 1, 0, 1, 0, 0)] - D * L[(0, 0, 0, 1, 0, 0)])) == 0)
    check("B*C*E3 = eq(0,1,0,1,0,1) - D*eq(0,0,0,1,0,1)", sp.expand(B * C * E3 - (L[(0, 1, 0, 1, 0, 1)] - D * L[(0, 0, 0, 1, 0, 1)])) == 0)
    check("E4 = -eq(0,0,0,0,1,0)", sp.expand(E4 + L[(0, 0, 0, 0, 1, 0)]) == 0)
    cert = sp.expand(r[6] * E1 + A * r[0] * E2 - r[5] * E3 - D * r[5] * E4 - A * D**2 * r[5] * p[3] * r[2])
    check("division-free certificate of the binary support theorem", cert == 0)
    # companion (remark after Theorem 3.3): the A_4 certificate with unstaggered weights and rho = 1
    E5 = D * p[3] * r[2] - A * p[2] * r[3]
    cert = sp.expand(r[6] * E1 + A * r[0] * E2 - r[5] * E3 - D * r[5] * E4 - A * D * r[5] * E5
                     - A**2 * D * r[5] * p[2] * r[3])
    check("companion: E5 = eq(0,0,0,0,0,0) and A^2*D*r5*p2*r3 = r6*E1 + A*r0*E2 - r5*E3 - D*r5*E4 - A*D*r5*E5 "
          "(the A_4 certificate for unstaggered weights)", sp.expand(E5 - E[(0, 0, 0, 0, 0, 0)]) == 0 and cert == 0)

    # scalar factor
    Er = fyb_binary(p, Q_BIN, r, factor=rho)
    F1 = p[6] * r[5] - A * p[1] * r[0]
    F2 = p[6] * r[1] - D * p[5] * r[0]
    F3 = A * p[1] * r[1] - rho * D * p[5] * r[5]
    ok = (sp.expand(rho * B * F1 - (A * Er[(0, 0, 1, 0, 0, 0)] - Er[(0, 0, 1, 0, 1, 0)])) == 0
          and sp.expand(C * F2 - (Er[(0, 1, 0, 0, 0, 1)] - D * Er[(0, 0, 0, 0, 0, 1)])) == 0
          and sp.expand(F3 - Er[(0, 1, 0, 0, 1, 0)]) == 0
          and sp.expand(r[0] * F3 + r[1] * F1 - r[5] * F2 - (1 - rho) * D * r[0] * p[5] * r[5]) == 0)
    check("scalar-factor identities (remark on orientations and scalar factors)", ok)

    # the two toric branches satisfy all equations
    p0, r0, r1, r5, p1, p5 = sp.symbols("p0_ r0_ r1_ r5_ p1_ p5_", nonzero=True)
    br_r = {p[0]: p0, p[1]: D * p0 * r0 / r1, p[2]: C * p0 * r1 / r0, p[3]: B * p0 * r5 / r0,
            p[4]: B * C * p0 * r1 * r5 / r0**2, p[5]: A * p0 * r0 / r5, p[6]: A * D * p0 * r0**2 / (r1 * r5),
            r[0]: r0, r[1]: r1, r[2]: 0, r[3]: 0, r[4]: r0**3 / (r1 * r5), r[5]: r5, r[6]: r1 * r5 / r0}
    br_p = {p[0]: p0, p[1]: p1, p[2]: 0, p[3]: 0, p[4]: p0**3 / (p1 * p5), p[5]: p5, p[6]: p1 * p5 / p0,
            r[0]: r0, r[1]: D * p0 * r0 / p1, r[2]: C * p5 * r0 / p0, r[3]: B * p1 * r0 / p0,
            r[4]: B * C * p1 * p5 * r0 / p0**2, r[5]: A * p0 * r0 / p5, r[6]: A * D * p0**2 * r0 / (p1 * p5)}
    check("toric branch r2=r3=0 satisfies all 18 equations",
          all(sp.simplify(e.subs(br_r)) == 0 for e in E.values()))
    check("toric branch p2=p3=0 satisfies all 18 equations",
          all(sp.simplify(e.subs(br_p)) == 0 for e in E.values()))
    # up to corner factors with product 1, the branches are (pure face, delta_{b,c}) and (delta_{b,c}, pure face)
    faces7 = sorted(LABEL, key=LABEL.get)

    def cornered(a_, b_, c_, d_, scale, identity=False):
        return [scale * a_**f[0] * b_**f[1] * c_**f[2] * d_**f[3] * (int(f[1] == f[2]) if identity else 1)
                for f in faces7]

    al, be, al2, be2 = r1 / r0, r5 / r0, p1 / p0, p5 / p0
    forms = [(br_r, r, cornered(be, 1 / al, 1 / be, al, r0, True)), (br_r, p, cornered(A / be, B * be, C * al, D / al, p0)),
             (br_p, p, cornered(be2, 1 / al2, 1 / be2, al2, p0, True)), (br_p, r, cornered(A / be2, B * al2, C * be2, D / al2, r0))]
    check("toric branches: up to corner factors with product 1, (P, R) is (pure face, delta_{b,c}) for r2=r3=0 "
          "and (delta_{b,c}, pure face) for p2=p3=0",
          all(sp.simplify(br[w[k]] - v) == 0 for br, w, vals in forms for k, v in enumerate(vals)))
    # uniqueness on each branch: solve the system on the support with the free parameters fixed,
    # and compare the unique solution with (3.19), respectively (3.20)
    Pn = list(sp.symbols("q0:7", nonzero=True))
    Rn = list(sp.symbols("s0:7", nonzero=True))
    to_n = {**dict(zip(p, Pn)), **dict(zip(r, Rn))}

    def is_branch(sol, br, free):
        expected = {to_n[k]: sp.sympify(v).subs(free) for k, v in br.items() if to_n[k] in sol[0]}
        return len(sol) == 1 and set(expected) == set(sol[0]) and all(sp.simplify(sol[0][k] - v) == 0
                                                                      for k, v in expected.items())

    Rz = list(Rn); Rz[2] = 0; Rz[3] = 0
    sol = sp.solve(list(fyb_binary(Pn, Q_BIN, Rz).values()),
                   [Pn[1], Pn[2], Pn[3], Pn[4], Pn[5], Pn[6], Rn[4], Rn[6]], dict=True)
    check("support {r2,r3}: the unique solution in terms of (p0,r0,r1,r5) is (3.19)",
          is_branch(sol, br_r, {p0: Pn[0], r0: Rn[0], r1: Rn[1], r5: Rn[5]}))
    Pz = list(Pn); Pz[2] = 0; Pz[3] = 0
    sol = sp.solve(list(fyb_binary(Pz, Q_BIN, Rn).values()),
                   [Rn[1], Rn[2], Rn[3], Rn[4], Rn[5], Rn[6], Pn[4], Pn[6]], dict=True)
    check("support {p2,p3}: the unique solution in terms of (p0,p1,p5,r0) is (3.20)",
          is_branch(sol, br_p, {p0: Pn[0], p1: Pn[1], p5: Pn[5], r0: Rn[0]}))

    if run_groebner:
        print("   (Groebner cross-check of the near-full support theorem over Q(A,B,C,D))")
        allvars = list(p) + list(r)
        t = sp.Symbol("t")
        eqs = list(E.values())
        survivors = []
        for k in (0, 1, 2):
            for Z in itertools.combinations(range(14), k):
                zero = {allvars[i]: 0 for i in Z}
                sys_ = [e for e in (sp.expand(x.subs(zero)) for x in eqs) if e != 0]
                if any(len(sp.Add.make_args(e)) == 1 for e in sys_):
                    continue
                nz = [v for i, v in enumerate(allvars) if i not in Z]
                G = sp.groebner(sys_ + [1 - t * sp.Mul(*nz)], *nz, t, order="grevlex",
                                domain=sp.QQ.frac_field(A, B, C, D))
                if not (len(G.exprs) == 1 and G.exprs[0] == 1):
                    survivors.append(tuple(str(allvars[i]) for i in Z))
        check("Groebner: the only zero sets of size <= 2 are {p2,p3} and {r2,r3}",
              sorted(survivors) == [("p2", "p3"), ("r2", "r3")])


# ---------------------------------------------------------------------------
# Blocked alphabet
# ---------------------------------------------------------------------------
S3 = ["0", "L", "R"]
BITS = {"0": (0, 0), "L": (0, 1), "R": (1, 0)}


def H_ok(s, t):
    return BITS[s][1] * BITS[t][0] == 0


def V_ok(s, t):
    return BITS[s][0] * BITS[t][0] == 0 and BITS[s][1] * BITS[t][1] == 0


def blocked_admissible(a, b, c, d):
    return H_ok(a, b) and H_ok(c, d) and V_ok(a, c) and V_ok(b, d)


def patch_ok(a, b, c, d):
    top = BITS[a] + BITS[b]
    bot = BITS[c] + BITS[d]
    return (all(top[i] * top[i + 1] == 0 for i in range(3)) and all(bot[i] * bot[i + 1] == 0 for i in range(3))
            and all(top[i] * bot[i] == 0 for i in range(4)))


FACES = [f for f in itertools.product(S3, repeat=4) if blocked_admissible(*f)]
NOCC = {"0": 0, "L": 1, "R": 1}


def Q_blk(a, b, c, d):
    if not blocked_admissible(a, b, c, d):
        return 0
    return A ** NOCC[a] * B ** NOCC[b] * C ** NOCC[c] * D ** NOCC[d]


def section_blocked() -> None:
    print("\n== Blocked alphabet ==")
    check("41 admissible blocked faces", len(FACES) == 41)
    check("blocked admissibility == 2x4 binary hard-square admissibility",
          all(blocked_admissible(*f) == patch_ok(*f) for f in itertools.product(S3, repeat=4)))

    print("\n== Blocked face Yang--Baxter equation ==")

    def run(support):
        Pv = {f: sp.Symbol("P_" + "".join(f)) for f in support}
        Rv = {f: sp.Symbol("R_" + "".join(f)) for f in support}
        eqs = {}
        for (a, b, c, d, e, f) in itertools.product(S3, repeat=6):
            lhs = sum(Pv.get((f, g, a, b), 0) * Q_blk(e, d, f, g) * Rv.get((d, c, g, b), 0) for g in S3)
            rhs = sum(Rv.get((e, g, f, a), 0) * Q_blk(g, c, a, b) * Pv.get((e, d, g, c), 0) for g in S3)
            ex = sp.expand(lhs - rhs)
            if ex != 0:
                eqs[(a, b, c, d, e, f)] = ex
        return eqs

    sing = run(FACES)[("0", "0", "L", "L", "0", "L")]
    check("one-equation obstruction: B C P(L0;00) R(LL;00) = 0",
          sp.expand(sing - B * C * sp.Symbol("P_L000") * sp.Symbol("R_LL00")) == 0)
    # weights also allowed on the faces whose four edges satisfy the vertical rule (remark after Corollary 4.2)
    FV = [f for f in itertools.product(S3, repeat=4)
          if V_ok(f[0], f[1]) and V_ok(f[2], f[3]) and V_ok(f[0], f[2]) and V_ok(f[1], f[3])]
    union = sorted(set(FACES) | set(FV))
    eqsV = run(union)
    check("weights on F and F_V (35 faces, 49 in the union): the equation of the one-equation obstruction is "
          "still B C P(L0;00) R(LL;00) = 0",
          len(FV) == 35 and len(union) == 49
          and sp.expand(eqsV[("0", "0", "L", "L", "0", "L")] - B * C * sp.Symbol("P_L000") * sp.Symbol("R_LL00")) == 0)


# ---------------------------------------------------------------------------
# The hard-hexagon alphabet A_4
# ---------------------------------------------------------------------------
def staggered_eqs(P, Q, R, rho_by_alpha=(1, 1)):
    """Staggered system: alpha-equations have P^{1-a} Q^a R^{1-a} on the left and
    rho_a * R^a Q^{1-a} P^a on the right (binary boundary data)."""
    out = {}
    for alpha in (0, 1):
        be = 1 - alpha
        for bd in itertools.product((0, 1), repeat=6):
            a, b, c, d, e, f = bd
            lhs = sum(w_bin(P[be], f, g, a, b) * w_bin(Q[alpha], e, d, f, g) * w_bin(R[be], d, c, g, b) for g in (0, 1))
            rhs = sum(w_bin(R[alpha], e, g, f, a) * w_bin(Q[be], g, c, a, b) * w_bin(P[alpha], e, d, g, c) for g in (0, 1))
            out[(alpha, bd)] = sp.expand(lhs - rho_by_alpha[alpha] * rhs)
    return out


def section_a4() -> None:
    import math
    print("\n== The hard-hexagon alphabet A_4 ==")
    lam = math.pi / 5
    Sh = {h: math.sin(h * lam) for h in range(0, 6)}
    adjh = lambda x, y: abs(x - y) == 1
    H4 = [1, 2, 3, 4]

    def Wrsos(a, b, c, d, u):
        if not (adjh(a, b) and adjh(c, d) and adjh(a, c) and adjh(b, d)):
            return 0.0
        val = math.sin(lam - u) if b == c else 0.0
        if a == d:
            val += math.sin(u) * math.sqrt(Sh[b] * Sh[c]) / Sh[a]
        return val

    worst = 0.0
    for (u, v) in [(0.37, 0.21), (0.11, 0.43), (-0.2, 0.55), (0.6, -0.15)]:
        for (a, b, c, d, e, f) in itertools.product(H4, repeat=6):
            lhs = sum(Wrsos(f, g, a, b, u) * Wrsos(e, d, f, g, u + v) * Wrsos(d, c, g, b, v) for g in H4)
            rhs = sum(Wrsos(e, g, f, a, v) * Wrsos(g, c, a, b, u + v) * Wrsos(e, d, g, c, u) for g in H4)
            worst = max(worst, abs(lhs - rhs))
    check(f"(float) critical RSOS A_4 weights satisfy the face YBE in the paper's convention (max residual {worst:.1e})",
          worst < 1e-15)
    # transport to occupations: phi_1(1)=1, phi_1(0)=3, phi_0(1)=4, phi_0(0)=2
    phi = {1: {1: 1, 0: 3}, 0: {1: 4, 0: 2}}
    faces = sorted(LABEL, key=LABEL.get)

    def types(u):
        return {t: [Wrsos(phi[t][x], phi[1 - t][y], phi[1 - t][uu], phi[t][w], u) for (x, y, uu, w) in faces]
                for t in (0, 1)}
    res = staggered_eqs(types(0.37), types(0.58), types(0.21))
    worst = max(abs(float(x)) for x in res.values())
    check(f"(float) the staggered form of the face YBE holds for the transported RSOS weights (max residual {worst:.1e})",
          worst < 1e-15)
    # exact: the transported critical weights do not depend on the type, and at u = -lambda they present the
    # hard-hexagon model at its critical activity (Example 2.13). Each weight is a combination of sin(lambda-u)
    # and sin(u), so equality at u = 0 and u = lambda gives equality for all u.
    lam_s = sp.pi / 5
    S_s = {h: sp.sin(h * lam_s) for h in range(1, 5)}

    def Wrsos_s(a, b, c, d, u):
        if not (adjh(a, b) and adjh(c, d) and adjh(a, c) and adjh(b, d)):
            return sp.Integer(0)
        return (sp.sin(lam_s - u) * int(b == c)
                + sp.sin(u) * int(a == d) * sp.sqrt(S_s[b] * S_s[c]) / S_s[a])

    def types_s(u):
        return {t: [Wrsos_s(phi[t][x], phi[1 - t][y], phi[1 - t][uu], phi[t][w], u) for (x, y, uu, w) in faces]
                for t in (0, 1)}
    same = all(sp.simplify(x - y) == 0 for u in (0, lam_s) for x, y in zip(*types_s(u).values()))
    wv = [sp.radsimp(sp.simplify(x)) for x in types_s(-lam_s)[1]]
    zv = sp.simplify(wv[5] * wv[3] * wv[2] * wv[1] / wv[0]**4)
    check("transported critical RSOS weights are the same for both types, and at u = -lambda give w6 = 0, "
          "w4*w0 = w3*w2 and z = (11+5*sqrt(5))/2: critical hard hexagons (Example 2.13)",
          same and sp.simplify(wv[0]) != 0 and sp.simplify(wv[6]) == 0
          and sp.simplify(wv[4] * wv[0] - wv[3] * wv[2]) == 0 and sp.simplify(zv - (11 + 5 * sp.sqrt(5)) / 2) == 0)
    # bijection between admissible A_4 faces of each type and admissible binary faces
    ok = True
    for t in (0, 1):
        img = [(phi[t][x], phi[1 - t][y], phi[1 - t][uu], phi[t][w]) for (x, y, uu, w) in faces]
        adm = [F for F in itertools.product(H4, repeat=4)
               if adjh(F[0], F[1]) and adjh(F[2], F[3]) and adjh(F[0], F[2]) and adjh(F[1], F[3]) and F[0] % 2 == t % 2]
        ok &= sorted(img) == sorted(adm)
    check("admissible A_4 faces of each type correspond bijectively to the 7 admissible binary faces", ok)

    # exact: Lemma 2.12 from the heights. Equation (2.2) on A_4, with P, Q, R given type by type through (2.5),
    # vanishes unless all six sides of the hexagon join adjacent heights, and is then the equation of (2.6)
    # indexed by the parity alpha of a and the occupations of the six heights. The equations of (2.6) not reached
    # in this way, those with two occupied neighbours on the hexagon, vanish identically.
    occ4 = {1: 1, 2: 0, 3: 0, 4: 1}

    def on_heights(Wt):
        def W(a, b, c, d):
            if not (adjh(a, b) and adjh(c, d) and adjh(a, c) and adjh(b, d)):
                return 0
            return Wt[a % 2][LABEL[(occ4[a], occ4[b], occ4[c], occ4[d])]]
        return W

    def fyb_a4(P, Q, R, factor):
        Ph, Qh, Rh = on_heights(P), on_heights(Q), on_heights(R)
        out = {}
        for bd in itertools.product(H4, repeat=6):
            a, b, c, d, e, f = bd
            lhs = sum(Ph(f, g, a, b) * Qh(e, d, f, g) * Rh(d, c, g, b) for g in H4)
            rhs = sum(Rh(e, g, f, a) * Qh(g, c, a, b) * Ph(e, d, g, c) for g in H4)
            out[bd] = sp.expand(lhs - factor * rhs)
        return out

    def hexagon(bd):
        return all(adjh(bd[k], bd[(k + 1) % 6]) for k in range(6))

    def staggered_index(bd):
        return bd[0] % 2, tuple(occ4[h] for h in bd)

    gen = {X: {t: list(sp.symbols(f"{X}{t}_0:7")) for t in (0, 1)} for X in "PQR"}
    EA = fyb_a4(gen["P"], gen["Q"], gen["R"], rho)
    ES = staggered_eqs(gen["P"], gen["Q"], gen["R"], rho_by_alpha=(rho, rho))
    hexagons = [bd for bd in EA if hexagon(bd)]
    image = {staggered_index(bd): EA[bd] for bd in hexagons}
    ok = all(ex == 0 for bd, ex in EA.items() if not hexagon(bd))
    ok &= len(image) == len(hexagons) and set(image) <= set(ES) and len(ES) == 128
    ok &= all(sp.expand(ES[k] - image.get(k, 0)) == 0 for k in ES)
    check("Lemma 2.12 (exact, from the heights 1-4): for arbitrary P, Q, R on A_4 the 4^6 equations (2.2) are "
          "the 128 equations (2.6), and vanish when a side of the hexagon joins non-adjacent heights", ok)

    # staggered system with pure hard-square middle faces
    A0, B0, C0, D0, A1, B1, C1, D1, r0_, r1_ = sp.symbols("A0 B0 C0 D0 A1 B1 C1 D1 rho0 rho1")
    Q = {0: [sp.Integer(1), D0, C0, B0, B0 * C0, A0, A0 * D0], 1: [sp.Integer(1), D1, C1, B1, B1 * C1, A1, A1 * D1]}
    p = {t: list(sp.symbols(f"p{t}_0:7")) for t in (0, 1)}
    r = {t: list(sp.symbols(f"r{t}_0:7")) for t in (0, 1)}
    E = staggered_eqs(p, Q, r, rho_by_alpha=(r0_, r1_))
    check("staggered system: 36 of the 128 equations are not identically zero",
          sum(1 for v in E.values() if v != 0) == 36)
    # the normalizations s_0, s_1 (Section 3.2): with middle face (s_0 Q^0, s_1 Q^1), an equation of parity alpha
    # is s_alpha times the equation of (2.6) with s_0 = s_1 = 1 and rho_alpha = rho s_{1-alpha} / s_alpha
    s = sp.symbols("s0 s1")
    Qs = {t: [s[t] * x for x in Q[t]] for t in (0, 1)}
    EAs = fyb_a4(p, Qs, r, rho)
    Ered = staggered_eqs(p, Q, r, rho_by_alpha=(rho * s[1] / s[0], rho * s[0] / s[1]))
    nonzero = [bd for bd, ex in EAs.items() if ex != 0]
    ok = len(nonzero) == 36 and sum(1 for bd in nonzero if bd[0] % 2 == 0) == 18 and all(map(hexagon, nonzero))
    ok &= all(sp.expand(EAs[bd] - s[bd[0] % 2] * Ered[staggered_index(bd)]) == 0 for bd in EAs if hexagon(bd))
    check("middle face (s_0 Q^0, s_1 Q^1) on A_4: 36 nonzero equations (2.2), 18 of each parity alpha, each s_alpha "
          "times the equation of (2.6) with s_0 = s_1 = 1 and rho_alpha = rho s_{1-alpha}/s_alpha", ok)
    e = lambda al, *bd: E[(al, bd)]
    listed = [
        ((0, 0, 0, 1, 0, 0, 0), p[1][0] * r[1][3] + D0 * p[1][3] * r[1][4] - r0_ * B1 * p[0][1] * r[0][0]),
        ((0, 0, 0, 1, 0, 1, 0), A0 * p[1][0] * r[1][3] + A0 * D0 * p[1][3] * r[1][4] - r0_ * B1 * p[0][6] * r[0][5]),
        ((1, 0, 0, 0, 1, 0, 0), B1 * p[0][0] * r[0][5] - r1_ * (p[1][3] * r[1][0] + A0 * p[1][4] * r[1][3])),
        ((1, 0, 1, 0, 1, 0, 0), B1 * p[0][1] * r[0][6] - r1_ * D0 * (p[1][3] * r[1][0] + A0 * p[1][4] * r[1][3])),
        ((1, 0, 0, 0, 1, 0, 1), B1 * C1 * p[0][5] * r[0][5] - r1_ * (p[1][3] * r[1][2] + A0 * p[1][4] * r[1][4])),
        ((1, 0, 1, 0, 1, 0, 1), B1 * C1 * p[0][6] * r[0][6] - r1_ * D0 * (p[1][3] * r[1][2] + A0 * p[1][4] * r[1][4])),
        ((0, 0, 0, 0, 0, 1, 0), A0 * p[1][0] * r[1][0] + A0 * D0 * p[1][3] * r[1][2] - r0_ * p[0][5] * r[0][5]),
        ((0, 0, 0, 0, 0, 0, 0), p[1][0] * r[1][0] + D0 * p[1][3] * r[1][2] - r0_ * (p[0][0] * r[0][0] + A1 * p[0][2] * r[0][3])),
    ]
    check("the eight staggered equations displayed in Section 3.2 agree with the enumeration",
          all(sp.expand(e(*k) - v) == 0 for k, v in listed))
    G1 = p[0][6] * r[0][5] - A0 * p[0][1] * r[0][0]
    G2 = p[0][1] * r[0][6] - D0 * p[0][0] * r[0][5]
    G3 = p[0][6] * r[0][6] - D0 * p[0][5] * r[0][5]
    G4 = r0_ * p[0][5] * r[0][5] - A0 * p[1][0] * r[1][0] - A0 * D0 * p[1][3] * r[1][2]
    G5 = p[1][0] * r[1][0] + D0 * p[1][3] * r[1][2] - r0_ * p[0][0] * r[0][0] - r0_ * A1 * p[0][2] * r[0][3]
    ok = sp.expand(r0_ * B1 * G1 - (A0 * e(0, 0, 0, 1, 0, 0, 0) - e(0, 0, 0, 1, 0, 1, 0))) == 0
    ok &= sp.expand(B1 * G2 - (e(1, 0, 1, 0, 1, 0, 0) - D0 * e(1, 0, 0, 0, 1, 0, 0))) == 0
    ok &= sp.expand(B1 * C1 * G3 - (e(1, 0, 1, 0, 1, 0, 1) - D0 * e(1, 0, 0, 0, 1, 0, 1))) == 0
    ok &= sp.expand(G4 + e(0, 0, 0, 0, 0, 1, 0)) == 0 and sp.expand(G5 - e(0, 0, 0, 0, 0, 0, 0)) == 0
    check("G1..G5 are the stated combinations of the staggered equations", ok)
    cert = sp.expand(r0_ * (r[0][6] * G1 + A0 * r[0][0] * G2 - r[0][5] * G3) - D0 * r[0][5] * G4
                     - A0 * D0 * r[0][5] * G5 - r0_ * A0 * A1 * D0 * r[0][5] * p[0][2] * r[0][3])
    check("A_4 certificate: rho0 A0 A1 D0 r5 p2 r3 (type 0) is the stated combination", cert == 0)
    # symmetry: exchanging the types maps the system to itself
    swap = {**{x: y for x, y in zip(p[0], p[1])}, **{x: y for x, y in zip(p[1], p[0])},
            **{x: y for x, y in zip(r[0], r[1])}, **{x: y for x, y in zip(r[1], r[0])},
            A0: A1, A1: A0, B0: B1, B1: B0, C0: C1, C1: C0, D0: D1, D1: D0, r0_: r1_, r1_: r0_}
    ok = all(sp.expand(E[(1 - al, bd)].xreplace(swap) - E[(al, bd)]) == 0 for (al, bd) in E)
    check("exchanging the two types (and alpha) maps the staggered system to itself", ok)
    # orientations (remark after Theorem 3.3): a symmetry g of the square, acting on the positions (a, b, c, d) of a
    # face, maps a staggered pure hard-square face W to the face weight W o g, which is again one; the types are
    # exchanged exactly when g moves the top-left corner to b or c.
    stag = on_heights(Qs)
    D4 = [(0, 1, 2, 3), (2, 0, 3, 1), (3, 2, 1, 0), (1, 3, 0, 2), (1, 0, 3, 2), (2, 3, 0, 1), (0, 2, 1, 3), (3, 1, 2, 0)]
    ok = True
    for g in D4:
        for t in (0, 1):
            w = [stag(*[F[i] for i in g]) for F in ((phi[t][x], phi[1 - t][y], phi[1 - t][uu], phi[t][ww])
                                                     for (x, y, uu, ww) in faces)]
            ok &= w[0] == s[t if g[0] in (0, 3) else 1 - t] and all(x != 0 for x in w)
            ok &= sp.expand(w[6] * w[0] - w[5] * w[1]) == 0 and sp.expand(w[4] * w[0] - w[3] * w[2]) == 0
    check("the eight symmetries of the square map staggered pure hard-square faces to staggered pure hard-square "
          "faces, exchanging the types when the top-left corner moves to b or c", ok)


def section_a4_tori() -> None:
    print("\n== Fibres of the occupation map, fibre sums, and transfer matrices on A_4 ==")

    # (1) Height configurations on A_4 versus hard squares on T_L: fibres are orbits of the reflection sigma.
    occ = {1: 1, 2: 0, 3: 0, 4: 1}
    sigma = {1: 4, 2: 3, 3: 2, 4: 1}
    for Lt, n_hs in ((2, 7), (4, 743)):
        sites = [(i, j) for i in range(Lt) for j in range(Lt)]
        nbrs = {(i, j): [((i + 1) % Lt, j), ((i - 1) % Lt, j), (i, (j + 1) % Lt), (i, (j - 1) % Lt)] for (i, j) in sites}
        idx = {s: k for k, s in enumerate(sites)}
        pairs = [(idx[s], idx[t]) for s in sites for t in nbrs[s][::2]]
        hs = [bits for bits in itertools.product((0, 1), repeat=Lt * Lt) if not any(bits[i] and bits[j] for i, j in pairs)]
        homs, h = [], {}

        def rec(k):
            if k == len(sites):
                homs.append(dict(h))
                return
            s = sites[k]
            for v in (1, 2, 3, 4):
                if all(t not in h or abs(v - h[t]) == 1 for t in nbrs[s]):
                    h[s] = v
                    rec(k + 1)
                    del h[s]

        rec(0)
        fib = {}
        for hh in homs:
            fib.setdefault(tuple(occ[hh[s]] for s in sites), []).append(hh)
        ok = len(hs) == n_hs and len(homs) == 2 * n_hs and set(fib) == set(hs)
        ok &= all(len(v) == 2 and v[1] == {s: sigma[x] for s, x in v[0].items()} for v in fib.values())
        check(f"T_{Lt}: {n_hs} hard-square configurations, {2 * n_hs} height configurations on A_4, "
              "every fibre an orbit of sigma of size 2", ok)

    # (2) Fibre sums: weights of the configurations used in the proof of Proposition 2.10, on T_6 and T_8.
    # Sites are (i, j) with j increasing downwards, so the bottom row of a face is at j + 1 and the paper's
    # e_2 is (0, -1) here. The second particle or pair is the first translated by (3, 1), the paper's 3e_1 - e_2,
    # which works on every T_L with L >= 6 even.
    w = {t: list(sp.symbols(f"w{t}_0:7")) for t in (0, 1)}

    def lift_weight(occset, odd_parity, Lt):
        """Face (i,j) has top row (i,j),(i+1,j) and bottom row (i,j+1),(i+1,j+1); its type is 1 iff
        its top-left site carries an odd height, i.e. iff (i+j) % 2 == odd_parity."""
        tot = sp.Integer(1)
        for i in range(Lt):
            for j in range(Lt):
                face = tuple(int(s in occset) for s in ((i, j), ((i + 1) % Lt, j), (i, (j + 1) % Lt), ((i + 1) % Lt, (j + 1) % Lt)))
                k = LABEL.get(face)
                if k is None:
                    return sp.Integer(0)
                tot *= w[1 if (i + j) % 2 == odd_parity else 0][k]
        return tot.subs({w[0][0]: 1, w[1][0]: 1})

    a1 = w[1][5] * w[0][3] * w[0][2] * w[1][1]
    a0 = w[0][5] * w[1][3] * w[1][2] * w[0][1]
    X = {t: w[t][6] / (w[t][5] * w[t][1]) for t in (0, 1)}
    Y = {t: w[t][4] / (w[t][3] * w[t][2]) for t in (0, 1)}
    cases = [
        ({(0, 0)}, a1, a0),
        ({(0, 0), (3, 1)}, a1 ** 2, a0 ** 2),
        ({(0, 0), (1, 1)}, a1 ** 2 * X[1], a0 ** 2 * X[0]),
        ({(0, 0), (1, 1), (3, 1), (4, 2)}, a1 ** 4 * X[1] ** 2, a0 ** 4 * X[0] ** 2),
        ({(2, 0), (1, 1)}, a1 ** 2 * Y[0], a0 ** 2 * Y[1]),
        ({(2, 0), (1, 1), (5, 1), (4, 2)}, a1 ** 4 * Y[0] ** 2, a0 ** 4 * Y[1] ** 2),
    ]
    ok = all(sp.simplify(lift_weight(S, 0, Lt) - e0) == 0 and sp.simplify(lift_weight(S, 1, Lt) - e1) == 0
             for Lt in (6, 8) for S, e0, e1 in cases)
    check("fibre sums on A_4 (L=6, 8): single particles, pairs and diagonal pairs have the weights used "
          "in the proof of the fibre-sum proposition", ok)

    # (3) Two-row transfer matrix of a staggered pure hard-square face versus (C_M Z_M)^2. T[par][x, y] is the weight
    # of the strip with bottom row x and top row y whose first height has parity par (Remark 2.11(1)).
    A0, B0, C0, A1, B1, C1, z = sp.symbols("A0 B0 C0 A1 B1 C1 z")
    D0 = z / (A0 * B1 * C1)
    D1 = z / (A1 * B0 * C0)
    Qt = {0: [1, D0, C0, B0, B0 * C0, A0, A0 * D0], 1: [1, D1, C1, B1, B1 * C1, A1, A1 * D1]}
    for M in (4, 6):
        rows = [r for r in itertools.product((0, 1), repeat=M) if all(not (r[k] and r[(k + 1) % M]) for k in range(M))]

        def strip(bot, top, par):
            tot = sp.Integer(1)
            for k in range(M):
                kk = LABEL.get((top[k], top[(k + 1) % M], bot[k], bot[(k + 1) % M]))
                if kk is None:
                    return sp.Integer(0)
                tot *= Qt[(k + par) % 2][kk]
            return tot

        n = len(rows)
        T = {par: sp.Matrix(n, n, lambda i, j: strip(rows[i], rows[j], par)) for par in (0, 1)}
        CZ = sp.Matrix(n, n, lambda i, j: z ** sum(rows[j]) * int(all(not (rows[i][k] and rows[j][k]) for k in range(M))))
        G1 = sp.diag(*[sp.prod([(Qt[(k + 1) % 2][2] * Qt[k % 2][1]) ** r[k] for k in range(M)]) for r in rows])
        diff = (T[1] * T[0] * G1 - G1 * CZ * CZ).applyfunc(sp.simplify)
        check(f"cylinder of circumference {M}: T_1 T_0 = G_1 (C_M Z_M)^2 G_1^(-1) for a staggered pure hard-square face",
              diff == sp.zeros(n, n))


def section_classification() -> None:
    print("\n== Classification of encodings without hidden heights (N <= 5 heights; slow) ==")
    Lt = 4
    nb = {}
    for i in range(Lt):
        for j in range(Lt):
            nb[(i, j)] = [((i + 1) % Lt, j), ((i - 1) % Lt, j), (i, (j + 1) % Lt), (i, (j - 1) % Lt)]
    sites = sorted(nb)
    even = frozenset(s for s in sites if (s[0] + s[1]) % 2 == 0)
    tests = [frozenset(), frozenset({(0, 0)}), even, frozenset(nb[(1, 1)]),
             frozenset({(0, 0), (0, 2)}), frozenset({(0, 0), (1, 1)}), frozenset({(0, 0), (2, 2)})]

    def count(occ, adj, U, O, cap):
        h, cnt = {}, 0

        def rec(k):
            nonlocal cnt
            if cnt > cap:
                return
            if k == len(sites):
                cnt += 1
                return
            s = sites[k]
            for x in (O if s in occ else U):
                if all(h.get(t) is None or adj[x][h[t]] for t in nb[s]):
                    h[s] = x
                    rec(k + 1)
                    del h[s]

        rec(0)
        return cnt

    classes = []
    examined = 0
    for N in range(2, 6):
        verts = list(range(N))
        for occmask in range(1, 2 ** N - 1):
            O = [v for v in verts if occmask >> v & 1]
            U = [v for v in verts if not occmask >> v & 1]
            cand = [(i, j) for i in verts for j in verts if i <= j and not (i in O and j in O)]
            for emask in range(2 ** len(cand)):
                Eg = [cand[k] for k in range(len(cand)) if emask >> k & 1]
                G = nx.Graph(); G.add_nodes_from(verts); G.add_edges_from(Eg)
                if not nx.is_connected(G):
                    continue
                examined += 1
                adj = [[0] * N for _ in verts]
                for i, j in Eg:
                    adj[i][j] = adj[j][i] = 1
                val, ok = None, True
                for cfg in tests:
                    c = count(cfg, adj, U, O, 20)
                    if c == 0 or c > 20 or (val is not None and c != val):
                        ok = False
                        break
                    val = c
                if ok:
                    for v in G:
                        G.nodes[v]["occ"] = v in O
                    if not any(nx.is_isomorphic(G, H, node_match=lambda a, b: a == b) for H in classes):
                        classes.append(G)
    # the two alphabets of Example 2.4, with their occupations and the loop at the empty height of the binary one
    binary = nx.Graph([(0, 0), (0, 1)])
    a4 = nx.path_graph([1, 2, 3, 4])
    for H, occupied in ((binary, {1}), (a4, {1, 4})):
        nx.set_node_attributes(H, {v: v in occupied for v in H}, "occ")
    iso = lambda G, H: nx.is_isomorphic(G, H, node_match=lambda a, b: a == b)
    check(f"{examined} labelled connected alphabets on 2 to 5 heights, with occupied and empty heights and no two "
          "occupied heights adjacent: exactly two isomorphism classes survive, the binary alphabet and A_4",
          len(classes) == 2 and any(iso(G, binary) for G in classes) and any(iso(G, a4) for G in classes))


# ---------------------------------------------------------------------------
# The binary alphabet with 1~1: all sixteen faces are admissible, and the middle face
# Q-hat is the pure hard-square face on the seven faces (2.4) and 0 on the other nine.
# ---------------------------------------------------------------------------
def section_full_binary() -> None:
    print("\n== Binary alphabet with 1~1 (all sixteen faces) ==")
    faces = list(itertools.product((0, 1), repeat=4))
    p = sp.symbols("p0:7")
    r = sp.symbols("r0:7")
    P = {f: p[LABEL[f]] if f in LABEL else sp.Symbol("p" + "".join(map(str, f))) for f in faces}
    R = {f: r[LABEL[f]] if f in LABEL else sp.Symbol("r" + "".join(map(str, f))) for f in faces}
    Qh = {f: Q_BIN[LABEL[f]] if f in LABEL else sp.Integer(0) for f in faces}

    def eq(bd, Qm=Qh):
        a, b, c, d, e, f = bd
        lhs = sum(P[(f, g, a, b)] * Qm[(e, d, f, g)] * R[(d, c, g, b)] for g in (0, 1))
        rhs = sum(R[(e, g, f, a)] * Qm[(g, c, a, b)] * P[(e, d, g, c)] for g in (0, 1))
        return sp.expand(lhs - rho * rhs)

    x, v = P[(0, 1, 0, 1)], P[(1, 0, 1, 0)]
    y, u = R[(0, 0, 1, 1)], R[(1, 1, 0, 0)]
    golden = fyb_binary(p, Q_BIN, r, factor=rho)
    six = [(0, 0, 1, 0, 0, 0), (0, 0, 1, 0, 1, 0), (0, 0, 0, 1, 0, 0), (0, 1, 0, 1, 0, 0), (0, 0, 0, 1, 0, 1),
           (0, 1, 0, 1, 0, 1)]
    listed = {
        (0, 0, 0, 0, 1, 0): A * D * p[3] * r[2] + A * p[0] * r[0] - rho * p[5] * r[5] - rho * A * v * u,
        (0, 0, 0, 0, 0, 0): D * p[3] * r[2] + p[0] * r[0] - rho * p[0] * r[0] - rho * A * p[2] * r[3],
        (1, 0, 1, 1, 0, 1): B * C * v * u - rho * B * C * x * y,
        (1, 0, 1, 1, 0, 0): B * p[2] * u - rho * B * C * x * r[1],
        (1, 0, 1, 0, 0, 1): C * v * r[3] - rho * B * C * p[1] * y,
        (1, 0, 1, 0, 0, 0): p[2] * r[3] + D * p[4] * r[4] - rho * B * C * p[1] * r[1],
    }
    ok = (all(sp.expand(eq(bd) - golden[bd]) == 0 for bd in six)
          and sp.expand(eq((0, 0, 0, 0, 1, 0)) - golden[(0, 0, 0, 0, 1, 0)] + rho * A * v * u) == 0
          and all(sp.expand(eq(bd) - val) == 0 for bd, val in listed.items()))
    check("binary alphabet with 1~1: (3.1)-(3.6) unchanged, (3.7) gains -A*rho*P(1010)*R(1100), and the other "
          "equations of Proposition 4.3 agree with the enumeration", ok)

    e1, e2, e3, e4, e5, e6 = (eq(bd) for bd in six)
    e7, f0 = eq((0, 0, 0, 0, 1, 0)), eq((0, 0, 0, 0, 0, 0))
    g, h2, h3, h4 = (eq(bd) for bd in [(1, 0, 1, 1, 0, 1), (1, 0, 1, 1, 0, 0), (1, 0, 1, 0, 0, 1), (1, 0, 1, 0, 0, 0)])
    S1 = C * r[6] * (A * e1 - e2) + rho * A * C * r[0] * (e4 - D * e3) - rho * r[5] * (e6 - D * e5)
    Y = B * C * D * r[5] * (A * f0 - e7) - S1 - rho * A * D * r[5] * g
    X = rho * B * C * x * y * h4 - (C * v * r[3] * h2 + rho * B * C * x * r[1] * h3 - p[2] * r[3] * g)
    cert = sp.expand(rho * A * r[5] * X - p[4] * r[4] * Y - rho * A**2 * B * C * D * p[2] * p[4] * r[3] * r[4] * r[5])
    check("binary alphabet with 1~1: division-free certificate rho*A^2*B*C*D*p2*p4*r3*r4*r5 (Proposition 4.3)", cert == 0)

    # the proof uses only the zeros of Q-hat on five faces (remark after Proposition 4.3): the twelve equations do
    # not change when Q-hat is arbitrary on the other four faces with two adjacent occupied sites, and each of the
    # five zeros enters one of them
    twelve = six + [(0, 0, 0, 0, 1, 0)] + [bd for bd in listed if bd != (0, 0, 0, 0, 1, 0)]
    five = [(1, 1, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 0)]
    others = [f for f in faces if f not in LABEL and f not in five]
    free = {f: sp.Symbol("q" + "".join(map(str, f))) for f in faces}
    ok = len(others) == 4 and all(sp.expand(eq(bd, {**Qh, **{f: free[f] for f in others}}) - eq(bd)) == 0 for bd in twelve)
    ok &= all(any(sp.expand(eq(bd, {**Qh, f: free[f]}) - eq(bd)) != 0 for bd in twelve) for f in five)
    check("binary alphabet with 1~1: the twelve equations use only the seven hard-square values of Q-hat and its "
          "zeros on 1100, 0011, 0101, 0111, 1110", ok)
    # the pure hard-square face extended by zero only on the four faces with a = c = 1 and on 1100 kills exactly the
    # configurations with two adjacent occupied sites, here on the tori T_2, T_3, T_4
    alt = {f for f in faces if f[0] == f[2] == 1} | {(1, 1, 0, 0)}
    ok = len(alt) == 5
    for Lt in (2, 3, 4):
        for bits in itertools.product((0, 1), repeat=Lt * Lt):
            x = lambda i, j: bits[(i % Lt) * Lt + j % Lt]
            killed = any((x(i, j), x(i + 1, j), x(i, j + 1), x(i + 1, j + 1)) in alt for i in range(Lt) for j in range(Lt))
            adjacent = any(x(i, j) and (x(i + 1, j) or x(i, j + 1)) for i in range(Lt) for j in range(Lt))
            ok &= killed == adjacent
    check("binary alphabet with 1~1: zeros on the faces with a = c = 1 and on 1100 kill exactly the configurations "
          "with two adjacent occupied sites (T_2, T_3, T_4)", ok)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--groebner", action="store_true", help="also run the Groebner cross-check of the near-full support theorem")
    parser.add_argument("--classify", action="store_true",
                        help="also run the brute-force check of the classification of encodings (slow)")
    args = parser.parse_args()
    section_binary(args.groebner)
    section_a4()
    section_a4_tori()
    section_blocked()
    section_full_binary()
    if args.classify:
        section_classification()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED:")
        for f in FAILURES:
            print("  -", f)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
