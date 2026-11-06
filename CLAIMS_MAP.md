# Claims map

Each finite statement of the paper *A local Yang–Baxter obstruction for the
hard-square model* that the script checks.

- **Numbering** follows the paper after its focus revision (2026-11-01).
  Results are also named, because numbers may change between versions.
- **Exactness.** Checks are exact (rational or integer arithmetic) unless
  marked *(float)*.
- **Role of the script.** Every proof in the paper can be read without the
  script, which only guards against arithmetic slips.

Run `python verify_certificates.py --groebner --classify`. The section
headers below match the headers printed by the script.

## Binary seven-weight faces

| Paper statement | Check |
|---|---|
| (2.4) the seven admissible binary faces | exactly seven admissible binary faces |
| (3.1)–(3.7), equations with a pure hard-square middle face | 18 nonzero boundary equations; displayed equations agree with the enumeration |
| Theorem 3.1 (Binary certificate), identity (3.8) | the combinations giving `E1`–`E4`; division-free certificate |
| Remark after Theorem 3.3, companion `r5 p2 r3 = 0` for unstaggered weights | `E5 = eq(0,0,0,0,0,0)` and the identity for `A^2 D r5 p2 r3` |
| Remark 3.2 (Orientations and scalar factors) | scalar-factor identities |
| Theorem 3.7 (Near-full support), branches (3.19)–(3.20) | both toric branches satisfy all 18 equations; on each support the unique solution is (3.19), respectively (3.20); up to corner factors with product 1 they are (pure face, `delta_{b,c}`) and (`delta_{b,c}`, pure face) |
| Theorem 3.7, cross-check | `--groebner`: the only zero sets of size at most 2 are `{p2,p3}` and `{r2,r3}` |

## The hard-hexagon alphabet `A_4`

| Paper statement | Check |
|---|---|
| Example 2.13 (Calibration) | *(float)* critical RSOS weights satisfy the face YBE in the paper's convention; exact: transported to occupations they do not depend on the type, and at `u = -lambda` they give critical hard hexagons, `z = (11+5 sqrt 5)/2` |
| Lemma 2.12 (Staggered form) | exact, from the heights 1–4: for arbitrary `P, Q, R` on `A_4` the equations (2.2) are the equations (2.6); bijection of `A_4` faces of each type with the binary faces; *(float)* staggered form holds for the transported RSOS weights |
| Section 3.2, normalizations `s_0, s_1` | with middle face `(s_0 Q^0, s_1 Q^1)`, each of the 36 nonzero equations (2.2) is `s_alpha` times the equation of (2.6) with `rho_alpha = rho s_{1-alpha}/s_alpha` |
| (3.9)–(3.16), staggered equations | 36 of 128 equations nonzero; displayed equations agree with the enumeration |
| Theorem 3.3 (Certificate on `A_4`), identity (3.17) | `G1`–`G5` are the stated combinations; the certificate identity; symmetry under exchange of types |
| Remark after Theorem 3.3, orientations | the eight symmetries of the square map staggered pure hard-square faces to staggered pure hard-square faces, exchanging the types when the top-left corner moves to `b` or `c` |

## Fibres of the occupation map, fibre sums, and transfer matrices on `A_4`

| Paper statement | Check |
|---|---|
| Example 2.4 (2) (the two lifts) | `T_2`, `T_4`: 7 and 743 hard-square configurations, 14 and 1486 height configurations on `A_4`; every fibre is a `sigma`-orbit of size 2 |
| Proposition 2.10 (Fibre sums) | weights on `T_6` and `T_8` of the configurations used in the proof |
| Remark 2.11 (1) (transfer matrices) | `T_1 T_0 G_1 = G_1 (C_M Z_M)^2` for circumference 4 and 6, with symbolic corner factors |

## Blocked alphabet

| Paper statement | Check |
|---|---|
| Section 4.1 (admissible blocked faces) | 41 faces; blocked admissibility is admissibility of the 2×4 binary pattern |

## Blocked face Yang–Baxter equation

| Paper statement | Check |
|---|---|
| Theorem 4.1 (One-equation obstruction) | `B C P(L0;00) R(LL;00) = 0` |
| Remark after Corollary 4.2, weights on `F ∪ F_V` | 35 faces in `F_V`, 49 in the union; the equation of Theorem 4.1 is unchanged |

## Binary alphabet with 1~1 (all sixteen faces)

| Paper statement | Check |
|---|---|
| Proposition 4.3 (Binary alphabet with `1~1`), the twelve equations of its proof | (3.1)–(3.6) unchanged; (3.7) gains `-A rho P(1010) R(1100)`; the other five equations agree with the enumeration |
| Proposition 4.3, certificate | division-free identity for `rho A^2 B C D p2 p4 r3 r4 r5` |
| Remark after Proposition 4.3 | the twelve equations use only the seven hard-square values of `Q-hat` and its zeros on five faces; the alternative zero set (faces with `a = c = 1` and `1100`) kills exactly the configurations with two adjacent occupied sites on `T_2`, `T_3`, `T_4` |

## Classification (illustration)

| Paper statement | Check |
|---|---|
| Proposition 2.5 (Classification) | `--classify`: of the 86 828 labelled connected alphabets on 2 to 5 heights, with occupied and empty heights and no two occupied heights adjacent, exactly two isomorphism classes survive; they are compared by isomorphism, occupations included, with the binary alphabet and `A_4`. This illustrates the proposition, which is proved in the text. |

## Exploratory: transfer matrices commuting with hard squares (not in the paper)

These checks are printed by `uv run python -m project.commuting` (add `--long` for the
rows marked so), not by `verify_certificates.py`. They answer the question of Section 5.2
of the J. Phys. A version; see
[`docs/simulations/commuting-small-n.md`](docs/simulations/commuting-small-n.md).

| Statement | Check |
|---|---|
| `T_N(Q) = D_CD C_N D_AB`; only `z` matters; corner gauges; the shift commutes | `N=2..6: T_N(Q) = D_CD C_N D_AB, conjugation to Q0, corner gauges, shift commutes` |
| For `N = 2, 3` the commuting faces form positive-dimensional sets | `N=2`, `N=3`: positive-dimensional solutions on the full support besides Q0 |
| For `N = 4` only multiples of `T_N(Q)`, `S`, `S^{-1}`; `Q` is a double point | `N=4`: solutions only on the supports of Q0, delta_{b,c} and delta_{a,d}; the only point is Q0, not reduced |
| For `N = 5, 6` (and `7`, `--long`) exactly those three, reduced, for every `z != 0` | `N=5`, `N=6`, `N=7`: exactly W ~ Q0, delta_{b,c}, delta_{a,d} |
| Tangent space at `Q` is trivial for every `z != 0` with `N = 5, 6` | invariant factors over Q[z] are `[1, 1, z**2, z**2]` |
| Tangent space at `Q` is trivial for `5 <= N <= 9` (`12` with `--long`) at five activities | `z=...: tangent dimensions ...` |
| Positive control: hard hexagons have a curve of commuting partners on Baxter's surface | the full support carries a curve of commuting partners, on Baxter's surface F = 0; the curve is the same for every N computed |
| Independent cross-check, `python -m project.commuting_crosscheck` | `z=..., N=4,5,6`: only the classes of T_N(Q0), S, S^-1; Q0 double at N=4, reduced from N=5; Jacobian ranks |
| A_4, `python -m project.commuting_a4`: every presentation reduces to `Q0 (+) Q0`; six trivial symmetries | `M=4: the six trivial symmetries multiply T_M by lam^M, ...` |
| A_4: near a presentation, only the trivial directions commute, for `M = 6, 8` (`10` with `--long`) at five activities | `M=...: tangent dimension ... (only the six trivial directions)` |
| A_4: the same for every `z != 0` with `M = 6` | `over Q[z], the gcd of the 8 x 8 minors ... divides z**...` |
| A_4 control: hard hexagons on both types have one more, unstaggered, direction | `hard hexagons ... one more direction, and it is unstaggered` |
