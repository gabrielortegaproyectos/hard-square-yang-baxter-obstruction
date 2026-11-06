"""Small-size commutation with the hard-square transfer matrix (src/project/commuting.py)."""
from fractions import Fraction

import pytest
import sympy as sp

from project import commuting as cm


@pytest.mark.parametrize("N", [2, 3, 4, 5])
def test_reduction_gauges_and_shift(N):
    assert cm.check_reduction(N)


def test_row_counts_are_lucas_numbers():
    assert [len(cm.hs_rows(N)) for N in range(2, 9)] == [3, 4, 7, 11, 18, 29, 47]


def test_evident_solutions_commute():
    N = 6
    K = cm.transfer_matrix(cm.Q0, N)
    for face in (cm.SHIFT, cm.SHIFT_INV, cm.Q0):
        T = cm.transfer_matrix(face, N)
        assert (K * T - T * K).expand() == sp.zeros(*K.shape)


def test_n4_set_is_the_evident_one_but_q0_is_double():
    comp = cm.components(4, cm.Q0, symbolic_z=True)
    assert set(comp) == set(cm.EXPECTED)
    gens = [cm.W[2], cm.W[4], cm.W[5], cm.W[6], cm.z]
    gb = sp.groebner(comp[cm.FULL][2], *gens, order="lex")
    assert gb.reduce(sp.expand((cm.W[4] - 1) ** 2))[1] == 0
    assert gb.reduce(cm.W[4] - 1)[1] != 0


def test_n5_only_evident_solutions_for_every_z():
    assert cm.matches_expected(cm.components(5, cm.Q0, symbolic_z=True))


def test_n5_tangent_trivial_for_every_z():
    assert cm.tangent_invariant_factors(5) == [1, 1, cm.z ** 2, cm.z ** 2]


@pytest.mark.parametrize("zval", [Fraction(2), Fraction(-1), Fraction(-1, 4)])
def test_tangent_dimensions(zval):
    assert [cm.tangent_dimension(N, zval) for N in range(2, 8)] == [5, 5, 4, 3, 3, 3]


def test_hard_hexagon_control_finds_baxter_curve():
    hh = [sp.sympify(x).subs(cm.z, 2) for x in cm.HARD_HEXAGONS]
    fixed, free, basis, dim = cm.components(4, hh, symbolic_z=False)[cm.FULL]
    assert dim == 1
    assert cm.on_baxter_surface(fixed, free, basis)


@pytest.mark.slow
def test_n6_only_evident_solutions_for_every_z():
    assert cm.matches_expected(cm.components(6, cm.Q0, symbolic_z=True))


# ---------------------------------------------------------------------------
# Independent cross-check and extension beyond the seven weights
# ---------------------------------------------------------------------------
from project import commuting_crosscheck as cc  # noqa: E402
from project import commuting_extended as ce  # noqa: E402


def test_crosscheck_n5_new_activity():
    supports, classes, reduced, _ = cc.check_size(5, sp.Integer(3))
    assert classes == {(0, 1, 2, 3, 4, 5, 6): {"Q0"}, (0, 1, 4, 5, 6): {"S"}, (0, 2, 3, 4, 6): {"S^-1"}}
    assert all(reduced.values())


def test_crosscheck_jacobian_ranks():
    assert {N: cc.tangent_rank(N, sp.Rational(-2, 7)) for N in (4, 5)} == {4: 3, 5: 4}


@pytest.mark.slow
def test_crosscheck_full():
    assert cc.main() == 0


@pytest.mark.parametrize("N", [3, 4, 5])
def test_nine_weight_evident_partners(N):
    z = sp.symbols("z")
    K = ce.tm9(ce.q9(z), N)
    assert ce.tm9(ce.EVIDENT["I"], N) == sp.eye(K.rows)
    assert ce.tm9([1, 1, 1, 1, 0, z, z, z, 1], N) == K * ce.tm9(ce.EVIDENT["S"], N)


def test_nine_weights_n4_has_sporadic_solution():
    cl = ce.classify9(4, sp.Integer(2))
    assert ce.EXPECTED_CLASSES <= set(cl) and len(set(cl) - ce.EXPECTED_CLASSES) == 1


@pytest.mark.slow
def test_nine_weights_n6_only_evident():
    cl = ce.classify9(6, sp.Integer(2))
    assert set(cl) == ce.EXPECTED_CLASSES
    assert all(not set(ce.VERTICAL) <= set(S) for S in cl["K"])


@pytest.mark.parametrize("N", [3, 4])
def test_sixteen_faces_block_and_decoupled_family(N):
    assert all(ce.check_sixteen(N))


def test_compatibility_determinants():
    assert all(abs(ce.compatibility_det(N)) == 1 for N in range(2, 10))
