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
