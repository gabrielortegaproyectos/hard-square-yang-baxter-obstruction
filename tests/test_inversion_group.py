"""The inversion group of binary face weights (src/project/inversion_group.py)."""
from fractions import Fraction

import pytest
import sympy as sp

from project import inversion_group as ig


def test_matrices_and_explicit_inverses():
    assert ig.check_blocks()
    assert ig.check_explicit_inverses()
    assert ig.check_determinants()


def test_action_on_coordinates_and_involutions():
    assert ig.check_coordinate_action()
    assert ig.check_involutions()


def test_the_equality_test_is_not_vacuous():
    point = ig.coordinates(ig.W)
    assert not ig.same(ig.coordinates(ig.invert(ig.W, "bc")), ig.i2(point))
    assert not ig.same(ig.i2(ig.i1(point)), point)


def test_period_five_and_order_five():
    assert ig.check_lyness_period()
    assert ig.check_step()
    assert ig.check_order_five()


def test_dihedral_group_of_order_ten():
    assert len(ig.group_elements()) == 10


def test_baxter_surface_is_preserved():
    assert ig.check_baxter_surface()


def test_exact_rational_point():
    """One fixed rational face weight, followed through the ten elements of the group."""
    weights = [Fraction(n, d) for n, d in ((3, 2), (-5, 7), (4, 3), (7, 5), (-2, 9), (11, 4), (5, 6))]
    start = ig.coordinates([sp.Rational(w.numerator, w.denominator) for w in weights])
    current = [sp.Rational(w.numerator, w.denominator) for w in weights]
    seen = [start]
    for _ in range(5):
        for reading in ("bc", "ad"):
            current = ig.invert(current, reading)
            seen.append(ig.coordinates(current))
    assert seen[-1] == start
    assert len(set(seen[:-1])) == 10


def test_a4_single_steps_and_period_products():
    assert ig.check_a4_single_steps()
    assert ig.check_a4_period_products()


def test_a4_fifth_power_rescales_the_sublattice_activities():
    assert ig.check_a4_fifth_power()


@pytest.mark.slow
def test_fifth_power_is_a_gauge_transformation_of_the_weights():
    assert ig.check_weights_up_to_gauge()


@pytest.mark.slow
def test_module_runs_clean():
    assert ig.main() == 0
