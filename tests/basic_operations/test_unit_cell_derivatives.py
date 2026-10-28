"""Analytical derivatives of unit_cell.py over the unit cell parameters
against central finite differences, on a triclinic cell (with 90° angles
many wrong terms vanish)."""
import inspect

import numpy
import pytest

from cryspy.A_functions_base import unit_cell

UNIT_CELL_PARAMETERS = numpy.array([5.1, 6.3, 7.7, *numpy.radians([81.0, 97.5, 104.0])])
INDEX_HKL = numpy.array([[1, 0, 2], [0, 1, 1], [0, 0, -2]])
STEP = 1e-6

# derivatives found wrong on 2026-10-28; fix them and remove from this list
WRONG_DERIVATIVES = {
    "calc_eq_ccs_by_unit_cell_parameters", "calc_m_b_by_unit_cell_parameters",
    "calc_m_b_norm_by_unit_cell_parameters", "calc_m_g_norm_by_unit_cell_parameters",
    "calc_m_inv_b_by_unit_cell_parameters", "calc_m_inv_b_norm_by_unit_cell_parameters",
    "calc_m_inv_b_norm_m_by_unit_cell_parameters", "calc_m_inv_m_b_norm_by_unit_cell_parameters",
    "calc_m_inv_m_by_unit_cell_parameters", "calc_m_inv_m_norm_by_unit_cell_parameters",
    "calc_m_m_by_unit_cell_parameters", "calc_m_m_norm_by_unit_cell_parameters",
    "calc_m_reciprocal_g_by_unit_cell_parameters", "calc_q_ccs_by_unit_cell_parameters",
    "calc_reciprocal_by_unit_cell_parameters",
}

NAMES = sorted(name for name in dir(unit_cell)
               if name.startswith("calc_") and name.endswith("_by_unit_cell_parameters"))


def take_value(function, unit_cell_parameters, **kwargs):
    if "index_hkl" in inspect.signature(function).parameters:
        return function(INDEX_HKL, unit_cell_parameters, **kwargs)
    return function(unit_cell_parameters, **kwargs)


@pytest.mark.parametrize("name", [
    pytest.param(name, marks=pytest.mark.xfail(strict=True, reason="wrong derivative"))
    if name in WRONG_DERIVATIVES else name for name in NAMES])
def test_derivative_over_unit_cell_parameters(name):
    function = getattr(unit_cell, name)
    value, dder = take_value(function, UNIT_CELL_PARAMETERS, flag_unit_cell_parameters=True)
    l_numerical = []
    for ind in range(6):
        shift = numpy.zeros(6)
        shift[ind] = STEP
        value_plus = take_value(function, UNIT_CELL_PARAMETERS + shift)[0]
        value_minus = take_value(function, UNIT_CELL_PARAMETERS - shift)[0]
        l_numerical.append((numpy.asarray(value_plus) - numpy.asarray(value_minus)) / (2. * STEP))
    analytical = numpy.asarray(dder["unit_cell_parameters"])
    value = numpy.asarray(value)
    # the parameter axis is first for scalars and lists of reflections,
    # last for matrices stored as components
    numerical = numpy.stack(l_numerical, axis=-1 if analytical.shape == value.shape + (6, ) else 0)
    assert numerical.shape == analytical.shape
    assert numpy.allclose(analytical, numerical, rtol=1e-5, atol=1e-7)
