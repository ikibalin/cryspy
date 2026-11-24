"""Small bugs of the roadmap section "Errors": loops with a longer category
name, refine_all_variables with undefined attributes, 2D interpolation
without scipy.interpolate.interp2d (removed in scipy 1.14)."""
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.function_1_matrices import calc_bilinear_interpolation_on_grid
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_pd2d_background import Pd2dBackground
from cryspy.procedure_rhochi.rhochi_pd2d import calc_background

DIR_FUNCTIONAL = os.path.join(os.path.dirname(__file__), "..", "functional")

S_ANISO_FIRST = """
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
Fe 0.01

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe Fe 0.125 0.125 0.125
"""


def test_loop_skips_category_with_longer_name():
    # atom_site_aniso starts with atom_site: it was taken for atom_site and
    # failed with UnboundLocalError
    atom_site = AtomSiteL.from_cif(S_ANISO_FIRST)
    assert [item.label for item in atom_site.items] == ["Fe"]
    assert atom_site.items[0].fract_x == pytest.approx(0.125)


@pytest.mark.parametrize("f_name", [
    os.path.join("Ho2Ti2O7_5K_1T_powder_1d_test", "powder1_calc.rcif"),
    os.path.join("tof", "tof.rcif")])
def test_calc_rcif_loads(f_name):
    rcif = cryspy.file_to_globaln(os.path.join(DIR_FUNCTIONAL, f_name))
    assert len(rcif.items) > 0


def test_refine_all_variables_only_defined_attributes():
    rcif = cryspy.file_to_globaln(os.path.join(
        DIR_FUNCTIONAL, "PbSO4_unpol_powder_test", "main.rcif"))
    crystal = [block for block in rcif.items if isinstance(block, cryspy.Crystal)][0]
    crystal.atom_site.refine_all_variables()
    l_name = rcif.get_variable_names()
    assert any("fract_x" in str(name) for name in l_name)
    assert not any("u_iso_or_equiv" in str(name) for name in l_name)
    for name in l_name:
        rcif.get_variable_by_name(name)


def test_bilinear_interpolation_on_grid():
    x_grid = numpy.array([10., 0., 20.])  # any order
    y_grid = numpy.array([-5., 5.])
    z_grid = 2. * x_grid[:, None] - 3. * y_grid[None, :] + 1.
    x = numpy.array([0., 3., 17.5, 25.])
    y = numpy.array([-7., 0., 2.5])
    z = calc_bilinear_interpolation_on_grid(x_grid, y_grid, z_grid, x, y)
    assert z.shape == (4, 3)
    # a plane is reproduced exactly inside the grid; outside it the nearest
    # edge value is kept
    x_in, y_in = numpy.clip(x, 0., 20.), numpy.clip(y, -5., 5.)
    assert numpy.allclose(z, 2. * x_in[:, None] - 3. * y_in[None, :] + 1.)


def test_pd2d_background_as_in_refinement():
    background = Pd2dBackground(
        gamma=[2., 40., 80.], nu=[-5., 0., 10.],
        intensity=numpy.array([[1., 2., 3.], [4., 6., 5.], [9., 7., 8.]]))
    gamma = numpy.linspace(0., 90., 19)
    nu = numpy.linspace(-8., 12., 11)
    signal = background.interpolate_by_points(gamma, nu)
    signal_rhochi = calc_background(
        numpy.radians(gamma), numpy.radians(nu), numpy.radians(background.gamma),
        numpy.radians(background.nu), numpy.array(background.intensity))[0]
    assert numpy.allclose(signal, signal_rhochi)
