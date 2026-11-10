"""Ho2BaNiO5, k = (1/2, 0, 1/2): comparison with the FullProf example hobk.

FullProf (hobk.pcr, D1B, lambda = 2.524 Å; one phase holds the nuclear and
the magnetic structure, Jbt = 10, operations of the group of k with MSYM)
gives Chi2 = 29.6 with the moments
Ho (0.125, 0, 8.994) μB and Ni (0.617, 0, -1.302) μB along a, b, c.
CrysPy reaches the same agreement and close moments when it starts from a
wrong magnetic structure.
"""
import os

import numpy

import cryspy

F_RCIF = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))),
    "examples", "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5",
    "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5.rcif")
K_VECTOR = numpy.array([0.5, 0.0, 0.5])


def calc_moment(crystal, label):
    """Moment (3,) of the atom in the first cell, Mcos cos(2 pi k.r) +
    Msin sin(2 pi k.r) along a, b, c (μB)."""
    site = crystal.atom_site[label]
    angle = 2. * numpy.pi * numpy.dot(
        K_VECTOR, [site.fract_x, site.fract_y, site.fract_z])
    rows = crystal.atom_site_moment_fourier.take_rows(label)
    return numpy.array([rows[axis].cos * numpy.cos(angle)
                        + rows[axis].sin * numpy.sin(angle) for axis in "xyz"])


def set_moment(crystal, label, moment):
    """Write a moment of the first cell (k equivalent to -k)."""
    site = crystal.atom_site[label]
    angle = 2. * numpy.pi * numpy.dot(
        K_VECTOR, [site.fract_x, site.fract_y, site.fract_z])
    rows = crystal.atom_site_moment_fourier.take_rows(label)
    for axis, value in zip("xyz", moment):
        rows[axis].cos = value * numpy.cos(angle)
        rows[axis].sin = value * numpy.sin(angle)


def test_ho2banio5_moments_as_in_fullprof():
    obj = cryspy.load_file(F_RCIF)
    crystal = obj.crystal_ho2banio5
    set_moment(crystal, "Ho", (0.0, 0.0, 5.0))
    set_moment(crystal, "Ni", (0.3, 0.0, -0.5))

    res = cryspy.rhochi_rietveld_refinement(obj)

    n_points = 400
    assert res["chi_sq"] / n_points < 25.
    moment_ho, moment_ni = calc_moment(crystal, "Ho"), calc_moment(crystal, "Ni")
    assert abs(moment_ho[2] - 8.994) < 0.15
    assert abs(moment_ho[0] - 0.125) < 0.05
    assert abs(moment_ni[2] - (-1.302)) < 0.05
    assert abs(moment_ni[0] - 0.617) < 0.3


def test_ho2banio5_without_moments_is_worse():
    obj = cryspy.load_file(F_RCIF)
    crystal = obj.crystal_ho2banio5
    chi_sq_magnetic = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    for item in crystal.atom_site_moment_fourier.items:
        item.cos, item.sin = 0.0, 0.0
    chi_sq_nuclear = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    assert chi_sq_nuclear > 50. * chi_sq_magnetic
