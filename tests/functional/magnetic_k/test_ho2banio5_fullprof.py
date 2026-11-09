"""Ho2BaNiO5, k = (1/2, 0, 1/2): comparison with the FullProf example hobk.

FullProf (hobk.pcr, D1B, lambda = 2.524 Å; one phase holds the nuclear and
the magnetic structure, Jbt = 10, operations of the group of k with MSYM)
gives Chi2 = 29.6 with the moments
Ho (0.125, 0, 8.994) μB and Ni (0.617, 0, -1.302) μB along a, b, c.
CrysPy reaches the same agreement and close moments when it starts from a
wrong magnetic structure.
"""
import os

import cryspy

F_RCIF = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))),
    "examples", "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5",
    "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5.rcif")


def test_ho2banio5_moments_as_in_fullprof():
    obj = cryspy.load_file(F_RCIF)
    fourier = obj.crystal_ho2banio5.atom_site_moment_fourier
    fourier["Ho"].real_x, fourier["Ho"].real_z = 0.0, 5.0
    fourier["Ni"].real_x, fourier["Ni"].real_z = 0.3, -0.5

    res = cryspy.rhochi_rietveld_refinement(obj)

    n_points = 400
    assert res["chi_sq"] / n_points < 25.
    assert abs(fourier["Ho"].real_z - 8.994) < 0.15
    assert abs(fourier["Ho"].real_x - 0.125) < 0.05
    assert abs(fourier["Ni"].real_z - (-1.302)) < 0.05
    assert abs(fourier["Ni"].real_x - 0.617) < 0.3


def test_ho2banio5_without_moments_is_worse():
    obj = cryspy.load_file(F_RCIF)
    fourier = obj.crystal_ho2banio5.atom_site_moment_fourier
    chi_sq_magnetic = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    for item in fourier.items:
        item.real_x, item.real_z = 0.0, 0.0
    chi_sq_nuclear = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    assert chi_sq_nuclear > 50. * chi_sq_magnetic
