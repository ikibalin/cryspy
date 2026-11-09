"""Refinement of a single-crystal data set with excluded reflections."""
import copy
import os

import numpy

import cryspy
from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")


def _calc(excluded, twin_fraction):
    global_dict = cryspy.load_file(F_NAME).get_dictionary()
    dict_crystal = global_dict["crystal_tbte3"]
    dict_diffrn = copy.deepcopy(global_dict["diffrn_t1p8nuc"])
    dict_diffrn["flip_ratio_excluded"] = excluded
    dict_diffrn["twin_fraction"] = numpy.array([twin_fraction], dtype=float)
    dict_diffrn["flags_twin_fraction"] = numpy.array([True], dtype=bool)
    chi_sq, n_hkl, der_chi_sq, dder_chi_sq, names = calc_chi_sq_for_diffrn_by_dictionary(
        dict_diffrn, dict_crystal, dict_in_out={})
    ind_twin = [i for i, name in enumerate(names) if "twin_fraction" in name]
    return chi_sq, n_hkl, der_chi_sq[ind_twin], dder_chi_sq


def test_excluded_reflections_with_refined_parameter():
    n_hkl_all = cryspy.load_file(F_NAME).get_dictionary()["diffrn_t1p8nuc"]["index_hkl"].shape[1]
    excluded = numpy.zeros(n_hkl_all, dtype=bool)
    excluded[::3] = True
    twin_fraction = 0.35

    chi_sq, n_hkl, der_chi_sq, dder_chi_sq = _calc(excluded, twin_fraction)

    assert n_hkl == n_hkl_all - excluded.sum()
    assert numpy.all(numpy.isfinite(dder_chi_sq))
    step = 1e-6
    chi_sq_plus = _calc(excluded, twin_fraction + step)[0]
    chi_sq_minus = _calc(excluded, twin_fraction - step)[0]
    der_numerical = (chi_sq_plus - chi_sq_minus) / (2. * step)
    assert numpy.isclose(der_chi_sq[0], der_numerical, rtol=1e-5)
