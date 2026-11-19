"""Magnetic domains (loop magnetic_domain) combined with structural twins.

A magnetic domain "R, theta" is the crystal transformed by the magnetic
operation: the same kernel as a twin of jana_twin with T = R^T and the time
reversal theta. Every structural twin has all the magnetic domains, so the
domain (s, m) has T_m T_s, theta_s theta_m and the fraction f_s f_m.
"""
import warnings

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.twin_domains import form_domains_by_dictionary
from cryspy.C_item_loop_classes.cl_1_magnetic_domain import parse_magnetic_operation
from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

from test_twin_field_domains import (
    INDEX_HKL, TWIN_LAW, calc_flip_ratios, form_crystal_ordered, form_experiment,
    matrix_ub_first_domain)

MOMENT = (0.5, 1.5, -1.0)
# direct-space operation of the twin law TWIN_LAW (h_2 = R^T h_1)
OPERATION_ROTATION = "y,-x,z"
# a structural twin: twofold axis along a
STRUCTURAL_LAW = numpy.diag([1., -1., -1.])


def form_magnetic_domain(operation_xyz: str, fraction: float) -> str:
    return "\n".join(["loop_", "_magnetic_domain_id", "_magnetic_domain_operation_xyz",
                      "_magnetic_domain_fraction", f"2 {operation_xyz:} {fraction:}"])


def form_jana_twin(l_twin) -> str:
    """l_twin: [(matrix T, fraction, time reversal), ...]"""
    ls_out = ["loop_", "_jana_twin_id"] + [
        f"_jana_twin_matrix[{i:}][{j:}]" for i in (1, 2, 3) for j in (1, 2, 3)] + [
        "_jana_twin_fraction", "_jana_twin_time_reversal"]
    for i_twin, (matrix, fraction, time_reversal) in enumerate(l_twin):
        ls_out.append(f"{i_twin + 2:} " + " ".join(f"{val:.1f}" for val in matrix.flatten())
                      + f" {fraction:} {time_reversal:}")
    return "\n".join(ls_out)


def test_parse_magnetic_operation():
    rotation, translation, time_reversal = parse_magnetic_operation("-y,x,z+1/2,-1")
    assert numpy.allclose(rotation, [0, -1, 0, 1, 0, 0, 0, 0, 1])
    assert numpy.allclose(translation, [0, 0, 0.5])
    assert time_reversal == -1
    assert parse_magnetic_operation("x,y,z")[2] == 1
    with pytest.raises(cryspy.CrysPyError):
        parse_magnetic_operation("x,y,z,0")


@pytest.mark.filterwarnings("ignore:The net ordered moment")
@pytest.mark.parametrize("time_reversal", (1, -1))
def test_magnetic_domain_is_magnetic_twin(time_reversal):
    """magnetic_domain 'R, theta' gives the flip ratios of jana_twin (R^T, theta)."""
    matrix_ub = matrix_ub_first_domain()
    crystal = form_crystal_ordered(MOMENT)
    text = form_experiment("exp", matrix_ub, INDEX_HKL)
    sign = "+1" if time_reversal > 0 else "-1"
    fr_domain = calc_flip_ratios(
        text + "\n" + form_magnetic_domain(f"{OPERATION_ROTATION},{sign}", 0.3), crystal)
    fr_twin = calc_flip_ratios(form_experiment(
        "exp", matrix_ub, INDEX_HKL, twin_fraction=0.3, time_reversal=time_reversal), crystal)
    assert numpy.allclose(fr_domain, fr_twin, rtol=1e-12, atol=1e-12)


@pytest.mark.filterwarnings("ignore:The net ordered moment")
def test_twins_times_magnetic_domains():
    """jana_twin (S, 0.2) with magnetic_domain (R, -1, 0.3) is the list of
    the three domains T_m T_s with the fractions f_s f_m."""
    matrix_ub = matrix_ub_first_domain()
    crystal = form_crystal_ordered(MOMENT)
    text = form_experiment("exp", matrix_ub, INDEX_HKL)
    fr_combined = calc_flip_ratios("\n".join([
        text, form_jana_twin([(STRUCTURAL_LAW, 0.2, 1)]),
        form_magnetic_domain(f"{OPERATION_ROTATION},-1", 0.3)]), crystal)
    fr_list = calc_flip_ratios("\n".join([text, form_jana_twin([
        (TWIN_LAW, 0.8 * 0.3, -1), (STRUCTURAL_LAW, 0.2 * 0.7, 1),
        (TWIN_LAW @ STRUCTURAL_LAW, 0.2 * 0.3, -1)])]), crystal)
    assert numpy.allclose(fr_combined, fr_list, rtol=1e-12, atol=1e-12)


def test_form_domains_by_dictionary():
    dict_diffrn = {
        "twin_fraction": numpy.array([0.2]),
        "twin_matrices": STRUCTURAL_LAW.reshape(9, 1),
        "twin_time_reversal": numpy.array([1.]),
        "magnetic_domain_fraction": numpy.array([0.3]),
        "magnetic_domain_matrices": TWIN_LAW.reshape(9, 1),
        "magnetic_domain_time_reversal": numpy.array([-1.])}
    d_domain = form_domains_by_dictionary(dict_diffrn)
    assert numpy.allclose(d_domain["fraction"], [0.56, 0.24, 0.14, 0.06])
    assert numpy.allclose(d_domain["time_reversal"], [1, -1, 1, -1])
    assert numpy.allclose(d_domain["matrices"][:, 3], (TWIN_LAW @ STRUCTURAL_LAW).flatten())
    assert form_domains_by_dictionary({}) is None


def test_derivatives_by_fractions():
    """Analytical derivatives of chi_sq by the twin and magnetic domain
    fractions against central finite differences."""
    matrix_ub = matrix_ub_first_domain()
    crystal = form_crystal_ordered(MOMENT)
    text = "\n".join([
        form_experiment("exp", matrix_ub, INDEX_HKL), form_jana_twin([(STRUCTURAL_LAW, 0.2, 1)]),
        form_magnetic_domain(f"{OPERATION_ROTATION},-1", 0.3)])
    global_object = cryspy.str_to_globaln(crystal + "\n" + text)
    global_object.diffrn_exp.jana_twin.items[0].fraction_refinement = True
    global_object.diffrn_exp.magnetic_domain.items[0].fraction_refinement = True
    d_global = global_object.get_dictionary()
    dict_crystal, dict_diffrn = d_global["crystal_fe"], d_global["diffrn_exp"]

    def calc(dict_diffrn_current):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return calc_chi_sq_for_diffrn_by_dictionary(
                dict_diffrn_current, dict_crystal, {}, flag_calc_analytical_derivatives=True)

    chi_sq, n_point, der_chi_sq, dder_chi_sq, l_parameter_name = calc(dict_diffrn)
    l_name = [name[1] for name in l_parameter_name]
    assert "twin_fraction" in l_name and "magnetic_domain_fraction" in l_name
    step = 1e-6
    for key in ("twin_fraction", "magnetic_domain_fraction"):
        l_value = []
        for sign in (1., -1.):
            dict_shifted = dict(dict_diffrn)
            dict_shifted[key] = dict_diffrn[key] + sign * step
            l_value.append(calc(dict_shifted)[0])
        der_numerical = (l_value[0] - l_value[1]) / (2. * step)
        der_analytical = der_chi_sq[l_name.index(key)]
        assert der_analytical == pytest.approx(der_numerical, rel=1e-4)
    assert chi_sq > 0.


def test_define_magnetic_domains_k0_mnf2():
    """MnF2 (P4_2'/mnm', k = 0): the only other domain is the 180° one."""
    import os
    f_name = os.path.join(os.path.dirname(__file__), "..", "MEM_test", "mnf2_moments", "main.rcif")
    obj = cryspy.load_file(f_name)
    diffrn = obj.diffrn_exp1
    diffrn.items[:] = [item for item in diffrn.items
                       if type(item).__name__ != "MagneticDomainL"]
    report = diffrn.define_magnetic_domains(obj.crystal_mnf2, flag_refine=True)
    assert "2 magnetic domains" in report
    assert diffrn.magnetic_domain.operation_xyz == ["x,y,z,-1"]
    assert diffrn.magnetic_domain.fraction == [0.5]
    assert "Magnetic domains: 2" in obj.crystal_mnf2.report_magnetic_domains()
