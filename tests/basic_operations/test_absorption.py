"""Absorption correction of single crystal data (exptl_absorpt, exptl_crystal)."""
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.absorption import (
    calc_absorption_coefficients_by_composition, calc_beam_directions_by_angles,
    calc_beam_directions_normal_beam, calc_euler_zyz_by_rotation_matrix,
    calc_rotation_matrix_busing_levy, calc_rotation_matrix_by_euler_zyz,
    calc_transmission_ellipsoid, calc_transmission_sphere_zero_angle,
    take_neutron_cross_sections)
from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

DIR_FUNCTIONAL = os.path.join(os.path.dirname(os.path.dirname(__file__)), "functional")
F_TBTE3 = os.path.join(DIR_FUNCTIONAL, "twin_tbte3", "main.rcif")
F_MNBITE = os.path.join(DIR_FUNCTIONAL, "mnbi4te7_single_lambdaover2", "main.rcif")


def calc_transmission_monte_carlo(mu, semi_axes, matrix_axes, e_in, e_out, n_point=200000):
    """Independent reference: random points in the ellipsoid, paths from the
    quadratic equation of the ellipsoid surface."""
    rng = numpy.random.default_rng(1)
    points = rng.normal(size=(n_point, 3))
    points *= (rng.random(n_point)**(1. / 3.) / numpy.linalg.norm(points, axis=1))[:, None]
    m_axes = numpy.asarray(matrix_axes)
    semi_axes = numpy.asarray(semi_axes)
    m_a = m_axes @ numpy.diag(1. / semi_axes**2) @ m_axes.T  # x^T A x = 1 on the surface
    x = points * semi_axes[None, :] @ m_axes.T

    def path(direction):
        aa = direction @ m_a @ direction
        bb = 2. * x @ m_a @ direction
        cc = numpy.einsum("ni,ij,nj->n", x, m_a, x) - 1.
        return (-bb + numpy.sqrt(bb**2 - 4. * aa * cc)) / (2. * aa)

    return numpy.exp(-mu * (path(-e_in) + path(e_out))).mean()


@pytest.mark.parametrize("mu_radius", [0.1, 0.5, 1., 2.])
def test_sphere_at_zero_angle_is_analytical(mu_radius):
    e_beam = numpy.array([[1.], [0.], [0.]])
    transmission = calc_transmission_ellipsoid(mu_radius, [1., 1., 1.], numpy.eye(3), e_beam, e_beam)
    assert transmission[0] == pytest.approx(float(calc_transmission_sphere_zero_angle(mu_radius)), rel=1e-3)


def test_sphere_back_scattering():
    # transmission of a sphere at 2theta = 180 deg against a Monte Carlo integral
    e_in = numpy.array([[1.], [0.], [0.]])
    transmission = calc_transmission_ellipsoid(1., [1., 1., 1.], numpy.eye(3), e_in, -e_in)
    reference = calc_transmission_monte_carlo(1., [1., 1., 1.], numpy.eye(3), e_in[:, 0], -e_in[:, 0])
    assert transmission[0] == pytest.approx(reference, rel=3e-3)
    # back scattering is absorbed less than forward scattering
    # (International Tables C, Table 6.3.3.1: A* decreases with theta)
    assert transmission[0] > float(calc_transmission_sphere_zero_angle(1.))


def test_ellipsoid_against_monte_carlo():
    matrix_axes = calc_rotation_matrix_by_euler_zyz([20., 35., -50.])
    semi_axes = [1.5, 0.6, 0.25]
    e_in = numpy.array([0.3, 0.9, -0.2])
    e_in /= numpy.linalg.norm(e_in)
    e_out = numpy.array([-0.5, 0.4, 0.7])
    e_out /= numpy.linalg.norm(e_out)
    transmission = calc_transmission_ellipsoid(0.8, semi_axes, matrix_axes, e_in[:, None], e_out[:, None])
    reference = calc_transmission_monte_carlo(0.8, semi_axes, matrix_axes, e_in, e_out)
    assert transmission[0] == pytest.approx(reference, rel=3e-3)


def test_euler_angles_round_trip():
    matrix = calc_rotation_matrix_by_euler_zyz([10., 40., 75.])
    assert numpy.allclose(calc_rotation_matrix_by_euler_zyz(calc_euler_zyz_by_rotation_matrix(matrix)), matrix)
    crystal = cryspy.ExptlCrystal(size_max=3., size_mid=1., size_min=0.5)
    crystal.set_axes_by_directions([0., 0., 1.], [1., 1., 0.])
    m_axes = crystal.get_axes_matrix()
    assert numpy.allclose(m_axes[:, 0], [0., 0., 1.], atol=1e-6)
    assert numpy.allclose(m_axes[:, 1], numpy.array([1., 1., 0.]) / numpy.sqrt(2.), atol=1e-6)


def test_mu_of_iron():
    # bcc Fe, a = 2.8665 A: mu = 1.20 1/cm at 1.798 A (NIST neutron
    # scattering lengths: sigma_abs 2.56 b, sigma_scat 11.62 b)
    mu_abs, mu_scat = calc_absorption_coefficients_by_composition(
        ["Fe"], [2.], 2.8665**3, 1.798)
    assert mu_abs + mu_scat == pytest.approx(0.1204, rel=5e-3)
    assert take_neutron_cross_sections("Fe3+") == pytest.approx((2.56, 11.62))
    assert take_neutron_cross_sections("Gd")[0] == pytest.approx(49700.)


def test_normal_beam_directions():
    u_matrix = cryspy.DiffrnOrientMatrix(ub_11=0.2, ub_12=0.05, ub_13=0.0, ub_21=-0.04, ub_22=0.21,
                                         ub_23=0.03, ub_31=0.01, ub_32=-0.02, ub_33=0.15)
    u_matrix.form_object()
    m_u = u_matrix.u
    wavelength = 1.2
    index_hkl = numpy.array([[1., 2., 1.], [0., 1., 3.], [2., -1., 1.]]).T
    q_lab = u_matrix.ub @ index_hkl
    q_ccs = m_u.T @ q_lab
    for side in ("left", "right"):
        e_in, e_out, flag = calc_beam_directions_normal_beam(q_ccs, m_u[2], wavelength, side)
        assert numpy.all(flag)
        assert numpy.allclose(e_in.T @ m_u[2], 0.)  # incident beam is horizontal
        assert numpy.allclose(numpy.linalg.norm(e_in + wavelength * q_ccs, axis=0), 1.)
        assert numpy.allclose(e_out, e_in + wavelength * q_ccs)
    # the same position of the detector as DiffrnOrientMatrix.calc_angle
    e_in, e_out, flag = calc_beam_directions_normal_beam(q_ccs, m_u[2], wavelength, "left")
    for i_hkl in range(index_hkl.shape[1]):
        gamma, nu, phi = u_matrix.calc_angle(*index_hkl[:, i_hkl], wavelength=wavelength)
        # rotate the crystal so that the incident beam is along x of the laboratory
        e_lab_in, e_lab_out = m_u @ e_in[:, i_hkl], m_u @ e_out[:, i_hkl]
        angle = numpy.arctan2(e_lab_in[1], e_lab_in[0])
        rot = numpy.array([[numpy.cos(angle), numpy.sin(angle), 0.],
                           [-numpy.sin(angle), numpy.cos(angle), 0.], [0., 0., 1.]])
        k_f = rot @ e_lab_out
        assert numpy.degrees(numpy.arcsin(k_f[2])) == pytest.approx(nu, abs=1e-6)
        assert numpy.degrees(numpy.arctan2(k_f[1], k_f[0])) % 360. == pytest.approx(gamma % 360., abs=1e-6)


def test_four_circle_directions():
    m_u = calc_rotation_matrix_by_euler_zyz([15., 25., 35.])
    wavelength = 1.0
    omega, chi, phi = numpy.array([10., -20.]), numpy.array([30., 50.]), numpy.array([45., 100.])
    m_r = calc_rotation_matrix_busing_levy(omega, chi, phi)
    # detector of Busing & Levy (k_i along y): k_f = (sin g, cos g, 0) at gamma = 40 deg
    k_f_bl = numpy.array([numpy.sin(numpy.radians(40.)), numpy.cos(numpy.radians(40.)), 0.])
    l_q, l_out = [], []
    for i in range(2):
        q_phi_bl = m_r[:, :, i].T @ (k_f_bl - numpy.array([0., 1., 0.])) / wavelength
        k_phi_bl = m_r[:, :, i].T @ k_f_bl
        to_ccsl = numpy.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
        l_q.append(m_u.T @ to_ccsl @ q_phi_bl)
        l_out.append(m_u.T @ to_ccsl @ k_phi_bl)
    q_ccs = numpy.array(l_q).T
    e_in, e_out, misfit = calc_beam_directions_by_angles(q_ccs, m_u.ravel(), wavelength, omega, chi, phi)
    assert numpy.allclose(misfit, 0., atol=1e-10)
    assert numpy.allclose(e_out, numpy.array(l_out).T)


def test_refinement_model_is_multiplied_by_transmission():
    obj = cryspy.load_file(F_TBTE3)
    diffrn = obj.diffrn_t1p8nuc
    cryspy.rhochi_no_refinement(obj)
    intensity_0 = numpy.array(diffrn.diffrn_refln.intensity_calc)
    diffrn.add_items([cryspy.ExptlCrystal(size_max=3., size_mid=1., size_min=0.5,
                                          ellipsoid_euler_beta=40.)])
    text = cryspy.apply_absorption_correction(obj)
    assert "T_min" in text
    cryspy.rhochi_no_refinement(obj)
    intensity_1 = numpy.array(diffrn.diffrn_refln.intensity_calc)
    flag = intensity_0 > 1e-3 * intensity_0.max()
    ratio = intensity_1[flag] / intensity_0[flag]
    t_min, t_max = diffrn.exptl_absorpt.correction_t_min, diffrn.exptl_absorpt.correction_t_max
    assert 0. < t_min < t_max < 1.
    # rounded intensity_calc of the rcif items
    assert ratio.min() >= t_min - 1e-3 and ratio.max() <= t_max + 1e-3
    # the rcif text is read back
    s_cif = diffrn.exptl_absorpt.to_cif() + "\n" + diffrn.exptl_crystal.to_cif()
    assert "_exptl_absorpt_correction_type gaussian" in s_cif
    assert "_exptl_crystal_size_min 0.5" in s_cif
    items = cryspy.str_to_items(s_cif)
    assert {type(item).__name__ for item in items} == {"ExptlAbsorpt", "ExptlCrystal"}


def test_flip_ratios_change_only_by_lambda_half():
    obj = cryspy.load_file(F_MNBITE)
    diffrn = [item for item in obj.items if isinstance(item, cryspy.Diffrn)][0]
    cryspy.rhochi_no_refinement(obj)
    fr_0 = numpy.array(diffrn.diffrn_refln.fr_calc)
    diffrn.add_items([cryspy.ExptlAbsorpt(correction_type="gaussian", coefficient_mu=1.,
                                          coefficient_mu_half=1.),
                      cryspy.ExptlCrystal(size_max=4., size_mid=2., size_min=0.3)])
    cryspy.rhochi_no_refinement(obj)
    # the same mu at lambda/2: the transmission cancels exactly
    assert numpy.allclose(numpy.array(diffrn.diffrn_refln.fr_calc), fr_0, atol=2e-5)
    diffrn.exptl_absorpt.coefficient_mu_half = 0.5
    cryspy.rhochi_no_refinement(obj)
    assert numpy.abs(numpy.array(diffrn.diffrn_refln.fr_calc) - fr_0).max() > 1e-3


@pytest.mark.parametrize("f_name", [F_TBTE3, F_MNBITE])
def test_derivatives_of_absorption_parameters(f_name):
    """Intensities with twins, flip ratios with the lambda/2 contamination."""
    obj = cryspy.load_file(f_name)
    diffrn = [item for item in obj.items if isinstance(item, cryspy.Diffrn)][0]
    crystal = [item for item in obj.items if isinstance(item, cryspy.Crystal)][0]
    diffrn.add_items([cryspy.ExptlAbsorpt(correction_type="gaussian", coefficient_mu=1.,
                                          coefficient_mu_half=0.4),
                      cryspy.ExptlCrystal(size_max=3., size_mid=1., size_min=0.5,
                                          ellipsoid_euler_beta=40.)])
    for name in obj.get_variable_names():
        obj.set_variable_by_name(name[:-1] + ((f"{name[-1][0]:}_refinement", name[-1][1]), ), False)
    diffrn.exptl_crystal.size_min_refinement = True
    diffrn.exptl_absorpt.coefficient_mu_refinement = True
    dict_global = obj.get_dictionary()
    dict_diffrn, dict_crystal = dict_global[diffrn.get_name()], dict_global[crystal.get_name()]
    chi_sq, n_hkl, der_chi_sq, dder_chi_sq, l_name = calc_chi_sq_for_diffrn_by_dictionary(
        dict_diffrn, dict_crystal, dict_in_out={}, flag_calc_analytical_derivatives=True)
    l_checked = []
    for i_name, name in enumerate(l_name):
        key, index = name[1], name[2][0]
        if not key.startswith("absorption_"):
            continue
        l_checked.append(key)
        value = numpy.array(dict_diffrn[key], dtype=float)
        step = 1e-4 * abs(value[index])
        l_chi = []
        for sign in (1., -1.):
            dict_shift = dict(dict_diffrn)
            value_shift = value.copy()
            value_shift[index] += sign * step
            dict_shift[key] = value_shift
            l_chi.append(calc_chi_sq_for_diffrn_by_dictionary(dict_shift, dict_crystal, dict_in_out={})[0])
        der_numerical = (l_chi[0] - l_chi[1]) / (2. * step)
        assert abs(der_numerical) > 1e-6
        assert der_chi_sq[i_name] == pytest.approx(der_numerical, rel=1e-3)
    assert sorted(l_checked) == ["absorption_mu", "absorption_size"]


def test_estimate_and_apply_by_assistant():
    obj = cryspy.load_file(F_TBTE3)
    result = cryspy.call_assistant_tool(obj, "estimate_absorption", {})
    assert result["missing"][0]["name"] == "crystal_size"
    result = cryspy.call_assistant_tool(obj, "estimate_absorption", {"crystal_size": 2.})
    assert "NOT needed" in result["text"]
    assert "mu = 0.059" in result["text"]
    result = cryspy.call_assistant_tool(obj, "apply_absorption_correction", {})
    assert result["missing"][0]["name"] == "sizes"
    result = cryspy.call_assistant_tool(obj, "apply_absorption_correction", {"sizes": [3., 1., 0.5]})
    assert "T_min" in result["text"]
    assert obj.diffrn_t1p8nuc.exptl_absorpt.correction_type == "gaussian"
    assert "absorption_correction" in cryspy.call_assistant_tool(obj, "list_recipes", {})["text"]
    assert "Absorption correction" in cryspy.get_help("is an absorption correction needed")


def test_estimate_says_needed_for_strong_absorber():
    obj = cryspy.load_file(F_TBTE3)
    diffrn = obj.diffrn_t1p8nuc
    diffrn.add_items([cryspy.ExptlAbsorpt(correction_type="none", coefficient_mu=2.),
                      cryspy.ExptlCrystal(size_max=3., size_mid=1., size_min=0.3)])
    report = cryspy.estimate_absorption_correction(obj)
    assert "NEEDED" in report
    assert "ellipsoid 3 x 1 x 0.3 mm" in report
    # nothing is changed by the estimate
    assert diffrn.exptl_absorpt.correction_type == "none"


def test_density_and_mu_half():
    obj = cryspy.load_file(F_TBTE3)
    text = cryspy.calc_absorption_coefficient(obj)
    assert "density 7.749 g/cm^3" in text  # 4 TbTe3 in 464.3 A^3
    item = obj.diffrn_t1p8nuc.exptl_absorpt
    # the true absorption halves at lambda/2, the scattering does not
    assert 0.5 * item.coefficient_mu < item.coefficient_mu_half < item.coefficient_mu
    # an element starting with D is not deuterium
    from cryspy.procedure_rhochi.rhochi_absorption import calc_absorption_coefficient_of_crystal
    crystal = obj.crystal_tbte3
    crystal.atom_site.items[0].type_symbol = "Dy"
    d_mu = calc_absorption_coefficient_of_crystal(crystal, 1.)
    assert d_mu["density"] > 7. and d_mu["resonant"] == ["Dy"]


def test_estimate_without_crystal_and_missing_inputs():
    obj = cryspy.load_file(F_TBTE3)
    diffrn = obj.diffrn_t1p8nuc
    diffrn.add_items([cryspy.ExptlAbsorpt(correction_type="sphere", coefficient_mu=1.)])
    with pytest.raises(ValueError, match="size_rad"):
        cryspy.rhochi_no_refinement(obj)
    obj.items.remove(obj.crystal_tbte3)
    assert "cannot be estimated" in cryspy.estimate_absorption_correction(obj, crystal_size=1.)


def test_refinement_of_the_shape():
    obj = cryspy.load_file(F_TBTE3)
    diffrn = obj.diffrn_t1p8nuc
    diffrn.add_items([cryspy.ExptlCrystal(size_max=3., size_mid=1., size_min=0.5)])
    cryspy.apply_absorption_correction(obj)
    chi_sq_0 = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    for name in obj.get_variable_names():
        obj.set_variable_by_name(name[:-1] + ((f"{name[-1][0]:}_refinement", name[-1][1]), ), False)
    diffrn.phase.scale_refinement = True
    diffrn.exptl_crystal.size_min_refinement = True
    result = cryspy.rhochi_rietveld_refinement(obj)
    assert result["chi_sq"] < chi_sq_0
    assert "diffrn_t1p8nuc.exptl_crystal.size_min" in result["parameter_path"]
    assert diffrn.exptl_crystal.size_min_sigma > 0.
    assert diffrn.exptl_crystal.size_min != 0.5
