import numpy
import scipy.special

from cryspy.A_functions_base.powder_diffraction_tof import (
    calc_d_min_max_by_time_thermal_neutrons,
    calc_hpv_eta,
    calc_lorentz_factor,
    calc_peak_shape_function,
    tof_Jorgensen,
    tof_Jorgensen_VonDreele,
    tof_non_convoluted_pseudo_voigt,
)
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile


def test_calc_d_min_max_by_time_thermal_neutrons_supports_zero_dtt2():
    time = numpy.array([3000.0, 19000.0], dtype=float)
    zero = numpy.array([2.921], dtype=float)
    dtt1 = numpy.array([6167.247], dtype=float)
    dtt2 = numpy.array([0.0], dtype=float)

    d_min_max = calc_d_min_max_by_time_thermal_neutrons(time, zero, dtt1, dtt2)

    expected = numpy.array([
        (time.min() - zero[0]) / dtt1[0],
        (time.max() - zero[0]) / dtt1[0],
    ], dtype=float)
    numpy.testing.assert_allclose(d_min_max, expected)

    tof_parameters = TOFParameters(
        zero=float(zero[0]),
        dtt1=float(dtt1[0]),
        dtt2=float(dtt2[0]),
        ttheta_bank=145.0,
    )
    numpy.testing.assert_allclose(tof_parameters.calc_d_min_max(time), expected)


def test_tof_lorentz_factor_uses_bank_angle():
    ttheta_bank = numpy.deg2rad(numpy.array([60.0, 100.0], dtype=float))

    actual, dder = calc_lorentz_factor(ttheta_bank, flag_ttheta=True)

    numpy.testing.assert_allclose(actual, numpy.sin(ttheta_bank))
    numpy.testing.assert_allclose(dder["ttheta"], numpy.cos(ttheta_bank))


def _expected_lorentz_profile(alpha, beta, gamma, time, time_hkl):
    norm = 0.5 * alpha * beta / (alpha + beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d - time_hkl_2d

    z1_2d = alpha[:, numpy.newaxis] * delta_2d + (0.5j * alpha * gamma)[:, numpy.newaxis]
    z2_2d = -beta[:, numpy.newaxis] * delta_2d + (0.5j * beta * gamma)[:, numpy.newaxis]

    with numpy.errstate(over="ignore", invalid="ignore"):
        term_1 = numpy.exp(z1_2d) * scipy.special.exp1(z1_2d)
        term_2 = numpy.exp(z2_2d) * scipy.special.exp1(z2_2d)

    term_1[numpy.isnan(term_1)] = 0.0
    term_1[numpy.isinf(term_1)] = 0.0
    term_2[numpy.isnan(term_2)] = 0.0
    term_2[numpy.isinf(term_2)] = 0.0

    return -2.0 * norm[:, numpy.newaxis] * (
        numpy.imag(term_1) + numpy.imag(term_2)
    ) / numpy.pi


def test_tof_jorgensen_von_dreele_matches_expected_lorentz_term():
    alpha = numpy.array([0.4], dtype=float)
    beta = numpy.array([0.2], dtype=float)
    sigma = numpy.array([10.0], dtype=float)
    gamma = numpy.array([5.0], dtype=float)
    time = numpy.linspace(-50.0, 50.0, 11)
    time_hkl = numpy.array([0.0], dtype=float)

    actual = tof_Jorgensen_VonDreele(alpha, beta, sigma, gamma, time, time_hkl)

    h_g_fwhm = sigma * numpy.sqrt(8.0 * numpy.log(2.0))
    h_pv, eta = calc_hpv_eta(h_g_fwhm, gamma)
    sigma_c = h_pv / numpy.sqrt(8.0 * numpy.log(2.0))
    gaussian = tof_Jorgensen(alpha, beta, sigma_c, time, time_hkl)
    lorentz = _expected_lorentz_profile(alpha, beta, h_pv, time, time_hkl)
    expected = (1.0 - eta)[:, numpy.newaxis] * gaussian + eta[:, numpy.newaxis] * lorentz

    numpy.testing.assert_allclose(actual, expected)
    assert numpy.isfinite(actual).all()
    assert numpy.all(actual >= 0.0)


def test_non_convoluted_pseudo_voigt_with_zero_gamma_is_gaussian():
    time = numpy.linspace(-3., 3., 41)
    time_hkl = numpy.array([0.25])
    sigma = numpy.array([0.9])
    gamma = numpy.array([0.])

    actual = tof_non_convoluted_pseudo_voigt(sigma, gamma, time, time_hkl)
    delta_t = time[:, numpy.newaxis] - time_hkl[numpy.newaxis, :]
    expected = numpy.exp(-0.5*numpy.square(delta_t/sigma)) / (
        numpy.sqrt(2.*numpy.pi)*sigma)

    numpy.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_non_convoluted_pseudo_voigt_uses_sigma_gamma_without_alpha_beta():
    time = numpy.linspace(-3., 3., 41)
    time_hkl = numpy.array([0.25])
    d = numpy.array([1.])
    sigmas = numpy.array([0.81, 0., 0.])
    gammas = numpy.array([0.35, 0., 0.])

    actual = calc_peak_shape_function(
        None, None, sigmas, d, time, time_hkl,
        gammas=gammas, peak_shape="non-conv-pseudo-Voigt")
    expected = tof_non_convoluted_pseudo_voigt(
        numpy.array([0.9]), numpy.array([0.35]), time, time_hkl)

    numpy.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_tof_profile_dictionary_non_convoluted_pseudo_voigt_has_no_alpha_beta():
    profile = TOFProfile(
        peak_shape="non-conv-pseudo-Voigt",
        sigma0=0.81,
        sigma1=0.,
        sigma2=0.,
        gamma0=0.35,
        gamma1=0.,
        gamma2=0.,
    )

    profile_dict = profile.get_dictionary()

    assert "profile_sigmas" in profile_dict
    assert "profile_gammas" in profile_dict
    assert "profile_alphas" not in profile_dict
    assert "profile_betas" not in profile_dict


def test_rhochi_tof_lorentz_factor_matches_cfml_convention(monkeypatch):
    from cryspy.procedure_rhochi import rhochi_tof

    time = numpy.array([1.0, 2.0, 3.0], dtype=float)
    ttheta_bank = numpy.deg2rad(100.0)
    d_hkl = numpy.array([2.0], dtype=float)
    multiplicity = numpy.array([2.0], dtype=float)
    iint_plus = numpy.array([3.0], dtype=float)
    iint_minus = numpy.array([5.0], dtype=float)
    phase_scale = numpy.array([7.0], dtype=float)

    monkeypatch.setattr(
        rhochi_tof,
        "calc_index_hkl_multiplicity_in_range",
        lambda *args, **kwargs: (
            numpy.array([[1], [0], [0]], dtype=int),
            multiplicity,
        ),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_sthovl_by_unit_cell_parameters",
        lambda *args, **kwargs: (0.5 / d_hkl, {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_f_nucl_by_dictionary",
        lambda *args, **kwargs: (numpy.array([1.0], dtype=float), {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_sft_ccs_by_dictionary",
        lambda *args, **kwargs: (numpy.zeros((3, 3, 1), dtype=float), {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_matrix_t",
        lambda *args, **kwargs: (numpy.zeros((3, 3, 1), dtype=float), {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_m1_m2_m1t",
        lambda *args, **kwargs: (numpy.zeros((3, 3, 1), dtype=float), {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_powder_iint_1d_para",
        lambda *args, **kwargs: (iint_plus, iint_minus, {}, {}),
    )
    monkeypatch.setattr(
        rhochi_tof,
        "calc_peak_shape_function",
        lambda *args, **kwargs: numpy.ones((time.size, d_hkl.size), dtype=float),
    )

    dict_tof = {
        "excluded_points": numpy.zeros(time.shape, dtype=bool),
        "time": time,
        "neutron_type": "thermal",
        "zero": numpy.array([0.0], dtype=float),
        "dtt1": numpy.array([1.0], dtype=float),
        "dtt2": numpy.array([0.0], dtype=float),
        "flags_zero": numpy.array([False], dtype=bool),
        "flags_dtt1": numpy.array([False], dtype=bool),
        "flags_dtt2": numpy.array([False], dtype=bool),
        "ttheta_bank": ttheta_bank,
        "phase_name": numpy.array(["phase"]),
        "phase_scale": phase_scale,
        "phase_ig": numpy.array([0.0], dtype=float),
        "flags_phase_scale": numpy.array([False], dtype=bool),
        "flags_phase_ig": numpy.array([False], dtype=bool),
        "profile_peak_shape": "Gauss",
        "profile_alphas": numpy.array([0.0, 0.0], dtype=float),
        "profile_betas": numpy.array([0.0, 0.0], dtype=float),
        "profile_sigmas": numpy.array([1.0, 0.0, 0.0], dtype=float),
        "flags_profile_alphas": numpy.array([False, False], dtype=bool),
        "flags_profile_betas": numpy.array([False, False], dtype=bool),
        "flags_profile_sigmas": numpy.array([False, False, False], dtype=bool),
        "profile_size_g": numpy.array([0.0], dtype=float),
        "profile_strain_g": numpy.array([0.0], dtype=float),
        "profile_size_l": numpy.array([0.0], dtype=float),
        "profile_strain_l": numpy.array([0.0], dtype=float),
        "flags_profile_size_g": numpy.array([False], dtype=bool),
        "flags_profile_strain_g": numpy.array([False], dtype=bool),
        "flags_profile_size_l": numpy.array([False], dtype=bool),
        "flags_profile_strain_l": numpy.array([False], dtype=bool),
        "signal_exp": numpy.stack(
            [numpy.zeros(time.shape, dtype=float), numpy.ones(time.shape, dtype=float)],
            axis=0,
        ),
        "type_name": "tof",
    }
    dict_crystal = {
        "name": "phase",
        "type_name": "crystal",
        "reduced_symm_elems": numpy.zeros((13, 1), dtype=int),
        "translation_elems": numpy.zeros((4, 1), dtype=int),
        "centrosymmetry": False,
        "unit_cell_parameters": numpy.ones(6, dtype=float),
        "flags_unit_cell_parameters": numpy.zeros(6, dtype=bool),
    }
    dict_in_out = {}

    rhochi_tof.calc_chi_sq_for_tof_by_dictionary(
        dict_tof,
        [dict_crystal],
        dict_in_out=dict_in_out,
    )

    expected_plus = (
        0.5 * phase_scale * iint_plus * multiplicity *
        numpy.power(d_hkl, 4) * numpy.sin(ttheta_bank)
    )
    expected_minus = (
        0.5 * phase_scale * iint_minus * multiplicity *
        numpy.power(d_hkl, 4) * numpy.sin(ttheta_bank)
    )
    dict_phase = dict_in_out["dict_in_out_phase"]

    numpy.testing.assert_allclose(
        dict_phase["iint_plus_with_factors"], expected_plus)
    numpy.testing.assert_allclose(
        dict_phase["iint_minus_with_factors"], expected_minus)
    numpy.testing.assert_allclose(
        dict_phase["signal_plus"], numpy.full(time.shape, expected_plus[0]))
    numpy.testing.assert_allclose(
        dict_phase["signal_minus"], numpy.full(time.shape, expected_minus[0]))


def test_tof_profile_dictionary_exports_size_strain_coefficients():
    profile = TOFProfile(
        peak_shape="pseudo-Voigt",
        sigma0=0.81,
        sigma1=0.,
        sigma2=0.,
        gamma0=0.35,
        gamma1=0.,
        gamma2=0.,
        alpha0=0.,
        alpha1=0.2971,
        beta0=0.04182,
        beta1=0.00224,
        size_g=0.25,
        strain_g=0.50,
        size_l=0.75,
        strain_l=1.00,
        size_g_refinement=True,
        strain_l_refinement=True,
    )

    profile_dict = profile.get_dictionary()

    numpy.testing.assert_allclose(profile_dict["profile_size_g"], [0.25])
    numpy.testing.assert_allclose(profile_dict["profile_strain_g"], [0.50])
    numpy.testing.assert_allclose(profile_dict["profile_size_l"], [0.75])
    numpy.testing.assert_allclose(profile_dict["profile_strain_l"], [1.00])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_size_g"], [True])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_strain_g"], [False])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_size_l"], [False])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_strain_l"], [True])


def test_tof_profile_dictionary_missing_size_strain_defaults_to_zero():
    profile = TOFProfile(
        peak_shape="pseudo-Voigt",
        sigma0=0.81,
        sigma1=0.,
        sigma2=0.,
        gamma0=0.35,
        gamma1=0.,
        gamma2=0.,
        alpha0=0.,
        alpha1=0.2971,
        beta0=0.04182,
        beta1=0.00224,
        size_g=None,
        strain_g=None,
    )

    profile_dict = profile.get_dictionary()

    numpy.testing.assert_allclose(profile_dict["profile_size_g"], [0.])
    numpy.testing.assert_allclose(profile_dict["profile_strain_g"], [0.])
    numpy.testing.assert_allclose(profile_dict["profile_size_l"], [0.])
    numpy.testing.assert_allclose(profile_dict["profile_strain_l"], [0.])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_size_g"], [False])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_strain_g"], [False])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_size_l"], [False])
    numpy.testing.assert_array_equal(profile_dict["flags_profile_strain_l"], [False])


def test_tof_size_strain_coefficients_match_sigma_gamma_paths():
    time = numpy.linspace(-3., 3., 41)
    time_hkl = numpy.array([0.25])
    d = numpy.array([1.2])
    alphas = numpy.array([0.4, 0.])
    betas = numpy.array([0.2, 0.])

    sigmas = numpy.array([0.81, 0.10, 0.20])
    size_g = numpy.array([0.05])
    strain_g = numpy.array([0.07])

    actual_gauss = calc_peak_shape_function(
        alphas, betas, sigmas, d, time, time_hkl,
        gammas=None, size_g=size_g, strain_g=strain_g,
        peak_shape="Gauss")
    expected_gauss = calc_peak_shape_function(
        alphas, betas,
        numpy.array([sigmas[0], sigmas[1] + strain_g[0], sigmas[2] + size_g[0]]),
        d, time, time_hkl, gammas=None, peak_shape="Gauss")

    gammas = numpy.array([0.35, 0.04, 0.03])
    size_l = numpy.array([0.06])
    strain_l = numpy.array([0.08])

    actual_pv = calc_peak_shape_function(
        alphas, betas, sigmas, d, time, time_hkl,
        gammas=gammas, size_g=size_g, strain_g=strain_g,
        size_l=size_l, strain_l=strain_l, peak_shape="pseudo-Voigt")
    expected_pv = calc_peak_shape_function(
        alphas, betas,
        numpy.array([sigmas[0], sigmas[1] + strain_g[0], sigmas[2] + size_g[0]]),
        d, time, time_hkl,
        gammas=numpy.array([gammas[0], gammas[1] + strain_l[0], gammas[2] + size_l[0]]),
        peak_shape="pseudo-Voigt")

    numpy.testing.assert_allclose(actual_gauss, expected_gauss)
    numpy.testing.assert_allclose(actual_pv, expected_pv)


def test_tof_cutoff_fwhm_truncates_far_points():
    alpha = numpy.array([0.4])
    beta = numpy.array([0.2])
    sigma = numpy.array([2.0])
    time = numpy.array([0.0, 100.0])
    time_hkl = numpy.array([0.0])

    full = tof_Jorgensen(alpha, beta, sigma, time, time_hkl)
    cut = tof_Jorgensen(alpha, beta, sigma, time, time_hkl, cutoff_fwhm=1.0)

    numpy.testing.assert_allclose(cut[0], full[0])
    assert full[1, 0] != 0.0
    assert cut[1, 0] == 0.0
    numpy.testing.assert_allclose(
        tof_Jorgensen(alpha, beta, sigma, time, time_hkl, cutoff_fwhm=0.0), full)
