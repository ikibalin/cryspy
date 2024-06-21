import numpy
import scipy
import scipy.interpolate
import scipy.special

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t

from cryspy.A_functions_base.unit_cell import \
    calc_sthovl_by_unit_cell_parameters, calc_matrix_t

from cryspy.A_functions_base.structure_factor import \
    calc_f_nucl_by_dictionary, calc_sft_ccs_by_dictionary, \
    calc_index_hkl_multiplicity_in_range

from cryspy.A_functions_base.integrated_intensity_powder_diffraction import \
    calc_powder_iint_1d_para

from cryspy.A_functions_base.preferred_orientation import calc_preferred_orientation_pd

from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_lorentz_factor

from cryspy.A_functions_base.powder_diffraction_tof import \
    calc_spectrum, calc_time_for_epithermal_neutrons_by_d, calc_time_for_thermal_neutrons_by_d, \
    calc_d_by_time_for_thermal_neutrons, calc_d_min_max_by_time_thermal_neutrons, \
    calc_d_min_max_by_time_epithermal_neutrons, calc_peak_shape_function

from cryspy.A_functions_base.powder_diffraction_tof_zcode import \
    calc_profile_by_zcode_parameters

from .rhochi_diffrn import get_flags
from .rhochi_pd import calc_background

na = numpy.newaxis


def calc_background_by_cosines(x, background_coefficients, x_min: float = 0., x_max: float = None, flag_background_coefficients: bool = False):
    if x_max is None:
        x_max = numpy.max(x)
    x_rel = ((x-x_min)/x_max)*numpy.pi
    val = numpy.linspace(
        0, background_coefficients.shape[0]-1, background_coefficients.shape[0])
    res = numpy.cos(numpy.expand_dims(val, axis=1) *
                    numpy.expand_dims(x_rel, axis=0))
    y = numpy.sum(numpy.expand_dims(
        background_coefficients, axis=1) * res, axis=0)
    dder = {}
    if flag_background_coefficients:
        dder["background_coefficients"] = res
    return y, dder


def calc_spectrum_incident(time, coefficients, type: str = "Maxwell", flag_coefficients: bool = False):
    exp = numpy.exp
    time_sq = numpy.square(time)
    time_4 = numpy.square(time_sq)
    a0, a1, a2, a3, a4 = coefficients[0], coefficients[1], coefficients[2], coefficients[3], coefficients[4]
    a5, a6, a7, a8 = coefficients[5], coefficients[6], coefficients[7], coefficients[8]
    if type == "Empirical-Exponents":
        res = a0 + a1 * exp(-a2 * time) + \
            a3 * exp(-a4 * time_sq) + \
            a5 * exp(-a6 * time*time_sq) + \
            a7 * exp(-a8 * time_4)
    else:  # type == "Maxwell"
        res = a0 + a1 * exp(-a2 * time_sq)/(time*time_4) + \
            a3 * exp(-a4 * time_sq) + \
            a5 * exp(-a6 * time*time_sq) + \
            a7 * exp(-a8 * numpy.square(time_sq))

    dder = {}
    if flag_coefficients:
        dder["coefficients"] = None
    return res, dder


def calc_chi_sq_for_tof_by_dictionary(
        dict_tof, dict_crystals, dict_in_out: dict = None, flag_use_precalculated_data: bool = False,
        flag_calc_analytical_derivatives: bool = False):
    """Calculate chi_sq for diffrn experiment.
    """
    if dict_in_out is None:
        flag_dict = False
        flag_use_precalculated_data = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
    dict_tof_keys = dict_tof.keys()

    phase_name = [hh["name"].lower() for hh in dict_crystals]

    excluded_points = dict_tof["excluded_points"]
    time = dict_tof["time"]
    dict_in_out["time"] = time
    neutron_type = dict_tof["neutron_type"]
    if neutron_type == "thermal":
        zero = dict_tof["zero"]
        dtt1 = dict_tof["dtt1"]
        dtt2 = dict_tof["dtt2"]
        d = calc_d_by_time_for_thermal_neutrons(time, zero, dtt1, dtt2)
        d_min_max = calc_d_min_max_by_time_thermal_neutrons(
            time, zero, dtt1, dtt2)
        flag_time = (
            numpy.any(dict_tof["flags_zero"]) or
            numpy.any(dict_tof["flags_dtt1"]) or
            numpy.any(dict_tof["flags_dtt2"]))
    else:  # epithermal
        zero = dict_tof["zero"]
        dtt1 = dict_tof["dtt1"]
        zerot = dict_tof["zerot"]
        dtt1t = dict_tof["dtt1t"]
        dtt2t = dict_tof["dtt2t"]
        flag_time = (
            numpy.any(dict_tof["flags_zero"]) or
            numpy.any(dict_tof["flags_dtt1"]) or
            numpy.any(dict_tof["flags_zerot"]) or
            numpy.any(dict_tof["flags_dtt1t"]) or
            numpy.any(dict_tof["flags_dtt2t"]))
        d_min_max = calc_d_min_max_by_time_epithermal_neutrons(
            time, zero, dtt1, zerot, dtt1t, dtt2t)
        raise AttributeError("Epithermal neutrons are not introudiced")

    sthovl_min = 0.5/d_min_max[1]
    sthovl_max = 0.5/d_min_max[0]

    if flag_dict:
        dict_in_out["d"] = d
        dict_in_out["excluded_points"] = excluded_points

    ttheta_bank = dict_tof["ttheta_bank"]
    lorentz_factor = calc_lorentz_factor(ttheta_bank, flag_ttheta=False)[0]

    if "beam_polarization" in dict_tof_keys:
        beam_polarization = dict_tof["beam_polarization"]
        flipper_efficiency = dict_tof["flipper_efficiency"]
        magnetic_field = dict_tof["magnetic_field"]
        flags_beam_polarization = dict_tof["flags_beam_polarization"]
        flags_flipper_efficiency = dict_tof["flags_flipper_efficiency"]
    else:
        beam_polarization, flipper_efficiency, magnetic_field = 0., 0., 0.
        flags_beam_polarization, flags_flipper_efficiency = False, False

    dict_tof_keys = dict_tof.keys()
    if "background_coefficients" in dict_tof_keys:
        background_coefficients = dict_tof["background_coefficients"]
        flags_background_coefficients = dict_tof["flags_background_coefficients"]

        flag_background_coefficients = numpy.any(flags_background_coefficients)
        if (flag_use_precalculated_data and ("signal_background" in dict_in_out) and not(flag_background_coefficients)):
            signal_background = dict_in_out["signal_background"]
        else:
            signal_background, dder_s_bkgr = calc_background_by_cosines(time, background_coefficients,
                                                                        flag_background_coefficients=(flag_background_coefficients and flag_calc_analytical_derivatives))
            dict_in_out["signal_background"] = signal_background
    elif "background_time" in dict_tof_keys:
        background_time = dict_tof["background_time"]
        background_intensity = dict_tof["background_intensity"]
        flags_background_intensity = dict_tof["flags_background_intensity"]
        flag_background_intensity = numpy.any(flags_background_intensity)
        if (flag_use_precalculated_data and ("signal_background" in dict_in_out) and not(flag_background_intensity)):
            signal_background = dict_in_out["signal_background"]
        else:
            signal_background, dder_s_bkgr = calc_background(
                time,
                background_time,
                background_intensity,
                flag_background_intensity=(flag_background_intensity and flag_calc_analytical_derivatives))
            dict_in_out["signal_background"] = signal_background
    else:
        signal_background = numpy.zeros_like(time)
        dict_in_out["signal_background"] = signal_background

    if "spectrum_incident_type" in dict_tof_keys:
        spectrum_incident_type = dict_tof["spectrum_incident_type"]
        spectrum_incident_coefficients = dict_tof["spectrum_incident_coefficients"]
        flags_spectrum_incident_coefficients = dict_tof["flags_spectrum_incident_coefficients"]
        flag_spectrum_incident_coefficients = numpy.any(
            flags_spectrum_incident_coefficients)
        if (flag_use_precalculated_data and ("spectrum_incident" in dict_in_out) and not(flag_spectrum_incident_coefficients)):
            spectrum_incident = dict_in_out["spectrum_incident"]
        else:
            spectrum_incident, dder_si = calc_spectrum_incident(
                time, spectrum_incident_coefficients, type=spectrum_incident_type)
        dict_in_out["spectrum_incident"] = spectrum_incident
    else:
        spectrum_incident = numpy.ones_like(time)
        dict_in_out["spectrum_incident"] = spectrum_incident

    phase_name = [hh["name"].lower() for hh in dict_crystals]

    pd_phase_name = dict_tof["phase_name"]
    phase_scale = dict_tof["phase_scale"]

    phase_ig = dict_tof["phase_ig"]  # IG_phase
    flags_phase_scale = dict_tof["flags_phase_scale"]
    flags_phase_ig = dict_tof["flags_phase_ig"]  # IG_phase

    profile_peak_shape = dict_tof["profile_peak_shape"]
    if profile_peak_shape == "pseudo-Voigt":
        profile_alphas = dict_tof["profile_alphas"]
        profile_betas = dict_tof["profile_betas"]
        profile_sigmas = dict_tof["profile_sigmas"]
        profile_gammas = dict_tof["profile_gammas"]
        flag_profile_shape = (
            numpy.any(dict_tof["flags_profile_alphas"]) or
            numpy.any(dict_tof["flags_profile_betas"]) or
            numpy.any(dict_tof["flags_profile_sigmas"]) or
            numpy.any(dict_tof["flags_profile_gammas"]))
    elif profile_peak_shape == "Gauss":
        profile_alphas = dict_tof["profile_alphas"]
        profile_betas = dict_tof["profile_betas"]
        profile_sigmas = dict_tof["profile_sigmas"]
        flag_profile_shape = (
            numpy.any(dict_tof["flags_profile_alphas"]) or
            numpy.any(dict_tof["flags_profile_betas"]) or
            numpy.any(dict_tof["flags_profile_sigmas"]))
    elif profile_peak_shape == "type0m":
        profile_alphas = dict_tof["profile_alphas"]
        profile_betas = dict_tof["profile_betas"]
        profile_sigmas = dict_tof["profile_sigmas"]
        profile_gammas = dict_tof["profile_gammas"]
        profile_rs = dict_tof["profile_rs"]
        flag_profile_shape = (
            numpy.any(dict_tof["flags_profile_alphas"]) or
            numpy.any(dict_tof["flags_profile_betas"]) or
            numpy.any(dict_tof["flags_profile_sigmas"]) or
            numpy.any(dict_tof["flags_profile_gammas"]) or
            numpy.any(dict_tof["flags_profile_rs"]))

    if "texture_name" in dict_tof_keys:
        flag_texture = True
        texture_name = dict_tof["texture_name"]
        texture_g1 = dict_tof["texture_g1"]
        texture_g2 = dict_tof["texture_g2"]
        texture_axis = dict_tof["texture_axis"]
        flags_texture_g1 = dict_tof["flags_texture_g1"]
        flags_texture_g2 = dict_tof["flags_texture_g2"]
        flags_texture_axis = dict_tof["flags_texture_axis"]
    else:
        flag_texture = False

    total_signal_plus = numpy.zeros_like(time)
    total_signal_minus = numpy.zeros_like(time)
    for p_name, p_scale, p_ig, flags_p_scale, flags_p_ig in zip(
            pd_phase_name, phase_scale, phase_ig, flags_phase_scale, flags_phase_ig):
        p_name = p_name.lower()
        flag_phase_texture = False
        if flag_texture:
            ind_texture = numpy.argwhere(texture_name == p_name)
            if ind_texture.shape[0] != 0:
                texture_g1 = texture_g1[ind_texture[0]]
                texture_g2 = texture_g2[ind_texture[0]]
                texture_axis = texture_axis[:, ind_texture[0]]
                flag_phase_texture = True
                flags_texture_g1 = flags_texture_g1[ind_texture[0]]
                flags_texture_g2 = flags_texture_g2[ind_texture[0]]
                flags_texture_axis = flags_texture_axis[:, ind_texture[0]]

        ind_phase = phase_name.index(p_name)
        dict_crystal = dict_crystals[ind_phase]
        dict_in_out_keys = dict_in_out.keys()
        if f"dict_in_out_{p_name:}" in dict_in_out_keys:
            dict_in_out_phase = dict_in_out[f"dict_in_out_{p_name:}"]
        else:
            dict_in_out_phase = {}
            dict_in_out[f"dict_in_out_{p_name:}"] = dict_in_out_phase

        dict_in_out_phase_keys = dict_in_out_phase.keys()

        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        translation_elems = dict_crystal["translation_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        unit_cell_parameters = dict_crystal["unit_cell_parameters"]
        flags_unit_cell_parameters = dict_crystal["flags_unit_cell_parameters"]
        flag_unit_cell_parameters = numpy.any(flags_unit_cell_parameters)

        if flag_unit_cell_parameters:
            sc_uc = dict_crystal["sc_uc"]
            v_uc = dict_crystal["v_uc"]
            unit_cell_parameters = numpy.dot(
                sc_uc, unit_cell_parameters) + v_uc

        if (flag_use_precalculated_data and
                ("index_hkl" in dict_in_out_phase_keys) and
                ("multiplicity_hkl" in dict_in_out_phase_keys) and not(flag_unit_cell_parameters)):
            index_hkl = dict_in_out_phase["index_hkl"]
            multiplicity_hkl = dict_in_out_phase["multiplicity_hkl"]
        else:
            if flag_phase_texture:
                reduced_symm_elems_p1 = numpy.array([[0], [0], [0], [1], [1], [0], [0], [
                                                    0], [1], [0], [0], [0], [1]], dtype=int)
                translation_elems_p1 = numpy.array(
                    [[0], [0], [0], [1]], dtype=int)
                index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                    sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems_p1, translation_elems_p1, centrosymmetry)
            else:
                index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                    sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems, translation_elems, centrosymmetry)

            if (("index_hkl" in dict_in_out_phase_keys) and flag_use_precalculated_data):
                if index_hkl.shape != dict_in_out_phase["index_hkl"].shape:
                    flag_use_precalculated_data = False
                else:
                    flag_use_precalculated_data = numpy.all(
                        numpy.logical_and(dict_in_out_phase["index_hkl"], index_hkl))

            dict_in_out_phase["index_hkl"] = index_hkl
            dict_in_out_phase["multiplicity_hkl"] = multiplicity_hkl

        flag_sthovl_hkl = flag_unit_cell_parameters
        sthovl_hkl, dder_sthovl_hkl = calc_sthovl_by_unit_cell_parameters(index_hkl,
                                                                          unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

        d_hkl = 0.5/sthovl_hkl
        if neutron_type == "thermal":
            time_hkl = calc_time_for_thermal_neutrons_by_d(
                d_hkl, zero, dtt1, dtt2)
        else:  # neutron_type == "epithermal"
            time_hkl = calc_time_for_epithermal_neutrons_by_d(
                d_hkl, zero, dtt1, zerot, dtt1t, dtt2t)

        wavelength_hkl = 2.*d_hkl*numpy.sin(0.5*ttheta_bank)
        wavelength_4_hkl = numpy.power(wavelength_hkl, 4)

        flag_time_hkl = flag_sthovl_hkl
        dict_in_out_phase["time_hkl"] = time_hkl
        dict_in_out_phase["d_hkl"] = d_hkl

        f_nucl, dder_f_nucl = calc_f_nucl_by_dictionary(
            dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
        flag_f_nucl = len(dder_f_nucl.keys()) > 0

        sft_ccs, dder_sft_ccs = calc_sft_ccs_by_dictionary(
            dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
        flag_sft_ccs = len(dder_sft_ccs.keys()) > 0

        flag_matrix_t = flag_unit_cell_parameters
        matrix_t, dder_matrix_t = calc_matrix_t(
            index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

        flag_tensor_sigma = flag_sft_ccs or flag_unit_cell_parameters
        tensor_sigma, dder_tensor_sigma = calc_m1_m2_m1t(
            matrix_t, sft_ccs, flag_m1=flag_sft_ccs, flag_m2=flag_unit_cell_parameters)

        flag_iint_plus_minus = flag_f_nucl or flag_tensor_sigma or flags_beam_polarization or flags_flipper_efficiency

        if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minu" in dict_in_out_phase_keys) and
                flag_use_precalculated_data and not(flag_iint_plus_minus)):
            iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
        else:
            iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_1d_para(
                f_nucl, tensor_sigma, beam_polarization, flipper_efficiency, magnetic_field,
                flag_f_nucl=flag_f_nucl, flag_tensor_sigma=flag_tensor_sigma,
                flag_polarization=flags_beam_polarization, flag_flipper=flags_flipper_efficiency)
            dict_in_out_phase["iint_plus"] = iint_plus
            dict_in_out_phase["iint_minus"] = iint_minus

        if flag_phase_texture:
            flag_texture_g1 = numpy.any(flags_texture_g1)
            flag_texture_g2 = numpy.any(flags_texture_g2)
            flag_texture_axis = numpy.any(flags_texture_axis)
            flag_hh = numpy.any(
                [flag_texture_g1, flag_texture_g2, flag_texture_axis])
            if (flag_use_precalculated_data and
                    ("preferred_orientation" in dict_in_out_phase_keys) and
                    not(flag_hh)):
                preferred_orientation = dict_in_out_phase["preferred_orientation"]
            else:
                preferred_orientation, dder_po = calc_preferred_orientation_tof(
                    index_hkl, texture_g1, texture_g2, texture_axis, unit_cell_parameters,
                    flag_texture_g1=flag_texture_g1 and flag_calc_analytical_derivatives,
                    flag_texture_g2=flag_texture_g2 and flag_calc_analytical_derivatives,
                    flag_texture_axis=flag_texture_axis and flag_calc_analytical_derivatives)
                dict_in_out_phase["preferred_orientation"] = preferred_orientation

        flag_profile_tof = (
            flag_unit_cell_parameters or flag_time or flag_profile_shape)
        if (("profile_tof" in dict_in_out_phase_keys) and
                flag_use_precalculated_data and not(flag_profile_tof)):
            profile_tof = dict_in_out_phase["profile_tof"]
        else:
            if profile_peak_shape == "pseudo-Voigt":
                profile_tof = calc_peak_shape_function(
                    profile_alphas, profile_betas, profile_sigmas,
                    d, time, time_hkl,
                    gammas=profile_gammas, size_g=0., size_l=0.,
                    strain_g=0., strain_l=0., peak_shape=profile_peak_shape)
            elif profile_peak_shape == "Gauss":
                profile_tof = calc_peak_shape_function(
                    profile_alphas, profile_betas, profile_sigmas,
                    d, time, time_hkl,
                    gammas=None, size_g=0., size_l=0.,
                    strain_g=0., strain_l=0., peak_shape=profile_peak_shape)
            elif profile_peak_shape == "type0m":
                time_2d = numpy.expand_dims(time, axis=1)
                time_hkl_2d = numpy.expand_dims(time_hkl, axis=0)
                d_hkl_2d = numpy.expand_dims(d_hkl, axis=0)
                delta_t_2d = time_2d - time_hkl_2d
                profile_sigmas_sq = numpy.square(profile_sigmas)
                profile_tof = calc_profile_by_zcode_parameters(
                    delta_t_2d, d_hkl_2d,
                    profile_sigmas_sq[0], profile_sigmas_sq[1], profile_sigmas_sq[2],
                    profile_gammas[0], profile_gammas[1], profile_gammas[2],
                    profile_rs[0], profile_rs[1], profile_rs[2],
                    profile_alphas[0], profile_alphas[1],
                    profile_betas[0], profile_betas[1], profile_betas[2])

            dict_in_out_phase["profile_tof"] = profile_tof

        # flags_p_scale
        iint_m_plus = iint_plus * multiplicity_hkl
        iint_m_minus = iint_minus * multiplicity_hkl

        dict_in_out_phase["iint_plus_with_factors"] = 0.5 * \
            p_scale * iint_m_plus*lorentz_factor*wavelength_4_hkl
        dict_in_out_phase["iint_minus_with_factors"] = 0.5 * \
            p_scale * iint_m_minus*lorentz_factor*wavelength_4_hkl
        if flag_texture:
            # 0.5 to have the same meaning for the scale factor as in FullProf
            signal_plus = 0.5 * p_scale * lorentz_factor * \
                (profile_tof * (iint_m_plus * wavelength_4_hkl *
                 preferred_orientation)[na, :]).sum(axis=1)  # sum over hkl
            signal_minus = 0.5 * p_scale * lorentz_factor * \
                (profile_tof * (iint_m_minus * wavelength_4_hkl *
                 preferred_orientation)[na, :]).sum(axis=1)
            dict_in_out_phase["iint_plus_with_factors"] *= preferred_orientation
            dict_in_out_phase["iint_minus_with_factors"] *= preferred_orientation
        else:
            signal_plus = 0.5 * p_scale * lorentz_factor * \
                (profile_tof * (iint_m_plus * wavelength_4_hkl)
                 [na, :]).sum(axis=1)
            signal_minus = 0.5 * p_scale * lorentz_factor * \
                (profile_tof * (iint_m_minus *
                 wavelength_4_hkl)[na, :]).sum(axis=1)

        dict_in_out_phase["signal_plus"] = signal_plus
        dict_in_out_phase["signal_minus"] = signal_minus
        total_signal_plus += signal_plus
        total_signal_minus += signal_minus

    total_signal_plus *= spectrum_incident
    total_signal_minus *= spectrum_incident
    if flag_dict:
        dict_in_out["signal_plus"] = total_signal_plus
        dict_in_out["signal_minus"] = total_signal_minus

    if ("signal_exp_plus" in dict_tof_keys) and ("signal_exp_minus" in dict_tof_keys):
        signal_exp_plus = dict_tof["signal_exp_plus"]
        signal_exp_minus = dict_tof["signal_exp_minus"]
        if flag_dict:
            dict_in_out["signal_exp_plus"] = signal_exp_plus
            dict_in_out["signal_exp_minus"] = signal_exp_minus
        flag_chi_sq_sum, flag_chi_sq_difference = True, True

        if "flag_chi_sq_sum" in dict_tof_keys:
            flag_chi_sq_sum = dict_tof["flag_chi_sq_sum"]

        if "flag_chi_sq_difference" in dict_tof_keys:
            flag_chi_sq_difference = dict_tof["flag_chi_sq_difference"]

        if flag_chi_sq_sum:
            signal_exp = signal_exp_plus[0, :] + signal_exp_minus[0, :]
            signal_sigma = numpy.sqrt(numpy.square(
                signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))

    else:
        signal_exp = dict_tof["signal_exp"][0, :]
        signal_sigma = dict_tof["signal_exp"][1, :]
        if flag_dict:
            dict_in_out["signal_exp"] = dict_tof["signal_exp"]
        flag_chi_sq_sum = True
        flag_chi_sq_difference = False

    chi_sq = 0.
    n_point = 0
    if flag_chi_sq_sum:
        in_points = numpy.logical_not(excluded_points)
        total_signal_sum = total_signal_plus + total_signal_minus + signal_background
        chi_sq_sum = ((numpy.square(
            (signal_exp - total_signal_sum)/signal_sigma)*in_points)).sum(axis=0)
        chi_sq += chi_sq_sum
        n_point += numpy.sum(in_points)

    if flag_chi_sq_difference:
        signal_exp_diff = signal_exp_plus[0, :] - signal_exp_minus[0, :]
        signal_sigma_diff = numpy.sqrt(numpy.square(
            signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))
        total_signal_diff = total_signal_plus - total_signal_minus
        chi_sq_diff = (numpy.square(
            (signal_exp_diff - total_signal_diff)/signal_sigma_diff)).sum(axis=0)
        chi_sq += chi_sq_diff
        n_point += signal_exp_diff.shape[0]
    if numpy.isnan(chi_sq):
        chi_sq = 1e30

    flags_pd = get_flags(dict_tof)
    l_flags_crystal = [get_flags(dict_crystal)
                       for dict_crystal in dict_crystals]

    l_parameter_name = []
    for way, flags in flags_pd.items():
        pd_type_name = dict_tof["type_name"]
        ind_1d = numpy.atleast_1d(numpy.argwhere(flags))  # .flatten()
        parameter_name = [(pd_type_name, ) + way + (tuple(ind_1d[ind, :]), )
                          for ind in range(ind_1d.shape[0])]
        l_parameter_name.extend(parameter_name)

    for flags_crystal, dict_crystal in zip(l_flags_crystal, dict_crystals):
        for way, flags in flags_crystal.items():
            crystal_type_name = dict_crystal["type_name"]
            ind_1d = numpy.atleast_1d(numpy.argwhere(flags))  # .flatten()
            parameter_name = [(crystal_type_name, ) + way + (tuple(ind_1d[ind, :]), )
                              for ind in range(ind_1d.shape[0])]
            l_parameter_name.extend(parameter_name)

    # if flag_background_coefficients:
    #     pass

    der_chi_sq = numpy.zeros((len(l_parameter_name), ), dtype=float)
    dder_chi_sq = numpy.zeros(
        (len(l_parameter_name), len(l_parameter_name)), dtype=float)

    return chi_sq, n_point, der_chi_sq, dder_chi_sq, l_parameter_name
