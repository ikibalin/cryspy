import numpy
import scipy
import scipy.interpolate

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t, calc_m_v

from cryspy.A_functions_base.unit_cell import \
    calc_sthovl_by_unit_cell_parameters, calc_matrix_t, calc_eq_ccs_by_unit_cell_parameters

from cryspy.A_functions_base.structure_factor import \
    calc_f_nucl_by_dictionary, calc_sft_ccs_by_dictionary, \
    calc_index_hkl_multiplicity_in_range, calc_f_m_perp_ordered_by_dictionary

from cryspy.A_functions_base.integrated_intensity_powder_diffraction import \
    calc_powder_iint_2d_para, calc_powder_iint_2d_ordered, calc_powder_iint_2d_mix

from cryspy.A_functions_base.preferred_orientation import \
    calc_preferred_orientation_pd2d, calc_gamma_nu_for_textured_peaks

from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_profile_pseudo_voight_2d, calc_lorentz_factor, calc_ttheta_phi_by_gamma_nu

from .rhochi_diffrn import get_flags

na = numpy.newaxis


def calc_background(gamma, nu, background_gamma, background_nu, background_intensity, flag_background_intensity: bool = False):
    # f = scipy.interpolate.interp2d( 
    #     background_gamma, background_nu, background_intensity.transpose(), kind="linear", fill_value=None) #FIXME: not sure about transpose
    # intensity = f(gamma, nu).transpose()
    # dder = {}
    # if flag_background_intensity:
    #     dder["background_intensity"] = None

    ga_p = background_gamma
    nu_p = background_nu
    signal_p = background_intensity
    ga = gamma

    if ga.min() < ga_p.min():
        ga_p = numpy.insert(ga_p, 0, ga.min(), axis=0)
        signal_p = numpy.insert(signal_p, 0, signal_p[0,:], axis=0)

    if ga.max() >= ga_p.max():
        ga_p = numpy.append(ga_p, ga.max()+0.001)
        signal_p = numpy.append(signal_p, signal_p[-1:,:], axis=0)

    if nu.min() < nu_p.min():
        nu_p = numpy.insert(nu_p, 0, nu.min(), axis=0)
        signal_p = numpy.insert(signal_p, 0, signal_p[:,0], axis=1)

    if nu.max() >= nu_p.max():
        nu_p = numpy.append(nu_p, nu.max()+0.001)
        signal_p = numpy.append(signal_p, signal_p[:,-1:], axis=1)

    ga_left = ga_p[:-1]
    ga_right = ga_p[1:]
    nu_left = nu_p[:-1]
    nu_right = nu_p[1:]

    flags_ga = numpy.logical_and(
        ga[:, na] >= ga_left[na, :], 
        ga[:, na] < ga_right[na, :])

    flags_nu = numpy.logical_and(
        nu[:, na] >= nu_left[na, :],
        nu[:, na] < nu_right[na, :])

    arg_ga_1 = numpy.argwhere(flags_ga)[:,1]
    arg_ga_2 = arg_ga_1 + 1
    arg_nu_1 = numpy.argwhere(flags_nu)[:,1]
    arg_nu_2 = arg_nu_1 + 1

    f_11 = signal_p[arg_ga_1[:, na],arg_nu_1[na, :]]
    f_12 = signal_p[arg_ga_1[:, na],arg_nu_2[na, :]]
    f_21 = signal_p[arg_ga_2[:, na],arg_nu_1[na, :]]
    f_22 = signal_p[arg_ga_2[:, na],arg_nu_2[na, :]]

    c_11 = (ga_p[arg_ga_2] - ga)[:, na] * (nu_p[arg_nu_2] - nu)[na, :]
    c_12 = (ga_p[arg_ga_2] - ga)[:, na] * (nu - nu_p[arg_nu_1])[na, :]
    c_21 = (ga - ga_p[arg_ga_1])[:, na] * (nu_p[arg_nu_2] - nu)[na, :]
    c_22 = (ga - ga_p[arg_ga_1])[:, na] * (nu - nu_p[arg_nu_1])[na, :]
    denom = (ga_p[arg_ga_2] - ga_p[arg_ga_1])[:, na] * (nu_p[arg_nu_2] - nu_p[arg_nu_1])[na, :]

    signal = (c_11 * f_11 + c_12 * f_12 + c_21 * f_21 + c_22 * f_22)/denom
    dder = {}
    if flag_background_intensity:
        dder_b = numpy.zeros(ga.shape + nu.shape + ga_p.shape + nu_p.shape, dtype=float)
        dder_b[:, :, arg_ga_1[:, na],arg_nu_1[na, :]] += (c_11/denom)[na, na, :, :]
        dder_b[:, :, arg_ga_1[:, na],arg_nu_2[na, :]] += (c_12/denom)[na, na, :, :]
        dder_b[:, :, arg_ga_2[:, na],arg_nu_1[na, :]] += (c_21/denom)[na, na, :, :]
        dder_b[:, :, arg_ga_2[:, na],arg_nu_2[na, :]] += (c_22/denom)[na, na, :, :]
        dder["background_intensity"] = dder_b
    return signal, dder 


def calc_chi_sq_for_pd2d_by_dictionary(
        dict_pd, dict_crystals, dict_in_out: dict = None, flag_use_precalculated_data: bool=False,
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
    dict_pd_keys = dict_pd.keys()

    phase_name = [hh["name"].lower() for hh in dict_crystals]

    excluded_points = dict_pd["excluded_points"]
    dict_in_out["excluded_points"] = excluded_points

    gamma = dict_pd["gamma"]
    nu = dict_pd["nu"]
    offset_gamma = dict_pd["offset_gamma"]
    offset_nu = dict_pd["offset_nu"]
    flags_offset_gamma = dict_pd["flags_offset_gamma"]
    flags_offset_nu = dict_pd["flags_offset_nu"]

    gamma_zs = gamma - offset_gamma
    nu_zs = nu - offset_nu

    dict_in_out["gamma"] = gamma_zs
    dict_in_out["nu"] = nu_zs

    wavelength = dict_pd["wavelength"]
    flags_wavelength = dict_pd["flags_wavelength"]

    magnetic_field = dict_pd["magnetic_field"]

    flag_polarized, flag_unpolarized = False, False
    if "signal_exp_plus" in dict_pd_keys:
        flag_polarized = True
    elif "signal_exp" in dict_pd_keys:
        flag_unpolarized = True

    if flag_polarized:
        beam_polarization = dict_pd["beam_polarization"]
        flipper_efficiency = dict_pd["flipper_efficiency"]

        flags_beam_polarization = dict_pd["flags_beam_polarization"]
        flags_flipper_efficiency = dict_pd["flags_flipper_efficiency"]
    elif flag_unpolarized:
        beam_polarization = numpy.array([0.,], dtype=float)
        flipper_efficiency = numpy.array([0.,], dtype=float)
        flags_beam_polarization = numpy.array([False,], dtype=bool)
        flags_flipper_efficiency = numpy.array([False,], dtype=bool)

    alpha_det = numpy.arccos(
        numpy.sin(nu_zs)[na, :]/numpy.sqrt(
            2-2*numpy.cos(gamma_zs)[:, na] * numpy.cos(nu_zs)[na, :]
        )
    ) 
    dict_in_out["alpha_detector"] = alpha_det

    flag_ttheta = flags_offset_gamma or flags_offset_nu
    ttheta_zs, phi_zs, dder_ttheta_zs, dder_phi_zs = calc_ttheta_phi_by_gamma_nu(
        gamma_zs[:, na], nu_zs[na, :], flag_gamma=flags_offset_gamma, flag_nu=flags_offset_nu)

    ttheta_min = calc_ttheta_phi_by_gamma_nu(gamma_zs.min(), 0., flag_gamma=False, flag_nu=False)[0]
    ttheta_max = calc_ttheta_phi_by_gamma_nu(gamma_zs.max(), nu_zs.max(), flag_gamma=False, flag_nu=False)[0]

    sthovl_min = numpy.sin(0.5*ttheta_min - numpy.pi/90.)/wavelength
    if sthovl_min < 0:
        sthovl_min = 0.00001
    sthovl_max = numpy.sin(0.5*ttheta_max + numpy.pi/90.)/wavelength
    
    background_gamma = dict_pd["background_gamma"]
    background_nu = dict_pd["background_nu"]
    background_intensity = dict_pd["background_intensity"]
    flags_background_intensity = dict_pd["flags_background_intensity"]

    flag_background_intensity = numpy.any(flags_background_intensity)
    if (flag_use_precalculated_data and ("signal_background" in dict_in_out) and not(flag_background_intensity)):
        signal_background = dict_in_out["signal_background"]
    else:
        signal_background, dder_s_bkgr = calc_background(gamma, nu, background_gamma, background_nu, background_intensity,
            flag_background_intensity= (flag_background_intensity and flag_calc_analytical_derivatives))
        dict_in_out["signal_background"] = signal_background

    pd_phase_name = dict_pd["phase_name"]
    pd_phase_scale = dict_pd["phase_scale"]
    pd_phase_resolution_parameters = dict_pd["phase_resolution_parameters"] # U_phase, V_phase, W_phase, X_phase, Y_phase
    pd_phase_ig = dict_pd["phase_ig"] # IG_phase
    flags_pd_phase_scale = dict_pd["flags_phase_scale"]
    flags_pd_phase_resolution_parameters = dict_pd["flags_phase_resolution_parameters"] # U_phase, V_phase, W_phase, X_phase, Y_phase
    flags_pd_phase_ig = dict_pd["flags_phase_ig"] # IG_phase

    resolution_parameters = dict_pd["resolution_parameters"] # U, V, W, X, Y
    asymmetry_parameters = dict_pd["asymmetry_parameters"] # p1, p2, p3, p4

    p_phi = dict_pd["resolution_phi_parameter"] 
    flag_p_phi = dict_pd["flags_resolution_phi_parameter"] 

    flags_resolution_parameters = dict_pd["flags_resolution_parameters"] 
    flags_asymmetry_parameters = dict_pd["flags_asymmetry_parameters"] 
    flag_asymmetry_parameters = numpy.any(flags_asymmetry_parameters)
    

    if "texture_name" in dict_pd_keys:
        flag_texture = True
        pd_texture_name = dict_pd["texture_name"]
        pd_texture_g1 = dict_pd["texture_g1"]
        pd_texture_g2 = dict_pd["texture_g2"]
        pd_texture_axis = dict_pd["texture_axis"]
        pd_flags_texture_g1 = dict_pd["flags_texture_g1"]
        pd_flags_texture_g2 = dict_pd["flags_texture_g2"]
        pd_flags_texture_axis = dict_pd["flags_texture_axis"]
    else:
        flag_texture = False

    lorentz_factor, dder_lf = calc_lorentz_factor(ttheta_zs, flag_ttheta=flag_ttheta)
    dict_in_out["lorentz_factor"] = lorentz_factor


    total_signal_plus = numpy.zeros_like(ttheta_zs)
    total_signal_minus = numpy.zeros_like(ttheta_zs)
    for p_name, p_scale, p_resolution, p_ig, flags_p_scale, flags_p_resolution, flags_p_ig in zip(pd_phase_name, 
            pd_phase_scale, pd_phase_resolution_parameters.transpose(), pd_phase_ig,
            flags_pd_phase_scale, flags_pd_phase_resolution_parameters.transpose(), flags_pd_phase_ig):
        p_name=p_name.lower()
        flag_phase_texture = False
        if flag_texture:
            ind_texture = numpy.argwhere(pd_texture_name==p_name)
            if ind_texture.shape[0] != 0:
                texture_g1 = pd_texture_g1[ind_texture[0]]
                texture_g2 = pd_texture_g2[ind_texture[0]]
                texture_axis = pd_texture_axis[:, ind_texture[0]]
                flag_phase_texture = True
                flags_texture_g1 = pd_flags_texture_g1[ind_texture[0]]
                flags_texture_g2 = pd_flags_texture_g2[ind_texture[0]]
                flags_texture_axis = pd_flags_texture_axis[:, ind_texture[0]]
        
        ind_phase = phase_name.index(p_name)
        dict_crystal = dict_crystals[ind_phase]
        dict_crystal_keys = dict_crystal.keys()
        dict_in_out_keys = dict_in_out.keys()
        if f"dict_in_out_{p_name:}" in dict_in_out_keys:
            dict_in_out_phase = dict_in_out[f"dict_in_out_{p_name:}"]
        else:
            dict_in_out_phase = {}
            dict_in_out[f"dict_in_out_{p_name:}"] = dict_in_out_phase

        dict_in_out_phase_keys = dict_in_out_phase.keys()

        if "reduced_symm_elems" in dict_crystal_keys:
            reduced_symm_elems  = dict_crystal["reduced_symm_elems"]
            translation_elems = dict_crystal["translation_elems"]
        elif "full_symm_elems" in dict_crystal_keys:
            full_symm_elems = dict_crystal["full_symm_elems"]
        elif "full_mcif_elems" in dict_crystal_keys:
            full_mcif_elems  = dict_crystal["full_mcif_elems"]

        unit_cell_parameters = dict_crystal["unit_cell_parameters"]
        flags_unit_cell_parameters = dict_crystal["flags_unit_cell_parameters"]
        flag_unit_cell_parameters = numpy.any(flags_unit_cell_parameters)

        if flag_unit_cell_parameters:
            sc_uc = dict_crystal["sc_uc"]
            v_uc = dict_crystal["v_uc"]
            unit_cell_parameters = numpy.dot(sc_uc, unit_cell_parameters) + v_uc            

        if (flag_use_precalculated_data and 
                ("index_hkl" in dict_in_out_phase_keys) and 
                ("multiplicity_hkl" in dict_in_out_phase_keys) and not(flag_unit_cell_parameters or flag_ttheta)):
            index_hkl = dict_in_out_phase["index_hkl"]
            multiplicity_hkl = dict_in_out_phase["multiplicity_hkl"]
        else:
            if flag_phase_texture:
                reduced_symm_elems_p1 = numpy.array([[0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1]], dtype=int)
                translation_elems_p1 = numpy.array([[0], [0], [0], [1]], dtype=int)
                index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                    sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems_p1, translation_elems_p1, False)
                # index_hkl = numpy.reshape(numpy.stack([index_hkl, -1*index_hkl], axis=2), (3, 2*index_hkl.shape[1]))
                # multiplicity_hkl = numpy.ones_like(index_hkl[0])
            else:
                if "reduced_symm_elems" in dict_crystal_keys:
                    centrosymmetry = dict_crystal["centrosymmetry"]
                    index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                        sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems, translation_elems, centrosymmetry)
                else:
                    translation_elems_p1 = numpy.array([[0], [0], [0], [1]], dtype=int)
                    index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                        sthovl_min, sthovl_max, unit_cell_parameters, full_mcif_elems[:13], translation_elems_p1, False)

            if (("index_hkl" in dict_in_out_phase_keys) and flag_use_precalculated_data):
                if index_hkl.shape != dict_in_out_phase["index_hkl"].shape:
                    flag_use_precalculated_data = False
                else: 
                    flag_use_precalculated_data = numpy.all(numpy.logical_and(dict_in_out_phase["index_hkl"], index_hkl))

            dict_in_out_phase["index_hkl"] = index_hkl
            dict_in_out_phase["multiplicity_hkl"] = multiplicity_hkl

        flag_sthovl_hkl = flag_unit_cell_parameters
        sthovl_hkl, dder_sthovl_hkl = calc_sthovl_by_unit_cell_parameters(index_hkl,
            unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

        flag_ttheta_hkl = flag_sthovl_hkl or flags_wavelength
        ttheta_hkl = 2*numpy.arcsin(wavelength * sthovl_hkl)
        dict_in_out_phase["ttheta_hkl"] = ttheta_hkl

        f_nucl, dder_f_nucl = calc_f_nucl_by_dictionary(
            dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
        flag_f_nucl = len(dder_f_nucl.keys()) > 0

        flag_para = False
        if "atom_para_index" in dict_crystal_keys:
            sft_ccs, dder_sft_ccs = calc_sft_ccs_by_dictionary(
                dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
            flag_sft_ccs  = len(dder_sft_ccs.keys()) > 0
            flag_matrix_t = flag_unit_cell_parameters
            matrix_t, dder_matrix_t = calc_matrix_t(
                index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

            flag_tensor_sigma = flag_sft_ccs or flag_unit_cell_parameters
            tensor_sigma, dder_tensor_sigma = calc_m1_m2_m1t(matrix_t, sft_ccs, flag_m1=flag_sft_ccs, flag_m2=flag_unit_cell_parameters)
            flag_para = True

        flag_ordered = False
        if "atom_ordered_index" in dict_crystal_keys:
            f_m_perp_o_ccs, dder_f_m_perp_o_ccs = calc_f_m_perp_ordered_by_dictionary(
                dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
            flag_f_m_perp_o =  len(dder_f_m_perp_o_ccs.keys()) > 0
            flag_ordered = True

            flag_matrix_t = flag_unit_cell_parameters
            matrix_t, dder_matrix_t = calc_matrix_t(
                index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
            f_m_perp_o, dder_f_m_perp_o = calc_m_v(matrix_t, f_m_perp_o_ccs, flag_m=flag_unit_cell_parameters, flag_v=flag_f_m_perp_o)
            

        if flag_para and not(flag_ordered):
            flag_iint_plus_minus = flag_f_nucl or flag_tensor_sigma or flags_beam_polarization or flags_flipper_efficiency
            if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minu" in dict_in_out_phase_keys) and
                    flag_use_precalculated_data and not(flag_iint_plus_minus)):
                iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
            else:
                iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_2d_para(
                    f_nucl, tensor_sigma, beam_polarization, flipper_efficiency, magnetic_field,
                    alpha_det, dict_in_out_phase, 
                    flag_f_nucl=flag_f_nucl and flag_calc_analytical_derivatives,
                    flag_tensor_sigma=flag_tensor_sigma and flag_calc_analytical_derivatives,
                    flag_polarization=flags_beam_polarization and flag_calc_analytical_derivatives,
                    flag_flipper=flags_flipper_efficiency and flag_calc_analytical_derivatives)

        elif not(flag_para) and flag_ordered:
            flag_iint_plus_minus = flag_f_nucl or flag_f_m_perp_o or flags_beam_polarization or flags_flipper_efficiency
            if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minus" in dict_in_out_phase_keys) and
                    flag_use_precalculated_data and not(flag_iint_plus_minus)):
                iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
            else:
                iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_2d_ordered(
                    f_nucl, f_m_perp_o, beam_polarization, flipper_efficiency,
                    alpha_det, dict_in_out_phase, 
                    flag_f_nucl=flag_f_nucl and flag_calc_analytical_derivatives,
                    flag_f_m_perp=flag_f_m_perp_o and flag_calc_analytical_derivatives,
                    flag_polarization=flags_beam_polarization and flag_calc_analytical_derivatives,
                    flag_flipper=flags_flipper_efficiency and flag_calc_analytical_derivatives)

        elif flag_para and flag_ordered:
            flag_iint_plus_minus = flag_f_nucl or flag_tensor_sigma or flag_f_m_perp_o or flags_beam_polarization or flags_flipper_efficiency
            if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minu" in dict_in_out_phase_keys) and
                    flag_use_precalculated_data and not(flag_iint_plus_minus)):
                iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
            else:
                iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_2d_mix(
                    f_nucl, tensor_sigma, f_m_perp_o, beam_polarization, flipper_efficiency, magnetic_field,
                    alpha_det, dict_in_out_phase, 
                    flag_f_nucl=flag_f_nucl and flag_calc_analytical_derivatives,
                    flag_tensor_sigma=flag_tensor_sigma and flag_calc_analytical_derivatives,
                    flag_f_m_perp_ordered=flag_f_m_perp_o and flag_calc_analytical_derivatives,
                    flag_polarization=flags_beam_polarization and flag_calc_analytical_derivatives,
                    flag_flipper=flags_flipper_efficiency and flag_calc_analytical_derivatives)
        else:
            iint_plus = numpy.square(numpy.abs(f_nucl))
            iint_minus = numpy.square(numpy.abs(f_nucl))

        dict_in_out_phase["iint_plus"] = iint_plus
        dict_in_out_phase["iint_minus"] = iint_minus

        if flag_phase_texture:
            flag_texture_g1 = numpy.any(flags_texture_g1)
            flag_texture_g2 = numpy.any(flags_texture_g2)
            flag_texture_axis = numpy.any(flags_texture_axis)
            flag_hh = numpy.any([flag_texture_g1, flag_texture_g2, flag_texture_axis])
            if (flag_use_precalculated_data and 
                    ("preferred_orientation" in dict_in_out_phase_keys) and
                    not(flag_hh)):
                preferred_orientation = dict_in_out_phase["preferred_orientation"]
            else:
                preferred_orientation, dder_po = calc_preferred_orientation_pd2d(alpha_det,
                    index_hkl, texture_g1, texture_g2, texture_axis, unit_cell_parameters,
                    flag_texture_g1=flag_texture_g1 and flag_calc_analytical_derivatives,
                    flag_texture_g2=flag_texture_g2 and flag_calc_analytical_derivatives,
                    flag_texture_axis=flag_texture_axis and flag_calc_analytical_derivatives)
                dict_in_out_phase["preferred_orientation"] = preferred_orientation

                # it is not necessary calculations but it is better to now (gamma,nu)_max of peaks

                eq_axis = calc_eq_ccs_by_unit_cell_parameters(
                    texture_axis, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
                if "eq_ccs" in dict_in_out_phase_keys:
                    eq_ccs = dict_in_out_phase["eq_ccs"]
                else:
                    eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters=unit_cell_parameters, flag_unit_cell_parameters=False)[0]
                    dict_in_out_phase["eq_ccs"] = eq_ccs
                ttheta_hkl = dict_in_out_phase["ttheta_hkl"]
                gamma_hkl, nu_hkl = calc_gamma_nu_for_textured_peaks(
                    eq_axis, eq_ccs, ttheta_hkl, texture_g1)
                dict_in_out_phase["gamma_hkl"] = gamma_hkl
                dict_in_out_phase["nu_hkl"] = nu_hkl
        
        flag_rp = numpy.any(flags_p_resolution) or numpy.any(flags_resolution_parameters)

        hh = resolution_parameters + p_resolution
        u, v, w, x, y = hh[0], hh[1], hh[2], hh[3], hh[4]
        p_1, p_2, p_3, p_4 = asymmetry_parameters[0], asymmetry_parameters[1], asymmetry_parameters[2], asymmetry_parameters[3]

        flag_profile_pv = flag_rp or flag_asymmetry_parameters or flag_p_phi
        if flag_use_precalculated_data and ("profile_pv" in dict_in_out_keys) and not(flag_profile_pv):
            profile_pv = dict_in_out_phase["profile_pv"]
        else:
            profile_pv, dder_pv = calc_profile_pseudo_voight_2d(ttheta_zs, phi_zs, ttheta_hkl, u, v, w, p_ig, x, y,
                p_1, p_2, p_3, p_4, 
                p_phi,
                flag_ttheta=flag_ttheta, flag_ttheta_hkl=flag_ttheta_hkl, flag_u=flag_rp,
                flag_v=flag_rp and flag_calc_analytical_derivatives,
                flag_w=flag_rp and flag_calc_analytical_derivatives,
                flag_i_g=flags_p_ig and flag_calc_analytical_derivatives,
                flag_x=flag_rp and flag_calc_analytical_derivatives,
                flag_y=flag_rp and flag_calc_analytical_derivatives, 
                flag_p_1=flag_asymmetry_parameters and flag_calc_analytical_derivatives,
                flag_p_2=flag_asymmetry_parameters and flag_calc_analytical_derivatives,
                flag_p_3=flag_asymmetry_parameters and flag_calc_analytical_derivatives,
                flag_p_4=flag_asymmetry_parameters and flag_calc_analytical_derivatives,
                flag_p_phi=flag_p_phi and flag_calc_analytical_derivatives)
            dict_in_out_phase["profile_pv"] = profile_pv


        # flags_p_scale
        iint_m_plus = iint_plus * multiplicity_hkl
        iint_m_minus = iint_minus * multiplicity_hkl


        if flag_phase_texture:
            # 0.5 to have the same meaning for the scale factor as in FullProf
            signal_plus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_plus * preferred_orientation).sum(axis=2) # sum over hkl
            signal_minus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_minus * preferred_orientation).sum(axis=2) 

        else:
            signal_plus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_plus).sum(axis=2) 
            signal_minus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_minus).sum(axis=2) 

        dict_in_out_phase["signal_plus"] = signal_plus
        dict_in_out_phase["signal_minus"] = signal_minus
        total_signal_plus += signal_plus
        total_signal_minus += signal_minus

    dict_in_out["signal_plus"] = total_signal_plus
    dict_in_out["signal_minus"] = total_signal_minus

    if flag_polarized:
        signal_exp_plus = dict_pd["signal_exp_plus"]
        signal_exp_minus = dict_pd["signal_exp_minus"]

        dict_in_out["signal_exp_plus"] = signal_exp_plus
        dict_in_out["signal_exp_minus"] = signal_exp_minus
        flag_chi_sq_sum = dict_pd["flag_chi_sq_sum"]
        flag_chi_sq_difference = dict_pd["flag_chi_sq_difference"]

    elif flag_unpolarized:
        signal_exp = dict_pd["signal_exp"]
        dict_in_out["signal_exp"] = signal_exp
        flag_chi_sq_sum = True
        flag_chi_sq_difference = False

    chi_sq = 0.
    n_point = 0
    if flag_chi_sq_sum:
        in_points = numpy.logical_not(excluded_points)
        if flag_polarized:
            signal_exp_sum = signal_exp_plus[0, :] + signal_exp_minus[0, :]
            signal_sigma_sum = numpy.sqrt(numpy.square(signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))
        elif flag_unpolarized:
            signal_exp_sum = signal_exp[0, :]
            signal_sigma_sum = signal_exp[1, :]
        nan_points_sum = numpy.logical_not(numpy.isnan(signal_exp_sum))
        in_points = numpy.logical_and(in_points, nan_points_sum)

        total_signal_sum = (total_signal_plus + total_signal_minus)
        diff_signal_sum = signal_exp_sum - total_signal_sum - signal_background
        inv_sigma_sq_sum = numpy.square(1./signal_sigma_sum)
        chi_sq_sum = numpy.sum(numpy.square(diff_signal_sum[in_points])*inv_sigma_sq_sum[in_points])
        chi_sq += chi_sq_sum
        n_point += numpy.sum(in_points)

    if flag_chi_sq_difference:
        signal_exp_diff = signal_exp_plus[0, :] - signal_exp_minus[0, :]
        signal_sigma_diff = numpy.sqrt(numpy.square(signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))
        total_signal_diff = total_signal_plus - total_signal_minus
        nan_points_diff = numpy.logical_not(numpy.isnan(signal_exp_diff))

        chi_sq_diff = numpy.sum(numpy.square((signal_exp_diff - total_signal_diff)[nan_points_diff]/signal_sigma_diff[nan_points_diff]))
        chi_sq += chi_sq_diff
        n_point += numpy.sum(nan_points_diff)

    if numpy.isnan(chi_sq):
        chi_sq = 1e30

    flags_pd2d = get_flags(dict_pd)
    l_flags_crystal = [get_flags(dict_crystal) for dict_crystal in dict_crystals]

    dder_signal_pd2d_sum = {}
    dder_signal_pd2d_difference = {}

    if flag_calc_analytical_derivatives:
        if flag_background_intensity:
            dder_signal_pd2d_sum['background_intensity'] = dder_s_bkgr['background_intensity']

    dder_signal_pd2d_sum_keys = dder_signal_pd2d_sum.keys()
    dder_signal_pd2d_difference_keys = dder_signal_pd2d_difference.keys()

    l_dder_plus_p = []
    l_parameter_name = []
    for way, flags in flags_pd2d.items():
        if flag_calc_analytical_derivatives:
            name = way[0]
            if name in dder_signal_pd2d_sum_keys:
                dder_plus_p = dder_signal_pd2d_sum[name][:, :, flags]
            l_dder_plus_p.append(dder_plus_p)

        pd_type_name = dict_pd["type_name"]
        ind_1d = numpy.atleast_1d(numpy.argwhere(flags)) #.flatten()
        parameter_name = [(pd_type_name, ) + way + (tuple(ind_1d[ind,:]), ) for ind in range(ind_1d.shape[0])]
        l_parameter_name.extend(parameter_name)

    for flags_crystal, dict_crystal in zip(l_flags_crystal, dict_crystals):
        for way, flags in flags_crystal.items():
            crystal_type_name = dict_crystal["type_name"]
            ind_1d = numpy.atleast_1d(numpy.argwhere(flags)) #.flatten()
            parameter_name = [(crystal_type_name, ) + way + (tuple(ind_1d[ind,:]), ) for ind in range(ind_1d.shape[0])]
            l_parameter_name.extend(parameter_name)


    der_chi_sq = numpy.zeros((len(l_parameter_name), ), dtype=float) 
    dder_chi_sq = numpy.zeros((len(l_parameter_name), len(l_parameter_name)), dtype=float)
    if flag_calc_analytical_derivatives:
        if len(l_dder_plus_p) > 0:
            dder_plus_p = numpy.concatenate(l_dder_plus_p, axis=2)
            if flag_chi_sq_sum:
                der_chi_sq = (((-2.*diff_signal_sum*inv_sigma_sq_sum)[:, :, na] * dder_plus_p).sum(axis=1)).sum(axis=0)
                dder_chi_sq = ((2*dder_plus_p[:, :, :, na] * dder_plus_p[:, :, na, :]).sum(axis=1)).sum(axis=0)

    return chi_sq, n_point, der_chi_sq, dder_chi_sq, l_parameter_name
    