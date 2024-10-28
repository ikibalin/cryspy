import numpy
import scipy
import scipy.interpolate

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t, calc_m_v

from cryspy.A_functions_base.unit_cell import \
    calc_sthovl_by_unit_cell_parameters, calc_matrix_t

from cryspy.A_functions_base.structure_factor import \
    calc_f_nucl_by_dictionary, \
    calc_f_charge_by_dictionary, \
    calc_sft_ccs_by_dictionary, \
    calc_index_hkl_multiplicity_in_range, \
    calc_f_m_perp_ordered_by_dictionary

from cryspy.A_functions_base.integrated_intensity_powder_diffraction import \
    calc_powder_iint_1d_para, calc_powder_iint_1d_ordered, calc_powder_iint_1d_mix

from cryspy.A_functions_base.preferred_orientation import calc_preferred_orientation_pd

from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_profile_pseudo_voight, calc_lorentz_factor

from .rhochi_diffrn import get_flags


na = numpy.newaxis

def calc_background(ttheta, background_ttheta, background_intensity, flag_background_intensity: bool = False):
    x_p = numpy.copy(background_ttheta)
    y_p = numpy.copy(background_intensity)
    x_min = ttheta.min()
    x_max = ttheta.max()
    if x_p.min() > x_min:
        y_0 = (y_p[1]-y_p[0])*(x_min - x_p[0])/(x_p[1]-x_p[0]) + y_p[0]
        x_p = numpy.insert(x_p, 0, x_min)
        y_p = numpy.insert(y_p, 0, y_0)
    if x_p.max() <= x_max:
        x_max = x_max + 1.
        y_last = (y_p[-1]-y_p[-2])*(x_max - x_p[-2])/(x_p[-1]-x_p[-2]) + y_p[-2]
        x_p = numpy.append(x_p, x_max)
        y_p = numpy.append(y_p, y_last)
    x_left = x_p[:-1]
    x_right = x_p[1:]
    flags = numpy.logical_and(ttheta[:, na] >= x_left[na, :], ttheta[:, na] < x_right[na, :])
    p0 = numpy.argwhere(flags)[:,1]
    p1 = p0 + 1 
    intensity = (y_p[p1]-y_p[p0]) * (ttheta-x_p[p0])/(x_p[p1]-x_p[p0]) + y_p[p0]
    # f = scipy.interpolate.interp1d(
    #     background_ttheta, background_intensity, kind="linear", fill_value="extrapolate")
    # intensity = f(ttheta)
    dder = {}
    if flag_background_intensity:
        ttheta_shift = ttheta[:, na] - background_ttheta[na, :]
        diff_b_tth = background_ttheta[1:]-background_ttheta[:-1]
        # y_n + (y_np1-y_np)*(x-x_n)/(x_np1-x_n)
        # 1  -(x-x_n)/(x_np1-x_n) + (x-x_nm1)/(x_n-x_nm1) 
        dder_bkgr = numpy.zeros((ttheta.size, background_ttheta.size), dtype=float)
        dder["background_intensity"] = dder_bkgr 
    return intensity, dder 


def calc_chi_sq_for_pd_by_dictionary(
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
    ttheta = dict_pd["ttheta"]
    offset_ttheta = dict_pd["offset_ttheta"]
    ttheta_zs = ttheta - offset_ttheta
    flags_offset_ttheta = dict_pd["flags_offset_ttheta"]
    if flag_dict:
        dict_in_out["ttheta"] = ttheta
        dict_in_out["ttheta_corrected"] = ttheta_zs  
        dict_in_out["excluded_points"] = excluded_points 


    wavelength = dict_pd["wavelength"]
    flags_wavelength = dict_pd["flags_wavelength"]
    radiation = dict_pd["radiation"]

    if "beam_polarization" in dict_pd_keys:
        beam_polarization = dict_pd["beam_polarization"]
        flipper_efficiency = dict_pd["flipper_efficiency"]
        magnetic_field = dict_pd["magnetic_field"]
        flags_beam_polarization = dict_pd["flags_beam_polarization"]
        flags_flipper_efficiency = dict_pd["flags_flipper_efficiency"]
    else:
        beam_polarization, flipper_efficiency, magnetic_field = 0., 0., 0.
        flags_beam_polarization, flags_flipper_efficiency = False, False

    sthovl_min = numpy.sin(0.5*ttheta_zs.min() - numpy.pi/90.)/wavelength
    if sthovl_min <= 0:
        sthovl_min = 0.0001
    sthovl_max = numpy.sin(0.5*ttheta_zs.max() + numpy.pi/90.)/wavelength
    if sthovl_max <= sthovl_min:
        sthovl_max = sthovl_min+0.01
        if sthovl_max >= 1.:
            sthovl_max =  0.99999/wavelength
        
    background_ttheta = dict_pd["background_ttheta"]
    background_intensity = dict_pd["background_intensity"]
    flags_background_intensity = dict_pd["flags_background_intensity"]

    flag_background_intensity = numpy.any(flags_background_intensity)
    if (flag_use_precalculated_data and ("signal_background" in dict_in_out) and not(flag_background_intensity)):
        signal_background = dict_in_out["signal_background"]
    else:
        signal_background, dder_s_bkgr = calc_background(ttheta, background_ttheta, background_intensity,
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
    if "asymmetry_parameters" in dict_pd_keys:
        asymmetry_parameters = dict_pd["asymmetry_parameters"] # p1, p2, p3, p4
        flags_asymmetry_parameters = dict_pd["flags_asymmetry_parameters"] 
        flag_asymmetry_parameters = numpy.any(flags_asymmetry_parameters)
        p_1, p_2, p_3, p_4 = asymmetry_parameters[0], asymmetry_parameters[1], asymmetry_parameters[2], asymmetry_parameters[3]

    else:
        p_1, p_2, p_3, p_4 = 0., 0., 0., 0.
        flag_asymmetry_parameters = False

    flags_resolution_parameters = dict_pd["flags_resolution_parameters"] 
    
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
    
    k = dict_pd["k"]
    cthm = dict_pd["cthm"]
    
    lorentz_factor, dder_lf = calc_lorentz_factor(ttheta_zs, k=k, cthm=cthm, flag_ttheta=flags_offset_ttheta)
    dict_in_out["lorentz_factor"] = lorentz_factor


    total_signal_plus = numpy.zeros_like(ttheta_zs)
    total_signal_minus = numpy.zeros_like(ttheta_zs)
    for p_name, p_scale, p_resolution, p_ig, flags_p_scale, flags_p_resolution, flags_p_ig in zip(pd_phase_name, 
            pd_phase_scale, pd_phase_resolution_parameters.transpose(), pd_phase_ig,
            flags_pd_phase_scale, flags_pd_phase_resolution_parameters.transpose(), flags_pd_phase_ig):
        p_name = p_name.lower()
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
        dict_in_out_keys = dict_in_out.keys()
        if f"dict_in_out_{p_name:}" in dict_in_out_keys:
            dict_in_out_phase = dict_in_out[f"dict_in_out_{p_name:}"]
        else:
            dict_in_out_phase = {}
            dict_in_out[f"dict_in_out_{p_name:}"] = dict_in_out_phase

        dict_in_out_phase_keys = dict_in_out_phase.keys()
        dict_crystal_keys = dict_crystal.keys()

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
                ("multiplicity_hkl" in dict_in_out_phase_keys) and not(flag_unit_cell_parameters or flags_offset_ttheta)):
            index_hkl = dict_in_out_phase["index_hkl"]
            multiplicity_hkl = dict_in_out_phase["multiplicity_hkl"]
        else:
            if flag_phase_texture:
                reduced_symm_elems_p1 = numpy.array([[0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1]], dtype=int)
                translation_elems_p1 = numpy.array([[0], [0], [0], [1]], dtype=int)
                index_hkl, multiplicity_hkl = calc_index_hkl_multiplicity_in_range(
                    sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems_p1, translation_elems_p1, False)
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
        ttheta_hkl = 2*numpy.arcsin(sthovl_hkl*wavelength)
        dict_in_out_phase["ttheta_hkl"] = ttheta_hkl + offset_ttheta
        if radiation[0].startswith("neutrons"):
            f_nucl, dder_f_nucl = calc_f_nucl_by_dictionary(
                dict_crystal, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
            flag_f_nucl = len(dder_f_nucl.keys()) > 0

            flag_para = False
            if (("atom_para_index" in dict_crystal_keys) and ("atom_para_susceptibility" in dict_crystal_keys)):
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
                    iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_1d_para(
                        f_nucl, tensor_sigma, beam_polarization, flipper_efficiency, magnetic_field,
                        flag_f_nucl=flag_f_nucl, flag_tensor_sigma=flag_tensor_sigma,
                        flag_polarization=flags_beam_polarization, flag_flipper=flags_flipper_efficiency)
            elif not(flag_para) and flag_ordered:
                flag_iint_plus_minus = flag_f_nucl or flag_f_m_perp_o or flags_beam_polarization or flags_flipper_efficiency
                if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minus" in dict_in_out_phase_keys) and
                        flag_use_precalculated_data and not(flag_iint_plus_minus)):
                    iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
                else:
                    iint_plus, dder_plus = calc_powder_iint_1d_ordered(
                        f_nucl, f_m_perp_o,
                        flag_f_nucl=flag_f_nucl and flag_calc_analytical_derivatives,
                        flag_f_m_perp=flag_f_m_perp_o and flag_calc_analytical_derivatives)
                    iint_minus = iint_plus
                    dder_minus = dder_plus
            elif flag_para and flag_ordered:
                flag_iint_plus_minus = flag_f_nucl or flag_tensor_sigma or flag_f_m_perp_o or flags_beam_polarization or flags_flipper_efficiency
                if (("iint_plus" in dict_in_out_phase_keys) and ("iint_minu" in dict_in_out_phase_keys) and
                        flag_use_precalculated_data and not(flag_iint_plus_minus)):
                    iint_plus, iint_minus = dict_in_out_phase["iint_plus"], dict_in_out_phase["iint_minus"]
                else:
                    iint_plus, iint_minus, dder_plus, dder_minus = calc_powder_iint_1d_mix(
                        f_nucl, f_m_perp_o, tensor_sigma, beam_polarization, flipper_efficiency, magnetic_field,
                        flag_f_nucl=flag_f_nucl and flag_calc_analytical_derivatives,
                        flag_f_m_perp_ordered=flag_f_m_perp_o and flag_calc_analytical_derivatives,
                        flag_tensor_sigma=flag_tensor_sigma and flag_calc_analytical_derivatives,
                        flag_polarization=flags_beam_polarization and flag_calc_analytical_derivatives,
                        flag_flipper=flags_flipper_efficiency and flag_calc_analytical_derivatives)
            else:
                iint_plus = numpy.square(numpy.abs(f_nucl))
                iint_minus = numpy.square(numpy.abs(f_nucl))

            dict_in_out_phase["iint_plus"] = iint_plus
            dict_in_out_phase["iint_minus"] = iint_minus
        elif radiation[0].startswith("X-rays"):
            f_charge, dder_f_charge = calc_f_charge_by_dictionary(
                dict_crystal, wavelength, dict_in_out_phase, flag_use_precalculated_data=flag_use_precalculated_data)
            flag_f_charge = len(dder_f_charge.keys()) > 0

            iint = numpy.square(numpy.abs(f_charge))
            # FIXME: preparation for XMD
            iint_plus = iint 
            iint_minus = iint

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
                preferred_orientation, dder_po = calc_preferred_orientation_pd(
                    index_hkl, texture_g1, texture_g2, texture_axis, unit_cell_parameters, 
                    flag_texture_g1=flag_texture_g1 and flag_calc_analytical_derivatives,
                    flag_texture_g2=flag_texture_g2 and flag_calc_analytical_derivatives,
                    flag_texture_axis=flag_texture_axis and flag_calc_analytical_derivatives)
                dict_in_out_phase["preferred_orientation"] = preferred_orientation
        
        flag_rp = numpy.any(flags_p_resolution) or numpy.any(flags_resolution_parameters)
        
        hh = resolution_parameters + p_resolution
        u, v, w, x, y = hh[0], hh[1], hh[2], hh[3], hh[4]
        
        profile_pv, dder_pv = calc_profile_pseudo_voight(ttheta_zs, ttheta_hkl, u, v, w, p_ig, x, y,
            p_1, p_2, p_3, p_4, 
            flag_ttheta=flags_offset_ttheta,
            flag_ttheta_hkl=flag_ttheta_hkl, flag_u=flag_rp,
            flag_v=flag_rp, flag_w=flag_rp, flag_i_g=flags_p_ig,
            flag_x=flag_rp, flag_y=flag_rp, 
            flag_p_1=flag_asymmetry_parameters, flag_p_2=flag_asymmetry_parameters,
            flag_p_3=flag_asymmetry_parameters, flag_p_4=flag_asymmetry_parameters)
        dict_in_out_phase["profile_pv"] = profile_pv


        # flags_p_scale
        iint_m_plus = iint_plus * multiplicity_hkl
        iint_m_minus = iint_minus * multiplicity_hkl
        lf = calc_lorentz_factor(ttheta_hkl, k=k, cthm=cthm, flag_ttheta=None)[0]
        dict_in_out_phase["iint_plus_with_factors"] = 0.5 * p_scale * lf * iint_m_plus
        dict_in_out_phase["iint_minus_with_factors"] = 0.5 * p_scale * lf * iint_m_minus
        if flag_texture:
            # 0.5 to have the same meaning for the scale factor as in FullProf
            signal_plus = 0.5 * p_scale * lorentz_factor * (profile_pv * (iint_m_plus * preferred_orientation)[na, :]).sum(axis=1) # sum over hkl
            signal_minus = 0.5 * p_scale * lorentz_factor * (profile_pv * (iint_m_minus * preferred_orientation)[na, :]).sum(axis=1) 
            dict_in_out_phase["iint_plus_with_factors"] *= preferred_orientation
            dict_in_out_phase["iint_minus_with_factors"] *= preferred_orientation
        else:
            signal_plus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_plus[na, :]).sum(axis=1) 
            signal_minus = 0.5 * p_scale * lorentz_factor * (profile_pv * iint_m_minus[na, :]).sum(axis=1) 
        
        dict_in_out_phase["signal_plus"] = signal_plus
        dict_in_out_phase["signal_minus"] = signal_minus
        total_signal_plus += signal_plus
        total_signal_minus += signal_minus

    if flag_dict:
        dict_in_out["signal_plus"] = total_signal_plus
        dict_in_out["signal_minus"] = total_signal_minus

    if ("signal_exp_plus" in dict_pd_keys) and ("signal_exp_minus" in dict_pd_keys):
        signal_exp_plus = dict_pd["signal_exp_plus"]
        signal_exp_minus = dict_pd["signal_exp_minus"]
        if flag_dict:
            dict_in_out["signal_exp_plus"] = signal_exp_plus 
            dict_in_out["signal_exp_minus"] = signal_exp_minus 
        flag_chi_sq_sum, flag_chi_sq_difference = True, True

        if "flag_chi_sq_sum" in dict_pd_keys:
            flag_chi_sq_sum = dict_pd["flag_chi_sq_sum"]

        if "flag_chi_sq_difference" in dict_pd_keys:
            flag_chi_sq_difference = dict_pd["flag_chi_sq_difference"]

        if flag_chi_sq_sum:
            signal_exp = signal_exp_plus[0, :] + signal_exp_minus[0, :]
            signal_sigma = numpy.sqrt(numpy.square(signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))

    else:
        signal_exp = dict_pd["signal_exp"][0,:]
        signal_sigma = dict_pd["signal_exp"][1,:]
        if flag_dict:
            dict_in_out["signal_exp"] = dict_pd["signal_exp"] 
        flag_chi_sq_sum = True
        flag_chi_sq_difference = False
        

    chi_sq = 0.
    n_point = 0
    if flag_chi_sq_sum:
        in_points = numpy.logical_not(excluded_points)
        total_signal_sum = total_signal_plus + total_signal_minus + signal_background
        chi_sq_sum = ((numpy.square((signal_exp - total_signal_sum)/signal_sigma)*in_points)).sum(axis=0)
        chi_sq += chi_sq_sum
        n_point += numpy.sum(in_points)
    
    if flag_chi_sq_difference:
        signal_exp_diff = signal_exp_plus[0, :] - signal_exp_minus[0, :]
        signal_sigma_diff = numpy.sqrt(numpy.square(signal_exp_plus[1, :]) + numpy.square(signal_exp_minus[1, :]))
        total_signal_diff = total_signal_plus - total_signal_minus
        chi_sq_diff = (numpy.square((signal_exp_diff - total_signal_diff)/signal_sigma_diff)).sum(axis=0)
        chi_sq += chi_sq_diff
        n_point += signal_exp_diff.shape[0]
    if numpy.isnan(chi_sq):
        chi_sq = 1e30


    flags_pd = get_flags(dict_pd)
    l_flags_crystal = [get_flags(dict_crystal) for dict_crystal in dict_crystals]

    l_parameter_name = []
    for way, flags in flags_pd.items():
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
    
    if flag_background_intensity:
        pass

    der_chi_sq = numpy.zeros((len(l_parameter_name), ), dtype=float) 
    dder_chi_sq = numpy.zeros((len(l_parameter_name), len(l_parameter_name)), dtype=float)

    return chi_sq, n_point, der_chi_sq, dder_chi_sq, l_parameter_name
    