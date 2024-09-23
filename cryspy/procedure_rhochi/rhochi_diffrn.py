import numpy

from cryspy.A_functions_base.unit_cell import \
    calc_volume_uc_by_unit_cell_parameters, calc_eq_ccs_by_unit_cell_parameters, \
    calc_sthovl_by_unit_cell_parameters

from cryspy.A_functions_base.structure_factor import \
    calc_f_nucl_by_dictionary, calc_sft_ccs_by_dictionary, calc_f_m_perp_by_sft

from cryspy.A_functions_base.extinction import \
    calc_extinction_sphere

from cryspy.A_functions_base.flip_ratio import \
    calc_intensities_by_structure_factors, calc_flip_ratio_by_iint, \
    calc_asymmetry_by_iint

na = numpy.newaxis



def get_flags(obj_dict):
    obj_dict_keys = obj_dict.keys()
    dict_out_variable = {}
    dict_out_refinement = {}
    for key in obj_dict_keys:
        if key.startswith("flags_"):
            label = key[key.find("_")+1:]
            flags = obj_dict[key]
            if numpy.any(flags):
                dict_out_variable[(label, )] = obj_dict[label]
                dict_out_refinement[(label, )] = flags
        if isinstance(obj_dict[key], dict):
            data_dict = obj_dict[key]
            data_dict_keys = data_dict.keys()
            for key_2 in data_dict_keys:
                if key_2.startswith("flags_"):
                    label_2 = key_2[key_2.find("_")+1:]
                    flags_2 = data_dict[key_2]
                    if numpy.any(flags_2):
                        dict_out_variable[(key, label_2)] = data_dict[label_2]
                        dict_out_refinement[(key, label_2)] = flags_2
            
    return dict_out_refinement


def calc_chi_sq_for_diffrn_by_dictionary(
        dict_diffrn, dict_crystal, dict_in_out: dict = None,
        flag_use_precalculated_data: bool=False,
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

    dict_diffrn_keys = dict_diffrn.keys()

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    index_hkl = dict_diffrn["index_hkl"]
    flags_unit_cell_parameters = dict_crystal["flags_unit_cell_parameters"]
    flag_unit_cell_parameters = numpy.any(flags_unit_cell_parameters)
    flag_volume_unit_cell = flag_unit_cell_parameters
    volume_unit_cell, dder_vuc = calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
    flag_sthovl = flag_unit_cell_parameters
    magnetic_field = dict_diffrn["magnetic_field"]
    flag_magnetic_field=False

    matrix_u = dict_diffrn["matrix_u"]
    wavelength = dict_diffrn["wavelength"]
    flags_wavelength = dict_diffrn["flags_wavelength"]
    flip_ratio_excluded = dict_diffrn["flip_ratio_excluded"] # it also work for intensity

    flag_unpolarized, flag_asymmetry = False, False
    if "flip_ratio_es" in dict_diffrn_keys:
        flip_ratio_es = dict_diffrn["flip_ratio_es"]
        beam_polarization = dict_diffrn["beam_polarization"]
        flags_beam_polarization = dict_diffrn["flags_beam_polarization"]
        flipper_efficiency = dict_diffrn["flipper_efficiency"]
        flags_flipper_efficiency = dict_diffrn["flags_flipper_efficiency"]
        if "flag_asymmetry" in dict_diffrn_keys:
            flag_asymmetry = dict_diffrn["flag_asymmetry"]
    elif "intensity_es" in dict_diffrn_keys:
        intensity_es = dict_diffrn["intensity_es"]
        phase_scale = dict_diffrn["phase_scale"]
        flags_phase_scale = dict_diffrn["flags_phase_scale"]
        beam_polarization = 1.
        flipper_efficiency = 1.
        flag_unpolarized = True
        flags_beam_polarization, flags_flipper_efficiency = False, False

    if "c_lambda2" in dict_diffrn_keys:
        c_lambda2 = dict_diffrn["c_lambda2"]
        flags_c_lambda2 = dict_diffrn["flags_c_lambda2"]
        index_2hkl = 2*index_hkl
    else:
        c_lambda2 = None
        flags_c_lambda2 = False
        index_2hkl = None

    flag_eq_ccs = flag_unit_cell_parameters
    if (flag_use_precalculated_data and not(flag_eq_ccs) and ("eq_ccs" in dict_in_out_keys)):
        eq_ccs = dict_in_out["eq_ccs"]
    else:
        eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["eq_ccs"] = eq_ccs

    flags_extinction_radius, flags_extinction_mosaicity = False, False
    if "extinction_model" in dict_diffrn_keys:
        extinction_model = dict_diffrn["extinction_model"]
        extinction_radius = dict_diffrn["extinction_radius"]
        flags_extinction_radius = dict_diffrn["flags_extinction_radius"]
        extinction_mosaicity = dict_diffrn["extinction_mosaicity"]
        flags_extinction_mosaicity = dict_diffrn["flags_extinction_mosaicity"]

        flag_sthovl = flag_unit_cell_parameters
        if (flag_use_precalculated_data and not(flag_sthovl) and ("sthovl" in dict_in_out_keys)):
            sthovl = dict_in_out["sthovl"]
        else:
            sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
            if flag_dict:
                dict_in_out["sthovl"] = sthovl
        
        flag_cos_2theta = flag_sthovl or flags_wavelength
        if (flag_use_precalculated_data and not(flag_cos_2theta) and ("cos_2theta" in dict_in_out_keys)):
            cos_2theta = dict_in_out["cos_2theta"]
        else:        
            cos_2theta = 1.-2*numpy.square(sthovl * wavelength)
            if flag_dict:
                dict_in_out["cos_2theta"] = cos_2theta

        flag_extincton = flags_extinction_radius or flags_extinction_mosaicity or flags_wavelength

        def func_extinction(f_sq, flag_f_sq: bool = False):
            return calc_extinction_sphere(
                f_sq, extinction_radius, extinction_mosaicity, volume_unit_cell, cos_2theta, wavelength,
                extinction_model, flag_f_sq=flag_f_sq, flag_radius=flags_extinction_radius,
                flag_mosaicity=flags_extinction_mosaicity,
                flag_volume_unit_cell=flag_volume_unit_cell,
                flag_cos_2theta=flag_cos_2theta,
                flag_wavelength=flags_wavelength)
    else:
        extinction_model = ""
        extinction_radius = None
        extinction_mosaicity = None
        func_extinction = None
        flag_extincton = False
    
    if (flag_use_precalculated_data and (f"dict_in_out_hkl_{dict_crystal['name']:}" in dict_in_out_keys)):
        dict_in_out_crystal_hkl = dict_in_out[f"dict_in_out_hkl_{dict_crystal['name']:}"]
        dict_in_out_crystal_hkl["index_hkl"] = index_hkl
    else:
        dict_in_out_crystal_hkl = {"index_hkl": index_hkl}
        dict_in_out[f"dict_in_out_hkl_{dict_crystal['name']:}"] = dict_in_out_crystal_hkl

    dict_in_out[f"dict_in_out_hkl_{dict_crystal['name']:}"] = dict_in_out_crystal_hkl
    f_nucl, dder_f_nucl = calc_f_nucl_by_dictionary(dict_crystal, dict_in_out_crystal_hkl, flag_use_precalculated_data=flag_use_precalculated_data)
    flag_f_nucl = len(dder_f_nucl.keys()) > 0

    sft_ccs, dder_sft_ccs = calc_sft_ccs_by_dictionary(dict_crystal, dict_in_out_crystal_hkl, flag_use_precalculated_data=flag_use_precalculated_data)
    flag_sft_ccs = len(dder_sft_ccs.keys()) > 0

    flag_mag_atom_susceptibility = "mag_atom_susceptibility" in dder_sft_ccs.keys()

    flag_f_m_perp = flag_sft_ccs or flag_magnetic_field or flag_eq_ccs
    
    if (flag_use_precalculated_data and not(flag_f_m_perp) and "f_m_perp" in dict_in_out_keys):
        f_m_perp = dict_in_out["f_m_perp"]
    else:
        f_m_perp, dder_f_m_perp = calc_f_m_perp_by_sft(
                sft_ccs, magnetic_field, eq_ccs, flag_sft_ccs=flag_sft_ccs, flag_magnetic_field=flag_magnetic_field, flag_eq_ccs=flag_eq_ccs)
        if flag_dict:
            dict_in_out["f_m_perp"] = f_m_perp

    if c_lambda2 is not None:
        if (flag_use_precalculated_data and (f"dict_in_out_2hkl_{dict_crystal['name']:}" in dict_in_out_keys)):
            dict_in_out_crystal_2hkl = dict_in_out[f"dict_in_out_2hkl_{dict_crystal['name']:}"]
            dict_in_out_crystal_2hkl["index_hkl"] = index_2hkl
        else:
            dict_in_out_crystal_2hkl = {"index_hkl": index_2hkl}
            dict_in_out[f"dict_in_out_2hkl_{dict_crystal['name']:}"] = dict_in_out_crystal_2hkl

        f_nucl_2hkl, dder_f_nucl_2hkl = calc_f_nucl_by_dictionary(dict_crystal, dict_in_out_crystal_2hkl, flag_use_precalculated_data=flag_use_precalculated_data)
        flag_f_nucl_2hkl = len(dder_f_nucl_2hkl.keys()) > 0
        if flag_dict:
            dict_in_out["f_nucl_2hkl"] = f_nucl_2hkl

        sft_ccs_2hkl, dder_sft_ccs_2hkl = calc_sft_ccs_by_dictionary(dict_crystal, dict_in_out_crystal_2hkl, flag_use_precalculated_data=flag_use_precalculated_data)
        flag_sft_ccs_2hkl = len(dder_sft_ccs_2hkl.keys()) > 0
        if flag_dict:
            dict_in_out["sft_ccs_2hkl"] = sft_ccs_2hkl

        flag_f_m_perp_2hkl = flag_sft_ccs_2hkl or flag_magnetic_field or flag_eq_ccs
        if (flag_use_precalculated_data and not(flag_f_m_perp) and "f_m_perp_2hkl" in dict_in_out_keys):
            f_m_perp_2hkl = dict_in_out["f_m_perp_2hkl"]
        else:
            f_m_perp_2hkl, dder_f_m_perp_2hkl = calc_f_m_perp_by_sft(
                sft_ccs_2hkl, magnetic_field, eq_ccs, flag_sft_ccs=flag_sft_ccs_2hkl, flag_magnetic_field=flag_magnetic_field, flag_eq_ccs=flag_eq_ccs)
            if flag_dict:
                dict_in_out["f_m_perp_2hkl"] = f_m_perp_2hkl

    else:
        f_nucl_2hkl = None
        flag_f_nucl_2hkl = False
        f_m_perp_2hkl = None
        flag_f_m_perp_2hkl = False


    flag_flip_ratio = flag_extincton or flags_beam_polarization or flags_flipper_efficiency or \
        flag_f_nucl or flag_f_m_perp or flags_c_lambda2 or flag_f_nucl_2hkl or flag_f_m_perp_2hkl
    if (flag_use_precalculated_data and not(flag_flip_ratio) 
            and "iint_plus" in dict_in_out_keys and "iint_minus" in dict_in_out_keys):
        iint_plus = dict_in_out["iint_plus"]
        iint_minus = dict_in_out["iint_minus"]
    else:
        axis_z = matrix_u[6:9]
        iint_plus, iint_minus, dder_plus, dder_minus = calc_intensities_by_structure_factors(
            beam_polarization, flipper_efficiency, f_nucl, f_m_perp, axis_z,
            func_extinction=func_extinction,
            c_lambda2=c_lambda2, f_nucl_2hkl=f_nucl_2hkl, f_m_perp_2hkl=f_m_perp_2hkl,
            flag_beam_polarization=flags_beam_polarization, flag_flipper_efficiency=flags_flipper_efficiency,
            flag_f_nucl=flag_f_nucl, flag_f_m_perp=flag_f_m_perp,
            flag_c_lambda2=flags_c_lambda2,
            flag_f_nucl_2hkl=flag_f_nucl_2hkl, flag_f_m_perp_2hkl=flag_f_m_perp_2hkl,
            dict_in_out=dict_in_out)
        if flag_dict:
            dict_in_out["iint_plus"] = iint_plus
            dict_in_out["iint_minus"] = iint_minus

    if flag_asymmetry:
        asymmetry_e = (flip_ratio_es[0] -1.)/(flip_ratio_es[0] + 1.)
        asymmetry_s = numpy.sqrt(2.)*flip_ratio_es[1] * numpy.sqrt(numpy.square(flip_ratio_es[0]) + 1.)/numpy.square(flip_ratio_es[0] + 1.)
        exp_value = numpy.stack([asymmetry_e, asymmetry_s], axis=0)

        model_exp, dder_model_exp = calc_asymmetry_by_iint(
            iint_plus, iint_minus, c_lambda2=None, iint_2hkl=None,
            flag_iint_plus=True, flag_iint_minus=True, 
            flag_c_lambda2=False, flag_iint_2hkl=False)
    elif flag_unpolarized:
        model_exp = (iint_plus + iint_minus)*phase_scale
        dict_in_out["intensity_calc"] = model_exp
        exp_value = intensity_es
        dder_model_exp = {
            "iint_plus": numpy.ones_like(iint_plus)*phase_scale,
            "iint_minus": numpy.ones_like(iint_minus)*phase_scale,
            }
    else:
        model_exp, dder_model_exp = calc_flip_ratio_by_iint(
            iint_plus, iint_minus, c_lambda2=None, iint_2hkl=None,
            flag_iint_plus=True, flag_iint_minus=True, 
            flag_c_lambda2=False, flag_iint_2hkl=False)
        exp_value = flip_ratio_es

    index_true = numpy.logical_not(flip_ratio_excluded)
    chi_sq = numpy.square((exp_value[0, index_true]-model_exp[index_true])/exp_value[1, index_true]).sum(axis=0)
    der_model_int_plus = dder_model_exp["iint_plus"]
    der_model_int_minus = dder_model_exp["iint_minus"]

    n_hkl = index_true.sum()
    if flag_dict:
        dict_in_out["chi_sq"] = chi_sq
        dict_in_out["n_hkl"] = n_hkl

    diff_exp_model = (exp_value[0, index_true]-model_exp[index_true])
    inv_sigma_sq_exp_value = 1./numpy.square(exp_value[1, index_true])
    
    dder_plus_crystal, dder_minus_crystal = {}, {}
    if flag_mag_atom_susceptibility:
        dder_f_m_perp_mas = (
            dder_f_m_perp['sft_ccs_real'][:,:, na, :, na]*dder_sft_ccs["mag_atom_susceptibility"].real[na, :, :, : , :]+
            dder_f_m_perp['sft_ccs_imag'][:,:, na, :, na]*dder_sft_ccs["mag_atom_susceptibility"].imag[na, :, :, : , :]).sum(axis=1)
        dder_plus_mas = (
            dder_plus["f_m_perp_real"][:, na, :, na] * dder_f_m_perp_mas.real + 
            dder_plus["f_m_perp_imag"][:, na, :, na] * dder_f_m_perp_mas.imag).sum(axis=0)
        dder_minus_mas = (
            dder_minus["f_m_perp_real"][:, na, :, na] * dder_f_m_perp_mas.real + 
            dder_minus["f_m_perp_imag"][:, na, :, na] * dder_f_m_perp_mas.imag).sum(axis=0)
        if c_lambda2 is not None:
            dder_f_m_perp_mas_2hkl = (
                dder_f_m_perp_2hkl['sft_ccs_real'][:,:, na, :, na]*dder_sft_ccs_2hkl["mag_atom_susceptibility"].real[na, :, :, : , :]+
                dder_f_m_perp_2hkl['sft_ccs_imag'][:,:, na, :, na]*dder_sft_ccs_2hkl["mag_atom_susceptibility"].imag[na, :, :, : , :]).sum(axis=1)
            dder_plus_mas += (
            dder_plus["f_m_perp_2hkl_real"][:, na, :, na] * dder_f_m_perp_mas_2hkl.real + 
            dder_plus["f_m_perp_2hkl_imag"][:, na, :, na] * dder_f_m_perp_mas_2hkl.imag).sum(axis=0)
            dder_minus_mas += (
            dder_minus["f_m_perp_2hkl_real"][:, na, :, na] * dder_f_m_perp_mas_2hkl.real + 
            dder_minus["f_m_perp_2hkl_imag"][:, na, :, na] * dder_f_m_perp_mas_2hkl.imag).sum(axis=0)
        dder_plus_crystal["mag_atom_susceptibility"] = numpy.swapaxes(dder_plus_mas, 0, 1)
        dder_minus_crystal["mag_atom_susceptibility"] = numpy.swapaxes(dder_minus_mas, 0, 1)
        
    dder_plus_diffrn, dder_minus_diffrn = {}, {}
    if flags_beam_polarization:
        dder_plus_diffrn["beam_polarization"] = dder_plus["beam_polarization"][:, na]
        dder_minus_diffrn["beam_polarization"] = dder_minus["beam_polarization"][:, na]
    if flags_flipper_efficiency:
        dder_minus_diffrn["flipper_efficiency"] = dder_minus["flipper_efficiency"][:, na]
    if flags_extinction_radius:
        dder_plus_diffrn["extinction_radius"] = dder_plus["radius"][:, na]
        dder_minus_diffrn["extinction_radius"] = dder_minus["radius"][:, na]
    if flags_extinction_mosaicity:
        dder_plus_diffrn["extinction_mosaicity"] = dder_plus["mosaicity"][:, na]
        dder_minus_diffrn["extinction_mosaicity"] = dder_minus["mosaicity"][:, na]
    if flags_c_lambda2:
        dder_plus_diffrn["c_lambda2"] = dder_plus["c_lambda2"][:, na]
        dder_minus_diffrn["c_lambda2"] = dder_minus["c_lambda2"][:, na]

    dder_plus_diffrn_keys = dder_plus_diffrn.keys()
    dder_minus_diffrn_keys = dder_minus_diffrn.keys()
    dder_plus_crystal_keys = dder_plus_crystal.keys()
    dder_minus_crystal_keys = dder_minus_crystal.keys()
    l_dder_plus_p = []
    l_dder_minus_p = []
    l_parameter_name = []
    der_chi_sq, dder_chi_sq = numpy.array([], dtype=float), numpy.array([[]], dtype=float)
    if True: # len(dder_plus_diffrn_keys) + len(dder_minus_diffrn_keys) + len(dder_plus_crystal_keys) + len(dder_minus_crystal_keys) > 0
        # Cacluation first and second derivatives over refinement parameters
        diffrn_type_name = dict_diffrn["type_name"]
        crystal_type_name = dict_crystal["type_name"]
        flags_diffrn = get_flags(dict_diffrn) 
        flags_crystal = get_flags(dict_crystal) 

        for way, flags in flags_diffrn.items():
            ind_1d = numpy.atleast_1d(numpy.argwhere(flags)) #.flatten()
            name = way[0]
            if ((name in dder_plus_diffrn_keys) and (name in dder_minus_diffrn_keys)):
                dder_plus_p = dder_plus_diffrn[name][:, flags] # numpy.take(dder_plus_diffrn[name], ind_1d, axis=1)
                dder_minus_p = dder_minus_diffrn[name][:, flags]
            elif name in dder_plus_diffrn_keys:
                dder_plus_p = dder_plus_diffrn[name][:, flags]
                dder_minus_p = numpy.zeros_like(dder_plus_p)
            elif name in dder_minus_diffrn_keys:
                dder_minus_p = dder_minus_diffrn[name][:, flags]
                dder_plus_p = numpy.zeros_like(dder_minus_p)
            elif name == "phase_scale":#FIXME
                dder_plus_p = numpy.zeros(shape=(index_true.size, 1), dtype=float)
                dder_minus_p = numpy.zeros(shape=(index_true.size, 1), dtype=float)
            else:
                raise AttributeError("It should not be like this.")
            parameter_name = [(diffrn_type_name, ) + way + (tuple(ind_1d[ind,:]), ) for ind in range(ind_1d.shape[0])]
            l_dder_plus_p.append(dder_plus_p)
            l_dder_minus_p.append(dder_minus_p)
            l_parameter_name.extend(parameter_name)

        for way, flags in flags_crystal.items():
            flag_der_hh = True
            ind_1d = numpy.atleast_1d(numpy.argwhere(flags)) #.flatten()
            name = way[0]
            if ((name in dder_plus_crystal_keys) and (name in dder_minus_crystal_keys)):
                dder_plus_p = dder_plus_crystal[name][:, flags]
                dder_minus_p = dder_minus_crystal[name][:, flags]
            elif name in dder_plus_crystal_keys:
                dder_plus_p = dder_plus_crystal[name][:, flags]
                dder_minus_p = numpy.zeros_like(dder_plus_p)
            elif name in dder_minus_crystal_keys:
                dder_minus_p = dder_minus_crystal[name][:, flags]
                dder_plus_p = numpy.zeros_like(dder_minus_p)
            else:
                flag_der_hh = False # raise AttributeError("It should not be like this.")

            parameter_name = [(crystal_type_name, ) + way + (tuple(ind_1d[ind,:]), ) for ind in range(ind_1d.shape[0])]
            l_parameter_name.extend(parameter_name)
            if flag_der_hh:
                l_dder_plus_p.append(dder_plus_p)
                l_dder_minus_p.append(dder_minus_p)

        if len(l_dder_plus_p) > 0:
            dder_plus_p = numpy.concatenate(l_dder_plus_p, axis=1)
            dder_minus_p = numpy.concatenate(l_dder_minus_p, axis=1)

            der_model_p = der_model_int_plus[:, na]*dder_plus_p + der_model_int_minus[:, na]*dder_minus_p
            der_chi_sq = -2.* (diff_exp_model[:, na] * inv_sigma_sq_exp_value[:, na] * der_model_p).sum(axis=0)
            dder_chi_sq = 2.* (inv_sigma_sq_exp_value[:, na, na]*(der_model_p[:, na, :]*der_model_p[:, :, na])).sum(axis=0) #w_diff is equal to zero

            # iint_minus_sq = numpy.square(iint_minus)
            # der_fr_p = (iint_minus[:, na] * dder_plus_p - iint_plus[:, na] * dder_minus_p) / iint_minus_sq[:, na]
            
            # dder_fr_p = (2*iint_minus[:, na, na] * dder_minus_p[:, :, na] * dder_minus_p[:, na, :] - 
            #    iint_minus_sq[:, na, na] * dder_plus_p[:, :, na] * dder_minus_p[:, na, :] - 
            #    dder_minus_p[:, :, na] * dder_plus_p[:, na, :]) / numpy.square(iint_minus_sq)[:, na, na]
    
            # der_chi_sq = -2.* (diff_exp_model[:, na] * inv_sigma_sq_exp_value[:, na] * der_fr_p).sum(axis=0)
            # dder_chi_sq = 2.* (inv_sigma_sq_exp_value[:, na, na]*(der_fr_p[:, na, :]*der_fr_p[:, :, na] - 
            #     0*diff_exp_model[:, na, na]*dder_fr_p)).sum(axis=0) #w_diff is equal to zero

    return chi_sq, n_hkl, der_chi_sq, dder_chi_sq, l_parameter_name
    