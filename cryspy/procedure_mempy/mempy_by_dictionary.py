import os
import numpy
import scipy
import scipy.optimize

from cryspy.A_functions_base.symmetry_elements import \
    calc_asymmetric_unit_cell_indexes

from cryspy.A_functions_base.mempy import \
    calc_mem_col, \
    calc_mem_chi, \
    calc_symm_elem_points_by_index_points, \
    get_uniform_density_col, \
    renormailize_density_col, \
    save_spin_density_into_file,\
    form_basins,\
    calc_point_susceptibility, \
    get_uniform_density_chi,\
    renormailize_density_chi, \
    calc_model_value_by_precalculated_data, \
    calc_chi_atoms

from cryspy.A_functions_base.unit_cell import \
    calc_volume_uc_by_unit_cell_parameters, \
    calc_sthovl_by_unit_cell_parameters, \
    calc_eq_ccs_by_unit_cell_parameters
    
from cryspy.A_functions_base.structure_factor import \
    calc_f_nucl_by_dictionary
from cryspy.A_functions_base.flip_ratio import \
    calc_iint, calc_flip_ratio_by_iint, \
    calc_asymmetry_by_iint
from cryspy.A_functions_base.extinction import \
    calc_extinction_sphere
from cryspy.A_functions_base.orbital_functions import \
    calc_density_spherical

from cryspy.A_functions_base.matrix_operations import \
    calc_vv_as_v1_v2_v1

from cryspy.A_functions_base.function_1_error_simplex import \
    error_estimation_simplex


def mempy_reconstruction_by_dictionary(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out,
        parameter_lambda:float=1.e-5, iteration_max:int=1000, parameter_lambda_min:float=1.e-9, delta_density:float=1.e-5):
    # **Input information about mem parameters**
    
    print("*******************************************")
    print("MEM reconstruction by CrysPy (module MEMPy)")
    print("*******************************************\n")

    print("MEM iteration parameters")
    print("------------------------")
    print(f" starting lambda parameter:     {parameter_lambda*1e6:.3f}*10^-6")
    print(f" maximal number of iterations:  {iteration_max:}")
    print(f" minimal lambda parameter:      {parameter_lambda_min*1e6:}*10^-6")
    print(f" delta_density:                 {delta_density*1e5:}*10^-5\n")

    dict_in_out_keys = dict_in_out.keys()
    print("Density reconstruction")
    print("----------------------")

    n_abc = dict_mem_parameters["points_abc"]
    print(f"Unit cell is devided on points {n_abc[0]:} x {n_abc[1]:} x {n_abc[2]:}.")
    channel_plus_minus = dict_mem_parameters["channel_plus_minus"]
    channel_chi = dict_mem_parameters["channel_chi"]

    if channel_plus_minus:
        magnetization_plus = dict_mem_parameters["magnetization_plus"]
        magnetization_minus = dict_mem_parameters["magnetization_minus"]
        file_spin_density = dict_mem_parameters["file_spin_density"]
        dict_in_out["magnetization_plus"] = magnetization_plus
        dict_in_out["magnetization_minus"] = magnetization_minus

    if channel_chi:
        flag_uniform_prior_density = dict_mem_parameters["flag_uniform_prior_density"]
        flag_only_magnetic_basins =  dict_mem_parameters["flag_only_magnetic_basins"]
        file_magnetization_density = dict_mem_parameters["file_magnetization_density"]

    flag_asymmetry =  dict_mem_parameters["flag_asymmetry"]
    gof_desired = dict_mem_parameters["gof_desired"]

    # **Input information about crystal**
    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    full_symm_elems = dict_crystal["full_symm_elems"]

    volume_unit_cell = calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=False)[0]

    reduced_symm_elems = dict_crystal["reduced_symm_elems"]
    centrosymmetry = dict_crystal["centrosymmetry"]
    if centrosymmetry:
        centrosymmetry_position = dict_crystal["centrosymmetry_position"]
    else:
        centrosymmetry_position = None
    translation_elems = dict_crystal["translation_elems"]

    atom_label = dict_crystal["atom_label"]
    atom_fract_xyz = dict_crystal["atom_fract_xyz"]
    atom_multiplicity = dict_crystal["atom_multiplicity"]

    if channel_chi:
        atom_para_label = dict_crystal["atom_para_label"]
        atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
        atom_para_sc_chi = dict_crystal["atom_para_sc_chi"]

    # **Index in asymmetric unit cell**
    print("Calculation of asymmetric unit cell...", end="\r")
    index_auc, point_multiplicity = calc_asymmetric_unit_cell_indexes(n_abc, full_symm_elems)
    symm_elem_auc = calc_symm_elem_points_by_index_points(index_auc, n_abc)
    print(f"Number of points in asymmetric unit cell is {index_auc.shape[1]:}.", end="\n")

    # **Basin devision**
    if channel_chi and flag_only_magnetic_basins:
        print("Devision of asymmetric unit cell on bassins...", end="\r")
        flag_atom_para = numpy.any(numpy.expand_dims(atom_label, axis=1) == numpy.expand_dims(atom_para_label, axis=0), axis=1)

        flag_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, atom_distance_auc_chi, atom_symm_elems_auc_chi = \
            form_basins(symm_elem_auc, full_symm_elems, unit_cell_parameters, atom_label[flag_atom_para],
                        atom_fract_xyz[:,flag_atom_para], atom_multiplicity[flag_atom_para], atom_para_label)
        dict_in_out["atom_multiplicity_channel_chi"] = atom_multiplicity_auc_chi
        print(f"Magnetic basins occupy entire unit cell.           \n(flag_only_magnetic_basins: {flag_only_magnetic_basins:})\n")
    elif channel_chi:
        print("Devision of asymmetric unit cell on bassins...", end="\r")
        flag_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, atom_distance_auc_chi, atom_symm_elems_auc_chi = \
            form_basins(symm_elem_auc, full_symm_elems, unit_cell_parameters, atom_label,
                        atom_fract_xyz, atom_multiplicity, atom_para_label)   
        dict_in_out["atom_multiplicity_channel_chi"] = atom_multiplicity_auc_chi
        print(f"Magnetic basins occupy area around magnetic atoms.          \n(flag_only_magnetic_basins: {flag_only_magnetic_basins:})\n")

    if channel_chi:
        index_auc_chi = index_auc[:, flag_chi]
        point_multiplicity_chi = point_multiplicity[flag_chi]
        dict_in_out["point_multiplicity_channel_chi"] = point_multiplicity_chi
        symm_elem_auc_chi = symm_elem_auc[:, flag_chi]
        dict_in_out["symm_elem_channel_chi"] = symm_elem_auc_chi

    if channel_plus_minus and channel_chi:
        flag_col = numpy.logical_not(flag_chi)

        index_auc_col = index_auc[:, flag_col]
        point_multiplicity_col = point_multiplicity[flag_col]
        symm_elem_auc_col = symm_elem_auc[:, flag_col]
        dict_in_out["point_multiplicity_channel_plus_minus"] = point_multiplicity_col
        dict_in_out["symm_elem_channel_plus_minus"] = symm_elem_auc_col

    elif channel_plus_minus:
        index_auc_col = numpy.copy(index_auc)
        point_multiplicity_col = numpy.copy(point_multiplicity)
        symm_elem_auc_col = numpy.copy(symm_elem_auc)
        dict_in_out["point_multiplicity_channel_plus_minus"] = point_multiplicity_col
        dict_in_out["symm_elem_channel_plus_minus"] = symm_elem_auc_col

    print(f"channel_plus_minus: {channel_plus_minus:}")
    print(f"channel_chi:        {channel_chi:}\n")

    if channel_plus_minus:
        print(f"Magnetization of unit cell: {magnetization_plus+magnetization_minus:.3f}  mu_B")
        print(f"(positive channel {magnetization_plus:.3f} mu_B, negative channel {magnetization_minus:.3f}  mu_B)")
        print(f"\nNumber of density points for channel_plus_minus is {index_auc_col.shape[1]}.")
    if channel_chi:
        print(f"Number of density points for channel_chi is {index_auc_chi.shape[1]}.")

    # **Susceptibility tensor $(3\times 3)$ for each point in magnetic basin**
    if channel_chi:
        print("Calculation of restriction on susceptibility...", end="\r")
        point_susceptibility = calc_point_susceptibility(
            unit_cell_parameters, atom_symm_elems_auc_chi, atom_label_auc_chi,
            atom_para_label, atom_para_susceptibility, atom_para_sc_chi, full_symm_elems, symm_elem_auc_chi)
        dict_in_out["susceptibility_channel_chi"] = point_susceptibility
        print(80*" ", end="\r")

    # **Prior density**
    number_unit_cell = numpy.prod(n_abc)
    print("\nCalculation of prior density...         ", end="\r")
    if channel_chi:
        if flag_uniform_prior_density:
            density_chi_prior = get_uniform_density_chi(point_multiplicity_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, volume_unit_cell, number_unit_cell)
            print("Prior density in channel chi is uniform.        ")
        else:
            density_chi_prior = numpy.zeros_like(atom_distance_auc_chi)
            for label in atom_para_label:
                flag_atom = atom_label_auc_chi==label

                dict_shell = dict_crystal[f"shell_{label:}"]
                kappa = float(dict_crystal["mag_atom_kappa"][dict_crystal["mag_atom_label"] == label])
                den_atom = calc_density_spherical(
                    atom_distance_auc_chi[flag_atom], dict_shell["core_population"], dict_shell["core_coeff"], dict_shell["core_zeta"],
                    dict_shell["core_n"], kappa)
                density_chi_prior[flag_atom] = den_atom
            density_chi_prior = renormailize_density_chi(density_chi_prior, point_multiplicity_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, volume_unit_cell, number_unit_cell)
            print("Prior density in channel chi is core.            ")
    if channel_plus_minus:
        density_col_prior = get_uniform_density_col(point_multiplicity_col, volume_unit_cell, number_unit_cell)
        print("Prior density in channel plus-minus is uniform.          ")

    # **Input information about experiments**
    flag_use_precalculated_data = False
    l_exp_value_sigma = []
    l_mem_chi, l_mem_col = [], []
    print(f"Number of experiments is {len(l_dict_diffrn):}.           ")
    for dict_diffrn in l_dict_diffrn:
        if "dict_in_out_"+dict_diffrn["type_name"] in dict_in_out_keys:
            diffrn_dict_in_out = dict_in_out["dict_in_out_"+dict_diffrn["type_name"]]
        else:
            diffrn_dict_in_out = {}
            dict_in_out["dict_in_out_"+dict_diffrn["type_name"]] = diffrn_dict_in_out

        index_hkl = dict_diffrn["index_hkl"]
        h_ccs = dict_diffrn["magnetic_field"]    
        eh_ccs = dict_diffrn["matrix_u"][6:]
        print(f"Preliminary calculation for experiment {dict_diffrn['name']:}...", end="\r")

        diffrn_dict_in_out["index_hkl"] = index_hkl
        diffrn_dict_in_out_keys = diffrn_dict_in_out.keys()
        if channel_plus_minus:
            if "dict_in_out_col" in diffrn_dict_in_out_keys:
                dict_in_out_col = diffrn_dict_in_out["dict_in_out_col"]
            else:
                dict_in_out_col = {}
                diffrn_dict_in_out["dict_in_out_col"] = dict_in_out_col
            
            mem_col = calc_mem_col(
                index_hkl, unit_cell_parameters, eh_ccs, full_symm_elems, symm_elem_auc_col, 
                volume_unit_cell, number_unit_cell, 
                point_multiplicity=point_multiplicity_col,
                dict_in_out=dict_in_out_col, flag_use_precalculated_data=flag_use_precalculated_data)
            
            diffrn_dict_in_out["mem_col"] = mem_col
            l_mem_col.append(mem_col)

        if channel_chi:
            if "dict_in_out_chi" in diffrn_dict_in_out_keys:
                dict_in_out_chi = diffrn_dict_in_out["dict_in_out_chi"]
            else:
                dict_in_out_chi = {}
                diffrn_dict_in_out["dict_in_out_chi"] = dict_in_out_chi
            mem_chi = calc_mem_chi(
                index_hkl, unit_cell_parameters, h_ccs, full_symm_elems, symm_elem_auc_chi,
                point_susceptibility, volume_unit_cell, number_unit_cell, 
                point_multiplicity=point_multiplicity_chi,
                dict_in_out=dict_in_out_chi, flag_use_precalculated_data=flag_use_precalculated_data)
            diffrn_dict_in_out["mem_chi"] = mem_chi
            l_mem_chi.append(mem_chi)

        f_nucl, dder = calc_f_nucl_by_dictionary(
        dict_crystal, diffrn_dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data)
        diffrn_dict_in_out["f_nucl"] = f_nucl


        flip_ratio_es  = dict_diffrn["flip_ratio_es"]
        if flag_asymmetry:
            asymmetry_e = (flip_ratio_es[0] -1.)/(flip_ratio_es[0] + 1.)
            asymmetry_s = numpy.sqrt(2.)*flip_ratio_es[1] * numpy.sqrt(numpy.square(flip_ratio_es[0]) + 1.)/numpy.square(flip_ratio_es[0] + 1.)
            asymmetry_es = numpy.stack([asymmetry_e, asymmetry_s], axis=0)
            l_exp_value_sigma.append(asymmetry_es)
        else:
            l_exp_value_sigma.append(flip_ratio_es)

        

    exp_value_sigma = numpy.concatenate(l_exp_value_sigma, axis=1)    
    if channel_plus_minus:
        mem_col = numpy.concatenate(l_mem_col, axis=1)
    if channel_chi:
        mem_chi = numpy.concatenate(l_mem_chi, axis=1)

    print(f"Total number of reflections is {exp_value_sigma.shape[1]: }.              ")
    if flag_asymmetry:
        print("Density reconstruction is based on asymmetry parameters.")
    else:
        print("Density reconstruction is based on flip ratios.         ")

    # **Preaparation to MEM itertion procedure**
    if channel_plus_minus:
        density_col = numpy.copy(density_col_prior)
        density_col_next = numpy.copy(density_col_prior)
    if channel_chi:
        density_chi = numpy.copy(density_chi_prior)
        density_chi_next = numpy.copy(density_chi_prior)

    # **MEM iteration**
    print("\nMEM iteration procedure")
    print("-----------------------")
    print(f"Desired GoF is {gof_desired:.2f}.")
    c_desired = gof_desired

    c_previous = numpy.inf
    if channel_plus_minus:
        der_c_den_col_previous = numpy.zeros_like(density_col_prior)
    if channel_chi:
        der_c_den_chi_previous = numpy.zeros_like(density_chi_prior)
    iteration = 0
    flag_next = True
    while flag_next:
        iteration += 1

        if channel_plus_minus:
            density_col = numpy.copy(density_col_next)
        if channel_chi:
            density_chi = numpy.copy(density_chi_next)
        l_model_value = []
        l_der_model_den_pm, l_der_model_den_chi = [], []
        for dict_diffrn in l_dict_diffrn:
            diffrn_dict_in_out = dict_in_out["dict_in_out_"+dict_diffrn['type_name']]
            index_hkl = diffrn_dict_in_out["index_hkl"]

            f_m_perp = numpy.zeros(index_hkl.shape, dtype=complex)
            if channel_plus_minus:
                mem_col_exp = diffrn_dict_in_out["mem_col"]
                hh = numpy.expand_dims(numpy.expand_dims(magnetization_plus * density_col[0] + magnetization_minus * density_col[1], axis=0), axis=1)
                f_m_perp_col = (hh*mem_col_exp).sum(axis=2)
                f_m_perp += f_m_perp_col

            if channel_chi:
                mem_chi_exp = diffrn_dict_in_out["mem_chi"]
                f_m_perp_chi = (density_chi*mem_chi_exp).sum(axis=2)
                f_m_perp += f_m_perp_chi


            beam_polarization = dict_diffrn["beam_polarization"]
            flipper_efficiency = dict_diffrn["flipper_efficiency"]
            matrix_u = dict_diffrn["matrix_u"]
            flip_ratio_es = dict_diffrn["flip_ratio_es"]

            f_nucl = diffrn_dict_in_out["f_nucl"]

            wavelength = dict_diffrn["wavelength"]
            sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
            cos_2theta = numpy.cos(2*numpy.arcsin(sthovl*wavelength))

            extinction_model = dict_diffrn["extinction_model"]
            extinction_radius = dict_diffrn["extinction_radius"]
            extinction_mosaicity = dict_diffrn["extinction_mosaicity"]

            func_extinction = lambda f_sq, flag_f_sq: calc_extinction_sphere(
                    f_sq, extinction_radius, extinction_mosaicity, volume_unit_cell, cos_2theta, wavelength,
                    extinction_model, flag_f_sq=False, flag_radius=False,
                    flag_mosaicity=False,
                    flag_volume_unit_cell=False,
                    flag_cos_2theta=False,
                    flag_wavelength=False)

            axis_z = matrix_u[6:9] 
            iint_plus, iint_minus, dder_plus, dder_minus = calc_iint(
                beam_polarization, flipper_efficiency, f_nucl, f_m_perp, axis_z, func_extinction = func_extinction,
                flag_beam_polarization = False, flag_flipper_efficiency = False,
                flag_f_nucl = False, flag_f_m_perp = True,
                dict_in_out = dict_in_out)

            diffrn_dict_in_out["flip_ratio"] = iint_plus/iint_minus

            der_int_plus_fm_perp_real = dder_plus["f_m_perp_real"]
            der_int_plus_fm_perp_imag = dder_plus["f_m_perp_imag"]
            der_int_minus_fm_perp_real = dder_minus["f_m_perp_real"]
            der_int_minus_fm_perp_imag = dder_minus["f_m_perp_imag"]

            if flag_asymmetry:
                model_exp, dder_model_exp = calc_asymmetry_by_iint(
                    iint_plus, iint_minus, c_lambda2=None, iint_2hkl=None,
                    flag_iint_plus=True, flag_iint_minus=True, 
                    flag_c_lambda2=False, flag_iint_2hkl=False)
            else:
                model_exp, dder_model_exp = calc_flip_ratio_by_iint(
                    iint_plus, iint_minus, c_lambda2=None, iint_2hkl=None,
                    flag_iint_plus=True, flag_iint_minus=True, 
                    flag_c_lambda2=False, flag_iint_2hkl=False)

            l_model_value.append(model_exp)
            der_model_int_plus = numpy.expand_dims(dder_model_exp["iint_plus"], axis=0)
            der_model_int_minus = numpy.expand_dims(dder_model_exp["iint_minus"], axis=0)

            if channel_plus_minus:
                der_model_den_pm_exp = (
                    (mem_col_exp.real*numpy.expand_dims(
                    der_model_int_plus*der_int_plus_fm_perp_real +
                    der_model_int_minus*der_int_minus_fm_perp_real, axis=2)
                    ).sum(axis=0) +
                    (mem_col_exp.imag*numpy.expand_dims(
                    der_model_int_plus*der_int_plus_fm_perp_imag +
                    der_model_int_minus*der_int_minus_fm_perp_imag, axis=2)
                    ).sum(axis=0))
                l_der_model_den_pm.append(der_model_den_pm_exp)

            if channel_chi:
                der_model_den_chi_exp = (
                    (mem_chi_exp.real*numpy.expand_dims(
                    der_model_int_plus*der_int_plus_fm_perp_real +
                    der_model_int_minus*der_int_minus_fm_perp_real, axis=2)
                    ).sum(axis=0) +
                    (mem_chi_exp.imag*numpy.expand_dims(
                    der_model_int_plus*der_int_plus_fm_perp_imag +
                    der_model_int_minus*der_int_minus_fm_perp_imag, axis=2)
                    ).sum(axis=0))
                l_der_model_den_chi.append(der_model_den_chi_exp)

        model_value = numpy.concatenate(l_model_value, axis=0)

        diff_value = (exp_value_sigma[0]-model_value)/exp_value_sigma[1]
        c = numpy.square(diff_value).sum(axis=0)/diff_value.shape[0]

        if channel_plus_minus:
            der_model_den_pm = numpy.concatenate(l_der_model_den_pm, axis=0)
            der_c_den_pm = (-2.)/diff_value.shape[0] * (
                numpy.expand_dims((diff_value/exp_value_sigma[1]),axis=1) *
                der_model_den_pm).sum(axis=0)

            der_c_den_col = numpy.stack([magnetization_plus * der_c_den_pm, magnetization_minus * der_c_den_pm], axis=0)

        if channel_chi:
            der_model_den_chi = numpy.concatenate(l_der_model_den_chi, axis=0)
            der_c_den_chi = (-2.)/diff_value.shape[0] * (
                numpy.expand_dims((diff_value/exp_value_sigma[1]),axis=1) *
                der_model_den_chi).sum(axis=0)

        if c > c_previous:
            parameter_lambda = 0.5 * parameter_lambda
            c = c_previous
            if channel_plus_minus:
                density_col = numpy.copy(density_col_previous)
                der_c_den_col = der_c_den_col_previous
            if channel_chi:
                density_chi = numpy.copy(density_chi_previous)
                der_c_den_chi = der_c_den_chi_previous
        else:
            c_previous = c
            parameter_lambda = 1.03 * parameter_lambda
            if channel_plus_minus:
                density_col_previous = numpy.copy(density_col)
                der_c_den_col_previous = der_c_den_col
            if channel_chi:
                density_chi_previous = numpy.copy(density_chi)
                der_c_den_chi_previous = der_c_den_chi

            print(f"Iteration {iteration:5}, lambda {parameter_lambda*1e6:.3f}*10^-6, chi_sq: {c:.2f}              ", end='\r')

        if channel_plus_minus:
            coeff = (parameter_lambda*number_unit_cell/(c_desired*volume_unit_cell))/point_multiplicity_col
            hh = (density_col+delta_density)*numpy.exp(-coeff*der_c_den_col)-delta_density
            hh = numpy.where(hh>0, hh, 0)
            density_col_next = renormailize_density_col(hh, point_multiplicity_col, volume_unit_cell, number_unit_cell)

        if channel_chi:
            coeff = (parameter_lambda*number_unit_cell/(c_desired*volume_unit_cell))*atom_multiplicity_auc_chi/point_multiplicity_chi
            hh = (density_chi+delta_density)*numpy.exp(-coeff*der_c_den_chi)-delta_density
            hh = numpy.where(hh>0, hh, 0)
            density_chi_next = renormailize_density_chi(hh, point_multiplicity_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, volume_unit_cell, number_unit_cell)


        if iteration >= iteration_max:
            flag_next = False
            print(f"Maximal number of iteration is reached ({iteration:}).         ", end='\n')
        if parameter_lambda < parameter_lambda_min:
            flag_next = False
            print(f"Minimal value of parameter lambda {parameter_lambda*1e6:.3f}*10^-6 is reached at iteration {iteration:}.            ", end='\n')
        if c <= c_desired:
            flag_next = False
            print(f"Desired value is reached at iteration {iteration:}.               ", end='\n')

    c_best = c_previous
    print(f"Chi_sq best is {c_best:.2f}")
    if channel_plus_minus:
        density_col_best = numpy.copy(density_col_previous)
        dict_in_out["density_channel_plus_minus"] = density_col_best
    if channel_chi:
        density_chi_best = numpy.copy(density_chi_previous)
        dict_in_out["density_channel_chi"] = density_chi

    # **Save to .den file**
    if channel_plus_minus and (file_spin_density is not None):
        spin_density = density_col_best * numpy.array([[magnetization_plus, ], [magnetization_minus, ]], dtype=float)
        save_spin_density_into_file(file_spin_density, index_auc_col, spin_density, n_abc, unit_cell_parameters,
                reduced_symm_elems, translation_elems, centrosymmetry, centrosymmetry_position)
        print(f"\nReconstructed spin density is written in file '{file_spin_density:}'.")
    if channel_chi and (file_magnetization_density is not None):
        spin_density = numpy.stack([density_chi_best, numpy.zeros_like(density_chi_best)], axis=0)
        save_spin_density_into_file(file_magnetization_density, index_auc_chi, spin_density, n_abc, unit_cell_parameters,
            reduced_symm_elems, translation_elems, centrosymmetry, centrosymmetry_position)
        print(f"\nReconstructed magnetization density is written in file '{file_magnetization_density:}'.")


def mempy_susceptibility_refinement(dict_channel_chi, dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out):
    print("****************************************")
    print("Susceptibility refinement (module MEMPy)")
    print("****************************************")

    number_points = numpy.prod(dict_mem_parameters["points_abc"])
    
    flag_asymmetry =  dict_mem_parameters["flag_asymmetry"]
    channel_plus_minus = dict_mem_parameters["channel_plus_minus"]
    channel_chi = dict_mem_parameters["channel_chi"]
    print(f"Channel plus/minus is {channel_plus_minus:}")
    print("ATTENTION: Channel plus/minus is not taken into account.")
    print(f"Channel chi is {channel_chi:}")

    print(f"Flag asymmetry is {flag_asymmetry:}")
    
    if channel_plus_minus:
        magnetization_plus = dict_mem_parameters["magnetization_plus"]
        magnetization_minus = dict_mem_parameters["magnetization_minus"]
    
    symm_elem_channel_chi = dict_channel_chi["symm_elem_channel_chi"]
    atom_multiplicity_channel_chi = dict_channel_chi["atom_multiplicity_channel_chi"]
    density_channel_chi = dict_channel_chi["density_channel_chi"]
    point_multiplicity_channel_chi = dict_channel_chi["point_multiplicity_channel_chi"]
    
    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    full_symm_elems = dict_crystal["full_symm_elems"]
    atom_fract_xyz = dict_crystal["atom_fract_xyz"]
    atom_para_sc_chi = dict_crystal["atom_para_sc_chi"]
    atom_para_index = dict_crystal["atom_para_index"]
    atom_para_label = dict_crystal["atom_para_label"]
    atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
    flags_atom_para_susceptibility = dict_crystal["flags_atom_para_susceptibility"]
    print(f"Number of refined parameters is {flags_atom_para_susceptibility.sum():}.")
    if flags_atom_para_susceptibility.sum() == 0:
        print("There is no refined susceptibility parameters.")
        return
    
    atom_para_fract_xyz = atom_fract_xyz[:, atom_para_index]
    
    n_atom_para = atom_para_susceptibility.shape[1]
    
    print("Preliminary calculations of chi atoms ...", end="\r")
    
    l_exp_value_sigma = []
    for dict_diffrn in l_dict_diffrn:
        flag_use_precalculated_data = False
        index_hkl = dict_diffrn["index_hkl"] 
        diffrn_dict_in_out = {"index_hkl": index_hkl}
    
        chi_atoms = calc_chi_atoms(
            unit_cell_parameters, number_points, full_symm_elems,
            index_hkl, atom_para_fract_xyz, atom_para_sc_chi, 
            symm_elem_channel_chi, point_multiplicity_channel_chi, density_channel_chi)
    
        diffrn_dict_in_out["chi_atoms"] = chi_atoms
    
        eq_ccs, dder = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)
        vp, dder = calc_vv_as_v1_v2_v1(eq_ccs)
        diffrn_dict_in_out["vp"] = vp
    
        f_nucl, dder = calc_f_nucl_by_dictionary(
            dict_crystal, diffrn_dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data)
        diffrn_dict_in_out["f_nucl"] = f_nucl

        dict_in_out["dict_in_out_"+dict_diffrn['type_name']] = diffrn_dict_in_out
    
        flip_ratio_es  = dict_diffrn["flip_ratio_es"]

        if flag_asymmetry:
            asymmetry_e = (flip_ratio_es[0] -1.)/(flip_ratio_es[0] + 1.)
            asymmetry_s = numpy.sqrt(2.)*flip_ratio_es[1] * numpy.sqrt(numpy.square(flip_ratio_es[0]) + 1.)/numpy.square(flip_ratio_es[0] + 1.)
            asymmetry_es = numpy.stack([asymmetry_e, asymmetry_s], axis=0)
            l_exp_value_sigma.append(asymmetry_es)
        else:
            l_exp_value_sigma.append(flip_ratio_es)
    exp_value_sigma = numpy.concatenate(l_exp_value_sigma, axis=1)


    
    def calc_chi_sq(param):
        atom_para_susceptibility[flags_atom_para_susceptibility] = param
        model_value = calc_model_value_by_precalculated_data(atom_para_susceptibility, unit_cell_parameters, flag_asymmetry, dict_in_out, l_dict_diffrn)
    
        chi_sq = numpy.square((model_value-exp_value_sigma[0])/exp_value_sigma[1]).sum()
        return chi_sq
    
    param_0 = atom_para_susceptibility[flags_atom_para_susceptibility]
    chi_sq_per_n = calc_chi_sq(param_0)/exp_value_sigma.shape[1]
    print(70*" ")
    print("Before susceptibility refinement")
    print("Susceptibility tensor:")
    for ind_at, label in enumerate(atom_para_label):
        print(f"{label:5}  {atom_para_susceptibility[0, ind_at]:.5f} {atom_para_susceptibility[1, ind_at]:.5f} {atom_para_susceptibility[2, ind_at]:.5f} {atom_para_susceptibility[3, ind_at]:.5f} {atom_para_susceptibility[4, ind_at]:.5f} {atom_para_susceptibility[5, ind_at]:.5f}")
        
    print(f"chi_sq_per_n is {chi_sq_per_n:.2f}.")
    print("Minimization procedure ...", end="\r")
    res = scipy.optimize.minimize(calc_chi_sq, param_0, method="Nelder-Mead")
    apss = None
    if "hess_inv" in res.keys():
        hess_inv = res["hess_inv"]
        dict_in_out["hess_inv"] = hess_inv
        sigma_p = numpy.sqrt(numpy.abs(numpy.diag(hess_inv)))
        atom_para_susceptibility_sigma = numpy.zeros_like(atom_para_susceptibility)
        atom_para_susceptibility_sigma[flags_atom_para_susceptibility] = sigma_p
        apss = (atom_para_sc_chi * numpy.expand_dims(atom_para_susceptibility_sigma, axis=0)).sum(axis=1)
        dict_in_out["atom_para_susceptibility_sigma"] = apss
    elif "final_simplex" in res.keys():
        n = exp_value_sigma.shape[1]
        m_error, dist_hh = error_estimation_simplex(
            res["final_simplex"][0], res["final_simplex"][1], calc_chi_sq)
        l_sigma = []
        for i, val_2 in zip(range(m_error.shape[0]), dist_hh):
            # slightly change definition, instead of (n-k) here is n
            error = (abs(m_error[i, i])*1./n)**0.5
            if m_error[i, i] < 0.:
                pass
                # warn("Negative diagonal elements of Hessian.", UserWarning)
            if val_2 > error:
                pass
                # warn("Minimum is not found.", UserWarning)
            l_sigma.append(max(error, val_2))
        sigma_p = numpy.array(l_sigma)
        atom_para_susceptibility_sigma = numpy.zeros_like(atom_para_susceptibility)
        atom_para_susceptibility_sigma[flags_atom_para_susceptibility] = sigma_p
        apss = (atom_para_sc_chi * numpy.expand_dims(atom_para_susceptibility_sigma, axis=0)).sum(axis=1)
        dict_in_out["atom_para_susceptibility_sigma"] = apss
        print(sigma_p)
    print(70*" ")
    chi_sq_per_n = calc_chi_sq(res.x)/exp_value_sigma.shape[1]
    
    atom_para_susceptibility[flags_atom_para_susceptibility] = res.x
    atom_para_susceptibility = (atom_para_sc_chi * numpy.expand_dims(atom_para_susceptibility, axis=0)).sum(axis=1)     
    dict_crystal["atom_para_susceptibility"] = atom_para_susceptibility
    print("After susceptibility refinement")
    print("Susceptibility tensor:")
    for ind_at, label in enumerate(atom_para_label):
        print(f"{label:5}  {atom_para_susceptibility[0, ind_at]:8.5f} {atom_para_susceptibility[1, ind_at]:8.5f} {atom_para_susceptibility[2, ind_at]:8.5f} {atom_para_susceptibility[3, ind_at]:8.5f} {atom_para_susceptibility[4, ind_at]:8.5f} {atom_para_susceptibility[5, ind_at]:8.5f}")
        if apss is not None:
            print(f"sigma  {apss[0, ind_at]:8.5f} {apss[1, ind_at]:8.5f} {apss[2, ind_at]:8.5f} {apss[3, ind_at]:8.5f} {apss[4, ind_at]:8.5f} {apss[5, ind_at]:8.5f}")

    print(f"chi_sq_per_n is {chi_sq_per_n:.2f}.")
    
    print(70*"*")
    print("End of MEMPy procedure for susceptibility refinement")
    print(70*"*")
    return 


def mempy_cycle_density_susceptibility(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out,
        parameter_lambda:float=1.e-5, iteration_max:int=1000, parameter_lambda_min:float=1.e-9, delta_density:float=1.e-5, n_cycle:int=10):
    print(70*"*")
    print("MEMPy: cycle iteration")
    print(70*"*")
    print(f"Number of cycles is {n_cycle:}")
    print(70*" ")
    for i_cycle in range(n_cycle):
        print(f"Cycle {i_cycle+1:}")
        print(len(f"Cycle {i_cycle+1:}")*"-")

        dict_in_out_den = {}
        mempy_reconstruction_by_dictionary(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out_den,
            parameter_lambda=parameter_lambda, iteration_max=iteration_max, parameter_lambda_min=parameter_lambda_min, delta_density=delta_density)
    
        dict_channel_chi = {
            'atom_multiplicity_channel_chi': dict_in_out_den['atom_multiplicity_channel_chi'],
            'point_multiplicity_channel_chi': dict_in_out_den['point_multiplicity_channel_chi'],
            'symm_elem_channel_chi': dict_in_out_den['symm_elem_channel_chi'],
            'susceptibility_channel_chi': dict_in_out_den['susceptibility_channel_chi'],
            'density_channel_chi': dict_in_out_den['density_channel_chi'],
        }
    
        dict_in_out_susc = {} 
        mempy_susceptibility_refinement(dict_channel_chi, dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out_susc)
        print(70*" ")
    
    dict_in_out["dict_in_out_den"] = dict_in_out_den
    dict_in_out["dict_in_out_susc"] = dict_in_out_susc
    return 