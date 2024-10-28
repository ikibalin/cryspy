import numpy
from .matrix_operations import calc_m1_m2_inv_m1, \
    calc_m_q_inv_m, calc_m_v, calc_vector_product_v1_v2_v1,\
    calc_mm_as_m_q_inv_m, \
    calc_mm_as_m1_m2_inv_m1, \
    calc_m_v, calc_m1_m2_inv_m1
from .symmetry_elements import calc_multiplicity_by_atom_symm_elems, calc_symm_flags
from .unit_cell import calc_eq_ccs_by_unit_cell_parameters, calc_m_m_by_unit_cell_parameters,\
    calc_volume_uc_by_unit_cell_parameters, \
    calc_sthovl_by_unit_cell_parameters
from .local_susceptibility import calc_chi_direct_norm_with_symmetry, calc_m_r_inv_m

from .symmetry_constraints import calc_sc_chi_full

from .extinction import calc_extinction_sphere
from .flip_ratio import \
    calc_iint, calc_flip_ratio_by_iint, \
    calc_asymmetry_by_iint

na = numpy.newaxis


def calc_symm_elem_points_by_index_points(index_point, n_abc):
    common_denominator = numpy.lcm.reduce(n_abc)
    coeff = numpy.floor_divide(numpy.lcm.reduce(n_abc), n_abc)
    denom = common_denominator*numpy.ones_like(index_point[0])
    symm_elem_points = numpy.stack([index_point[0]*coeff[0], index_point[1]*coeff[1], index_point[2]*coeff[2], denom], axis=0)
    return symm_elem_points


def calc_mem_col(
        index_hkl, unit_cell_parameters, eh_ccs, full_symm_elems, symm_elem_points,
        volume_unit_cell, number_unit_cell,
        point_multiplicity=None, dict_in_out: dict=None, flag_use_precalculated_data: bool = False):
    dict_in_out_keys = dict_in_out.keys()
    if (("index_hkl" in dict_in_out_keys) and flag_use_precalculated_data):
        if numpy.any(dict_in_out["index_hkl"] != index_hkl):
            flag_use_precalculated_data = False
            dict_in_out["index_hkl"] = index_hkl
    else:
        dict_in_out["index_hkl"] = index_hkl

    if (("eq_ccs" in dict_in_out_keys) and flag_use_precalculated_data):
        eq_ccs = dict_in_out["eq_ccs"]
    else:
        eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        dict_in_out["eq_ccs"] = eq_ccs

    if point_multiplicity is None:
        if (("point_multiplicity" in dict_in_out_keys) and flag_use_precalculated_data):
            point_multiplicity = dict_in_out["point_multiplicity"]
        else:
            point_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, symm_elem_points)
            dict_in_out["point_multiplicity"] = point_multiplicity
    else:
        dict_in_out["point_multiplicity"] = point_multiplicity

    eh_perp_ccs = calc_vector_product_v1_v2_v1(eq_ccs, numpy.expand_dims(eh_ccs, axis=1), flag_v1=False, flag_v2=False)[0]
    fract_points = symm_elem_points[:3]/numpy.expand_dims(symm_elem_points[3], axis=0)

    rr = calc_m_v(
        numpy.expand_dims(full_symm_elems[4:13], axis=1),
        numpy.expand_dims(fract_points, axis=2), flag_m=False, flag_v=False)[0]

    fract_b = full_symm_elems[:3]/numpy.expand_dims(full_symm_elems[3], axis=0)
    hh = numpy.expand_dims(rr + numpy.expand_dims(fract_b, axis=1), axis=1)
    index_hkl_3d = numpy.expand_dims(numpy.expand_dims(index_hkl, axis=2), axis=3)
    # [hkl, points, symm]
    phase_3d = numpy.exp(-2.*numpy.pi * 1j*(hh[0] * index_hkl_3d[0] + hh[1] * index_hkl_3d[1] + hh[2] * index_hkl_3d[2]))
    phase_2d = phase_3d.sum(axis=2)
    # [3, hkl, points]
    mem_col = 0.2695*(numpy.expand_dims(phase_2d*numpy.expand_dims(point_multiplicity, axis=0),axis=0) * 
        numpy.expand_dims(eh_perp_ccs, axis=2))/full_symm_elems.shape[1]*volume_unit_cell/number_unit_cell
    return mem_col


def calc_mem_chi(
        index_hkl, unit_cell_parameters, h_ccs, full_symm_elems, symm_elem_points,
        point_susceptibility, volume_unit_cell, number_unit_cell,
        point_multiplicity=None, dict_in_out: dict=None, flag_use_precalculated_data: bool = False):
    dict_in_out_keys = dict_in_out.keys()
    if (("index_hkl" in dict_in_out_keys) and flag_use_precalculated_data):
        if numpy.any(dict_in_out["index_hkl"] != index_hkl):
            flag_use_precalculated_data = False
            dict_in_out["index_hkl"] = index_hkl
    else:
        dict_in_out["index_hkl"] = index_hkl

    if (("eq_ccs" in dict_in_out_keys) and flag_use_precalculated_data):
        eq_ccs = dict_in_out["eq_ccs"]
    else:
        eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        dict_in_out["eq_ccs"] = eq_ccs

    if point_multiplicity is None:
        if (("point_multiplicity" in dict_in_out_keys) and flag_use_precalculated_data):
            point_multiplicity = dict_in_out["point_multiplicity"]
        else:
            point_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, symm_elem_points)
            dict_in_out["point_multiplicity"] = point_multiplicity
    else:
        dict_in_out["point_multiplicity"] = point_multiplicity

    m_r_inv_m = calc_m_r_inv_m(unit_cell_parameters, full_symm_elems, flag_unit_cell_parameters=False)[0]
    point_susceptibility_2d_ccs = calc_m1_m2_inv_m1(
        numpy.expand_dims(m_r_inv_m, axis=1),
        numpy.expand_dims(point_susceptibility, axis=2), flag_m1=False, flag_m2=False)[0]

    m_ccs = calc_m_v(point_susceptibility_2d_ccs, 
        numpy.expand_dims(numpy.expand_dims(h_ccs, axis=1), axis=2), flag_m=False, flag_v=False)[0]

    m_perp_ccs = calc_vector_product_v1_v2_v1(
        numpy.expand_dims(numpy.expand_dims(eq_ccs, axis=2), axis=3), 
        numpy.expand_dims(m_ccs, axis=1), flag_v1=False, flag_v2=False)[0]


    fract_points = symm_elem_points[:3]/numpy.expand_dims(symm_elem_points[3], axis=0)
    rr = calc_m_v(
        numpy.expand_dims(full_symm_elems[4:13], axis=1),
        numpy.expand_dims(fract_points, axis=2), flag_m=False, flag_v=False)[0]

    fract_b = full_symm_elems[:3]/numpy.expand_dims(full_symm_elems[3], axis=0)
    hh = numpy.expand_dims(rr + numpy.expand_dims(fract_b, axis=1), axis=1)
    index_hkl_3d = numpy.expand_dims(numpy.expand_dims(index_hkl, axis=2), axis=3)
    # [hkl, points, symm]
    phase_3d = numpy.exp(-2.*numpy.pi * 1j*(hh[0] * index_hkl_3d[0] + hh[1] * index_hkl_3d[1] + hh[2] * index_hkl_3d[2]))

    # [3, hkl, points]
    mem_chi = 0.2695*((numpy.expand_dims(phase_3d, axis=0) * m_perp_ccs).sum(axis=3))/full_symm_elems.shape[1] * \
        numpy.expand_dims(numpy.expand_dims(point_multiplicity, axis=0), axis=1)*volume_unit_cell/number_unit_cell
    return mem_chi


def renormailize_density_col(
        density_col, point_multiplicity,
        volume_unit_cell, points_unit_cell):
    coeff = float(points_unit_cell)/volume_unit_cell
    prod_den_mult = density_col*point_multiplicity
    norm_density_col = coeff*density_col/numpy.expand_dims(prod_den_mult.sum(axis=1), axis=1)
    return norm_density_col


def get_uniform_density_col(point_multiplicity, volume_unit_cell, points_unit_cell):
    density_col = numpy.ones((2, ) + point_multiplicity.shape, dtype=float)
    norm_density_col = renormailize_density_col(
        density_col, point_multiplicity, volume_unit_cell, points_unit_cell)
    return norm_density_col


def renormailize_density_chi(
        point_density, point_multiplicity, point_atom_label,
        point_atom_multiplicity,
        volume_unit_cell, points_unit_cell):
    na = numpy.newaxis
    coeff = float(points_unit_cell)/volume_unit_cell
    point_norm_density = numpy.copy(point_density)
    prod_den_mult = point_density*point_multiplicity/point_atom_multiplicity
    #FIXME: very slow solution. Redo it.
    atom_label = numpy.unique(point_atom_label)
    flag_2d = point_atom_label[na, :] == atom_label[:, na]

    for flag_1d in flag_2d:
        m_at = numpy.sum(prod_den_mult, where=flag_1d)
        if numpy.any(flag_1d):
            point_norm_density[flag_1d] *= coeff/float(m_at)
    return point_norm_density


def get_uniform_density_chi(point_multiplicity, point_atom_label, point_atom_multiplicity, volume_unit_cell, points_unit_cell):
    density_chi = numpy.ones(point_multiplicity.shape, dtype=float)
    norm_density_col = renormailize_density_chi(
        density_chi, point_multiplicity, point_atom_label,
        point_atom_multiplicity,
        volume_unit_cell, points_unit_cell)
    return norm_density_col


def save_spin_density_into_file(
        f_name: str, index_auc, spin_density, n_abc,
        unit_cell_parameters,
        reduced_symm_elems, translation_elems, centrosymmetry: bool, centrosymmetry_position):
    """Save to file.
    if centrosymmetry is False then centrosymmetry_position can be None
    """
    ls_out = []
    ls_out.append("Collinear density (created by CrysPy)")
    ls_out.append(f"{index_auc.shape[1]:}")

    for ind_xyz, sd in zip(index_auc.transpose(), spin_density.transpose()):
        ls_out.append(f"{ind_xyz[0]:4} {ind_xyz[1]:4} {ind_xyz[2]:4} {sd[0]+sd[1]:15.7f} {sd[0]:15.7f} {sd[1]:15.7f}")

    irad = 180./numpy.pi
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c, al = unit_cell_parameters[2], unit_cell_parameters[3]*irad
    be, ga = unit_cell_parameters[4]*irad, unit_cell_parameters[5]*irad

    ls_out.append(
        f"{a:10.5f}{b:10.5f}{c:10.5f}{al:10.5f}{be:10.5f}{ga:10.5f}")

    n_rs = reduced_symm_elems.shape[1]
    n_transl = translation_elems.shape[1]
    centr = int(centrosymmetry)
    ls_out.append(f"{n_abc[0]:5}{n_abc[1]:5}{n_abc[2]:5}{n_rs:5}{centr:5}{n_transl:5}")

    for s_e in reduced_symm_elems.transpose():
        r_11, r_12, r_13 = s_e[4], s_e[5], s_e[6]
        r_21, r_22, r_23 = s_e[7], s_e[8], s_e[9]
        r_31, r_32, r_33 = s_e[10], s_e[11], s_e[12]
        b_1, b_2, b_3 = s_e[0]/s_e[3], s_e[1]/s_e[3], s_e[2]/s_e[3]

        ls_out.append(
            f"{r_11:4}{r_21:4}{r_31:4}  {r_12:4}{r_22:4}{r_32:4}  \
{r_13:4}{r_23:4}{r_33:4}    {b_1:8.5f}{b_2:8.5f}{b_3:8.5f}")

    for orig in translation_elems.transpose():
        ls_out.append(
            f"{float(orig[0]/orig[3]):8.4f}{float(orig[1]/orig[3]):8.4f}\
{float(orig[2]/orig[3]):8.4f}")

    with open(f_name, "w") as fid:
        fid.write("\n".join(ls_out))

    return


def calc_index_atom_symmetry_closest_to_fract_xyz(
        fract_xyz, fract_atom_xyz, full_symm_elems, unit_cell_parameters):
    """
    Calculate index of atoms and applied symmetry to have closest atoms.

    Basins are defined as closest points to atom

    fract_xyz = [3, n_points]
    fract_atom_xyz = [3, n_atoms]

    Output:
        - n_atom_index = [n_points]: integers from 0 until n_atoms
        - n_symmetry = [n_points]: integers from 0 until n_symmetry
    """
    
    fract_xyz = numpy.mod(fract_xyz, 1.)
    fract_atom_xyz = fract_atom_xyz
    
    elem_r = full_symm_elems[4:13]
    fract_b = full_symm_elems[:3]/numpy.expand_dims(full_symm_elems[3], axis=0)
    rr_atom = calc_m_v(
        numpy.expand_dims(elem_r, axis=1),
        numpy.expand_dims(fract_atom_xyz, axis=2), flag_m=False, flag_v=False)[0]
    
    fract_atom_s = numpy.mod(
        rr_atom + numpy.expand_dims(fract_b, axis=1), 1.)
    
    diff_fract = (numpy.expand_dims(numpy.expand_dims(fract_xyz, axis=2), axis=3) -
                  numpy.expand_dims(fract_atom_s, axis=1))

    flag_g_half = numpy.abs(diff_fract) > 0.5
    flag_sign = numpy.where(diff_fract >= 0.0, -1.0, 1.0)
    diff_fract_1 = numpy.where(flag_g_half, diff_fract+flag_sign, diff_fract)
    
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    diff_position_1 = calc_m_v(m_m, diff_fract_1, flag_m=False, flag_v=False)[0]

    dist_sq_3d = numpy.square(diff_position_1).sum(axis=0)
    dist_sq_2d_over_sym = numpy.min(dist_sq_3d, axis=2)
    dist_sq_2d_over_at = numpy.min(dist_sq_3d, axis=1)
    
    ind_at = numpy.argmin(dist_sq_2d_over_sym, axis=1)
    ind_sym = numpy.argmin(dist_sq_2d_over_at, axis=1)
    dist_sq = numpy.min(dist_sq_2d_over_sym, axis=1)
    distance = numpy.sqrt(dist_sq)

    return ind_at, ind_sym, distance

def form_basins(
        symm_elem_auc, 
        full_symm_elems, unit_cell_parameters,
        atom_label, atom_fract_xyz, atom_multiplicity, mag_atom_label):

    # fractions in asymmetric unit cell
    point_fract_xyz_auc = symm_elem_auc[:3, :]/numpy.expand_dims(symm_elem_auc[3, :], axis=0)

    # separate on basins
    ind_at, ind_sym, point_atom_distance = calc_index_atom_symmetry_closest_to_fract_xyz(
            point_fract_xyz_auc, atom_fract_xyz, full_symm_elems, unit_cell_parameters)


    point_atom_label = atom_label[ind_at]
    point_atom_multiplicity = atom_multiplicity[ind_at]

    point_atom_symm_elems = full_symm_elems[:, ind_sym]

    flag_chi = numpy.sum(
        numpy.expand_dims(point_atom_label, axis=1) == numpy.expand_dims(mag_atom_label, axis=0),
        axis=1).astype(bool)

    atom_label_auc_chi = point_atom_label[flag_chi]
    atom_multiplicity_auc_chi = point_atom_multiplicity[flag_chi]
    atom_distance_auc_chi = point_atom_distance[flag_chi]
    atom_symm_elems_auc_chi = point_atom_symm_elems[:, flag_chi]

    return flag_chi, atom_label_auc_chi, atom_multiplicity_auc_chi, atom_distance_auc_chi, atom_symm_elems_auc_chi

def calc_point_susceptibility(
        unit_cell_parameters, atom_symm_elems_chi, atom_label_chi, atom_para_label, atom_para_susceptibility, atom_para_sc_chi,
        full_symm_elems, point_symm_elem):
    atom_para_susceptibility_averaged = (atom_para_sc_chi * numpy.expand_dims(atom_para_susceptibility,axis=0)).sum(axis=1)
    m_r_inv_m = calc_m_r_inv_m(unit_cell_parameters, atom_symm_elems_chi, flag_unit_cell_parameters=False)[0]
    flag_mag = numpy.expand_dims(atom_label_chi, axis=1) == numpy.expand_dims(atom_para_label, axis=0)
    ind_mag_atom = numpy.argmax(flag_mag, axis=1)
    point_susceptibility_a = atom_para_susceptibility_averaged[:, ind_mag_atom]
    point_susceptibility = calc_m_q_inv_m(m_r_inv_m, point_susceptibility_a, flag_m=False, flag_q=False)[0]

    flag_2d = calc_symm_flags(
        numpy.expand_dims(full_symm_elems, axis=2), 
        numpy.expand_dims(point_symm_elem, axis=1))

    l_scm_chi = []
    for ind in range(flag_2d.shape[1]):
        symm_elem_direct = full_symm_elems[:, flag_2d[:, ind]]
        # r_ccs = calc_m_r_inv_m(unit_cell_parameters, symm_elem_direct, flag_unit_cell_parameters=False)[0]
        scm_chi = calc_sc_chi_full(symm_elem_direct, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        l_scm_chi.append(scm_chi)
    point_scm_chi = numpy.stack(l_scm_chi, axis=-1)

    point_susceptibility_aver = (point_scm_chi * numpy.expand_dims(point_susceptibility,axis=0)).sum(axis=1)
    return point_susceptibility_aver


def calc_chi_atoms(
        unit_cell_parameters, number_points, full_symm_elems,
        index_hkl, atom_para_fract_xyz, atom_para_sc_chi, symm_elem_channel_chi,
        point_multiplicity_channel_chi, density_channel_chi):
    
    m_r_inv_m, dder = calc_m_r_inv_m(unit_cell_parameters, full_symm_elems)
    mim_1, dder = calc_mm_as_m1_m2_inv_m1(m_r_inv_m)
    m_m, dder = calc_m_m_by_unit_cell_parameters(unit_cell_parameters)
    
    
    # for SCM
    flag_2d = calc_symm_flags(
        numpy.expand_dims(full_symm_elems, axis=2), 
        numpy.expand_dims(symm_elem_channel_chi, axis=1))
    
    l_scm_chi = []
    for ind in range(flag_2d.shape[1]):
        symm_elem_direct = full_symm_elems[:, flag_2d[:, ind]]
        scm_chi = calc_sc_chi_full(symm_elem_direct, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        l_scm_chi.append(scm_chi)
    point_scm_chi = numpy.stack(l_scm_chi, axis=-1)
    
    point_fract_xyz_auc = symm_elem_channel_chi[:3, :]/numpy.expand_dims(symm_elem_channel_chi[3, :], axis=0)
    ind_at, ind_sym, point_atom_distance = calc_index_atom_symmetry_closest_to_fract_xyz(
            point_fract_xyz_auc, atom_para_fract_xyz, full_symm_elems, unit_cell_parameters)
    
    symm_elems_ai = full_symm_elems[:, ind_sym]
    
    m_r_ai_inv_m, dder = calc_m_r_inv_m(unit_cell_parameters, symm_elems_ai)
    mim_2, dder = calc_mm_as_m_q_inv_m(m_r_ai_inv_m)
    
    rr = calc_m_v(
        numpy.expand_dims(full_symm_elems[4:13], axis=1),
        numpy.expand_dims(point_fract_xyz_auc, axis=2), flag_m=False, flag_v=False)[0]

    fract_b = full_symm_elems[:3]/numpy.expand_dims(full_symm_elems[3], axis=0)
    hh = numpy.expand_dims(rr + numpy.expand_dims(fract_b, axis=1), axis=1)
    index_hkl_3d = numpy.expand_dims(numpy.expand_dims(index_hkl, axis=2), axis=3)
    # [hkl, points, symm]
    phase_3d = numpy.exp(-2.*numpy.pi * 1j*(hh[0] * index_hkl_3d[0] + hh[1] * index_hkl_3d[1] + hh[2] * index_hkl_3d[2]))    
    
    # sum over symmetry
    l_hhh = [(numpy.expand_dims(phase_3d[ind:(ind+1), :, :], axis=(0,1)) *
          numpy.expand_dims(mim_1, axis=(2,3))).sum(axis=4) for ind in range(phase_3d.shape[0])]
    
    hhh = numpy.concatenate(l_hhh, axis=2)
    
    hhh2 = numpy.expand_dims(point_multiplicity_channel_chi*density_channel_chi, axis=(0,1,2)) * hhh 
    
    
    hh1 = (numpy.expand_dims(mim_2, axis=(2,3)) * numpy.expand_dims(atom_para_sc_chi, axis=(0,4))).sum(axis=1)
    hh2 = (numpy.expand_dims(point_scm_chi, axis=(2,3)) * numpy.expand_dims(hh1, axis=0)).sum(axis=1)
    
    # sum over nu
    l_flag_at = [ind_at == i_at for i_at in range(ind_at.max()+1)]
    l_hhh3 = []
    for ind in range(hhh2.shape[2]):
        hhhh3 = ((numpy.expand_dims(hhh2[:, :, ind:(ind+1), :], axis=(3, 4)) * 
                  numpy.expand_dims(hh2, axis=(0, 2)))).sum(axis=1)
        l_hhh4 = [(hhhh3[:, :, :, ind_at, :]*numpy.expand_dims(flag_at, axis=(0,1,2))).sum(axis=3) for ind_at, flag_at in enumerate(l_flag_at)]
        hhh3 = numpy.stack(l_hhh4, axis=3)
        l_hhh3.append(hhh3)
    hhh3 = numpy.concatenate(l_hhh3, axis=1)    
    
    
    volume_uc, dder = calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters)
    coeff = 0.2695 * volume_uc/(full_symm_elems.shape[1]*number_points)
    chi_atoms = coeff*hhh3
    return chi_atoms


def calc_f_m_perp_ccs_by_chi_atoms(chi_atoms, magnetic_field, vp, atom_para_susceptibility):
    m_chi = (chi_atoms*numpy.expand_dims(atom_para_susceptibility, axis=(0,1))).sum(axis=(2,3))
    mv, dder = calc_m_v(m_chi, magnetic_field)
    f_m_perp_chi = (vp * numpy.expand_dims(mv, axis=0)).sum(axis=1)
    return f_m_perp_chi

def calc_model_value_by_precalculated_data(atom_para_susceptibility, unit_cell_parameters, flag_asymmetry, dict_in_out, l_dict_diffrn):
    flag_use_precalculated_data = False
    volume_unit_cell, dder = calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters)
    l_model_value = []
    for dict_diffrn in l_dict_diffrn:
        diffrn_dict_in_out = dict_in_out["dict_in_out_"+dict_diffrn['type_name']]
        vp = diffrn_dict_in_out["vp"]
        chi_atoms = diffrn_dict_in_out["chi_atoms"]
    
        f_nucl = diffrn_dict_in_out["f_nucl"]
    
        index_hkl = dict_diffrn["index_hkl"]
        magnetic_field = dict_diffrn["magnetic_field"]

        beam_polarization = dict_diffrn["beam_polarization"]
        flipper_efficiency = dict_diffrn["flipper_efficiency"]
        matrix_u = dict_diffrn["matrix_u"]
        flip_ratio_es = dict_diffrn["flip_ratio_es"]
    
        f_m_perp_plus_minus = numpy.zeros(index_hkl.shape, dtype=complex)

    
        f_m_perp_chi = calc_f_m_perp_ccs_by_chi_atoms(chi_atoms, magnetic_field, vp, atom_para_susceptibility)
    
        f_m_perp = f_m_perp_chi + f_m_perp_plus_minus


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
            dict_in_out = dict_in_out, flag_use_precalculated_data = flag_use_precalculated_data)

    
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
    model_value = numpy.concatenate(l_model_value, axis=0)
    return model_value





#FIXME delete
def calc_point_chi_direct(
        mag_atom_label, atom_chi, symm_elems, symm_flag_point_to_point,
        point_atom_label, point_atom_symm_elems, unit_cell_parameters):
    """Calculate susceptibility of point i in direct unit cell.
    mag_atom_label[n_mag_atom] is label of magnetic atoms
    atom_chi[6, n_mag_atom] is local susceptibility of magnetic atoms
    symm_elems[13, n_el_symm] is symmetry operations
    symm_flag_point_to_point[] is flag of symmetry operations which transform point to point
    point_atom_label[n_point] is the label of atom to whom the point is belongs
    point_atom_symm_elems[n_point, ] is the symmetry elements 

    Output:
    point_chi_direct[6, n_point] is the susceptibility parameters in direct unit cell
    """

    flag_2d = point_atom_label[:, na] == mag_atom_label[na, :]
    index_susc = numpy.argmax(flag_2d, axis=1)
    chi_recip_norm = atom_chi[:, index_susc]
    p_a_s_el_r = point_atom_symm_elems[4:, :]

    chi_s_direct, dder_ = calc_chi_direct_norm_with_symmetry(chi_recip_norm, p_a_s_el_r)

    symm_elems_r = symm_elems[4:, :]
    chi_s_direct_2d = calc_m1_m2_inv_m1(symm_elems_r[:, :, na], chi_s_direct[:, na, :], flag_m1=False, flag_m2=False)[0]
    hh = symm_flag_point_to_point.astype(int)
    n_at = (hh.sum(axis=0)).astype(float)
    point_chi_direct = (chi_s_direct_2d*hh[na, :, :]).sum(axis=1)/n_at
    
    return point_chi_direct


#FIXME delete
def calc_point_chi_ccs(
        mag_atom_label, atom_chi, symm_elems, symm_flag_point_to_point,
        point_atom_label, point_atom_symm_elems, unit_cell_parameters):
    """Calculate susceptibility of point i in cartesian coordinate system.
    mag_atom_label[n_mag_atom] is label of magnetic atoms
    atom_chi[6, n_mag_atom] is local susceptibility of magnetic atoms
    symm_elems[13, n_el_symm] is symmetry operations
    symm_flag_point_to_point[] is flag of symmetry operations which transform point to point
    point_atom_label[n_point] is the label of atom to whom the point is belongs
    point_atom_symm_elems[n_point, ] is the symmetry elements 

    Output:
    point_chi_direct[6, n_point] is the susceptibility parameters in cartesian coordinate system
    """

    point_chi_ccs = calc_point_chi_direct(
        mag_atom_label, atom_chi, symm_elems, symm_flag_point_to_point,
        point_atom_label, point_atom_symm_elems, unit_cell_parameters)
    # point_chi_ccs, dder_ = calc_susceptibility_ccs(point_chi_direct, unit_cell_parameters)
    return point_chi_ccs


#FIXME delete
def calc_point_moment_ccs_by_point_chi_direct(symm_elems, point_chi_ccs, field_ccs, unit_cell_parameters):
    """Calculate magnetic moments [3, n_points, n_symm_elems]
    """
    chi_s_ccs_2d = calc_m1_m2_inv_m1(symm_elems[4:, :, na], point_chi_ccs[:, na, :], flag_m1=False, flag_m2=False)[0]
    # chi_s_ccs_2d, dder_ = calc_susceptibility_ccs(chi_s_direct_2d, unit_cell_parameters)
    moment_ccs = calc_m_v(chi_s_ccs_2d, field_ccs, flag_m=False, flag_v=False)[0]
    return moment_ccs


