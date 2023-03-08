# -*- coding: utf-8 -*-
"""
Expressions for calculations structure factors

For details see documentation.
"""
import numpy

from .matrix_operations import calc_det_m, calc_m1_m2, calc_m1_m2_inv_m1, calc_m_v, calc_vector_product_v1_v2_v1, calc_m_q_inv_m
from .unit_cell import calc_eq_ccs_by_unit_cell_parameters, calc_m_m_by_unit_cell_parameters, calc_m_m_norm_by_unit_cell_parameters, calc_sthovl_by_unit_cell_parameters
from .debye_waller_factor import calc_dwf
from .symmetry_elements import calc_multiplicity_by_atom_symm_elems, calc_full_symm_elems_by_reduced, calc_equivalent_reflections
from .magnetic_form_factor import calc_form_factor
from .local_susceptibility import calc_m_r_inv_m

na = numpy.newaxis

def get_atom_symm_elems_by_atom_fract_xyz(atom_fract_xyz):
    ones = numpy.ones_like(atom_fract_xyz[0]).astype(int)
    atom_symm_elems = numpy.stack([
        (numpy.round(atom_fract_xyz[0]*3*10**6, decimals=0)).astype(int),
        (numpy.round(atom_fract_xyz[1]*3*10**6, decimals=0)).astype(int),
        (numpy.round(atom_fract_xyz[2]*3*10**6, decimals=0)).astype(int),
        ones*3*10**6], axis=0)
    atom_symm_elems[atom_symm_elems==999999] = 1000000 # 1/3
    atom_symm_elems[atom_symm_elems==1999998] = 2000000 # 2/3 as 0.666666
    atom_symm_elems[atom_symm_elems==2000001] = 2000000 # 2/3 as 0.666667
    return atom_symm_elems

 
def calc_f_m_perp_by_sft(
        sft_ccs, magnetic_field, eq_ccs,
        flag_sft_ccs: bool = False,
        flag_magnetic_field: bool = False,
        flag_eq_ccs: bool = False):
    """Calculate perpendicular component of magnetic structure factor by susceptibility factor tensor.
    All parameters are defined in Cartesian coordinate system (x||a*, z||c).
    """
    f_m, dder_f_m = calc_m_v(
        sft_ccs, magnetic_field, flag_m=flag_sft_ccs, flag_v=flag_magnetic_field)
    if flag_sft_ccs:
        dder_f_m["sft_ccs_real"] = dder_f_m.pop("m_real")
        dder_f_m["sft_ccs_imag"] = dder_f_m.pop("m_imag")
    if flag_magnetic_field:
        dder_f_m["magnetic_field"] = dder_f_m.pop("v")

    
    flag_f_m = flag_sft_ccs or flag_magnetic_field
    f_m_perp, dder_f_m_perp = calc_vector_product_v1_v2_v1(
        eq_ccs, f_m, flag_v1=flag_eq_ccs, flag_v2=flag_f_m)
    if flag_eq_ccs:
        dder_f_m_perp["eq_ccs"] = dder_f_m_perp.pop("v1")
    if flag_f_m:
        dder_f_m_perp["f_m_real"] = dder_f_m_perp.pop("v2_real")
        dder_f_m_perp["f_m_imag"] = dder_f_m_perp.pop("v2_imag")

    dder = {}
    if flag_sft_ccs:
        dder["sft_ccs_real"] = (
            numpy.expand_dims(dder_f_m_perp["f_m_real"], axis=2)*numpy.expand_dims(dder_f_m["sft_ccs_real"].real, axis=0) +
            numpy.expand_dims(dder_f_m_perp["f_m_imag"], axis=2)*numpy.expand_dims(dder_f_m["sft_ccs_real"].imag, axis=0)).sum(axis=1)
        dder["sft_ccs_imag"] = (
            numpy.expand_dims(dder_f_m_perp["f_m_real"], axis=2)*numpy.expand_dims(dder_f_m["sft_ccs_imag"].real, axis=0) +
            numpy.expand_dims(dder_f_m_perp["f_m_imag"], axis=2)*numpy.expand_dims(dder_f_m["sft_ccs_imag"].imag, axis=0)).sum(axis=1)
    if flag_magnetic_field:
        dder["magnetic_field"] = (
            numpy.expand_dims(dder_f_m_perp["f_m_real"], axis=2)*numpy.expand_dims(dder_f_m["magnetic_field"].real, axis=0) +
            numpy.expand_dims(dder_f_m_perp["f_m_imag"], axis=2)*numpy.expand_dims(dder_f_m["magnetic_field"].imag, axis=0)).sum(axis=1)
    if flag_eq_ccs:
        dder["eq_ccs"] = dder_f_m_perp["eq_css"]
    return f_m_perp, dder    


def calc_pr1(index_hkl, reduced_symm_elems, fract_xyz, flag_fract_xyz: bool = False):
    """Calculate PR1, dimensions [hkl, rs, atoms].
    For more details see documentation module "Structure factor".
    """
    index_hkl_exp = numpy.expand_dims(numpy.expand_dims(index_hkl, axis=2), axis=3)
    h, k, l = index_hkl_exp[0], index_hkl_exp[1], index_hkl_exp[2]
    reduced_symm_elems_exp = numpy.expand_dims(numpy.expand_dims(reduced_symm_elems, axis=1), axis=3)
    r_11, r_12, r_13 = reduced_symm_elems_exp[4], reduced_symm_elems_exp[5], reduced_symm_elems_exp[6]
    r_21, r_22, r_23 = reduced_symm_elems_exp[7], reduced_symm_elems_exp[8], reduced_symm_elems_exp[9]
    r_31, r_32, r_33 = reduced_symm_elems_exp[10], reduced_symm_elems_exp[11], reduced_symm_elems_exp[12]
    fract_xyz_exp = numpy.expand_dims(numpy.expand_dims(fract_xyz, axis=1), axis=2)
    x, y, z = fract_xyz_exp[0], fract_xyz_exp[1], fract_xyz_exp[2]

    hh = h*(r_11*x + r_12*y + r_13*z) + k*(r_21*x + r_22*y + r_23*z) + l*(r_31*x + r_32*y + r_33*z)
    res = numpy.exp(-2.*numpy.pi*1j*hh)
    dder = {}
    if flag_fract_xyz:
        dder["fract_xyz"] = numpy.stack([
            numpy.exp(-2.*numpy.pi*1j*(h*r_11 + k*r_21 + l*r_31)),
            numpy.exp(-2.*numpy.pi*1j*(h*r_11 + k*r_21 + l*r_31)),
            numpy.exp(-2.*numpy.pi*1j*(h*r_11 + k*r_21 + l*r_31))], axis=0)
    return res, dder


def calc_pr2(index_hkl, reduced_symm_elems):
    """Calculate PR2, dimensions, dimensions [hkl, rs].
    For more details see documentation module "Structure factor".
    """
    index_hkl_exp = numpy.expand_dims(index_hkl, axis=2)
    h, k, l = index_hkl_exp[0], index_hkl_exp[1], index_hkl_exp[2]
    reduced_symm_elems_exp = numpy.expand_dims(reduced_symm_elems, axis=1)
    b_1, b_2, b_3, b_d = reduced_symm_elems_exp[0], reduced_symm_elems_exp[1], reduced_symm_elems_exp[2], reduced_symm_elems_exp[3]
    hh = h*(b_1.astype(float)/b_d) + k*(b_2.astype(float)/b_d) + l*(b_3.astype(float)/b_d)
    res = numpy.exp(-2.*numpy.pi*1j*hh)
    return res


def calc_pr3(index_hkl, translation_elems):
    """Calculate PR3, dimensions [hkl,].
    For more details see documentation module "Structure factor".
    """
    index_hkl_exp = numpy.expand_dims(index_hkl, axis=2)
    h, k, l = index_hkl_exp[0], index_hkl_exp[1], index_hkl_exp[2]
    translation_elems_exp = numpy.expand_dims(translation_elems, axis=1)
    t_1, t_2, t_3, t_d = translation_elems_exp[0], translation_elems_exp[1], translation_elems_exp[2], translation_elems_exp[3] 
    hh = (h*t_1+k*t_2+l*t_3).astype(float)
    res =(numpy.exp(-2.*numpy.pi*1j*hh/t_d)).sum(axis=1)/translation_elems.shape[-1]
    return res


def calc_pr4(index_hkl, centrosymmetry_position=None):
    """Calculate PR4.
    For more details see documentation module "Structure factor".
    """
    h, k, l = index_hkl[0], index_hkl[1], index_hkl[2]
    if centrosymmetry_position is None:
        res = numpy.zeros_like(h)
    else:
        p_1, p_2, p_3 = centrosymmetry_position[0]/centrosymmetry_position[3], centrosymmetry_position[1]/centrosymmetry_position[3], centrosymmetry_position[2]/centrosymmetry_position[3]
        res = numpy.exp(-4.*numpy.pi * 1j * (h*p_1 + k*p_2 + l*p_3))
    return res


def calc_pr5(reduced_symm_elems, unit_cell_parameters, flag_unit_cell_parameters: bool=False):
    """Calculate PR5, dimensions [rs,].
    For more details see documentation module "Structure factor".
    """
    res, dder = calc_m_r_inv_m(unit_cell_parameters, reduced_symm_elems, flag_unit_cell_parameters=flag_unit_cell_parameters)
    return res, dder


def calc_f_asym_a_by_pr(
        atom_multiplicity, debye_waller, atom_occupancy, pr_1, pr_2,
        flag_debye_waller: bool = False, flag_atom_occupancy: bool = False, flag_pr_1: bool = False):
    """Calculate preliminary asymmetric structure factor by preliminary defined parameters.
    For more details see documentation module "Structure factor".
    """
    # dimension of pr_1 is [hkl, symmetry, a]

    res = (pr_2[:, :, na]*(pr_1*atom_multiplicity[na, na, :]*atom_occupancy[na, na, :]*debye_waller)).sum(axis=1)/pr_2.shape[-1]
    dder = {}
    # if flag_scat_length_neutron:
    #     dder["scat_length_neutron_real"] = (numpy.expand_dims(pr_2, axis=2)*\
    #         (pr_1*atom_multiplicity*atom_occupancy*debye_waller)).sum(axis=1)/pr_2.shape[-1] # FIXME: only for neutron diffraction
    #     dder["scat_length_neutron_imag"] = 1j*(numpy.expand_dims(pr_2, axis=2)*\
    #         (pr_1*atom_multiplicity*atom_occupancy*debye_waller)).sum(axis=1)/pr_2.shape[-1] # FIXME: only for neutron diffraction
    if flag_debye_waller:
        dder["debye_waller"] = (numpy.expand_dims(pr_2, axis=2)*\
            (pr_1*atom_multiplicity*atom_occupancy))/pr_2.shape[-1]
    if flag_atom_occupancy:
        dder["atom_occupancy"] = (numpy.expand_dims(pr_2, axis=2)*\
            (pr_1*atom_multiplicity*debye_waller))/pr_2.shape[-1]
    if flag_pr_1:
        dder["pr_1_real"] = (pr_2[:,:,na]*atom_multiplicity[na, na, :]*atom_occupancy[na, na, :]*debye_waller)/pr_2.shape[-1]
        dder["pr_1_imag"] = 1j*(pr_2[:,:,na]*atom_multiplicity[na, na, :]*atom_occupancy[na, na, :]*debye_waller)/pr_2.shape[-1]
    return res, dder

# Delete IT
# def calc_f_asym_by_pr(
#         atom_multiplicity, scat_length_neutron, debye_waller, atom_occupancy, pr_1, pr_2,
#         flag_scat_length_neutron: bool = False, flag_debye_waller: bool = False,
#         flag_atom_occupancy: bool = False, flag_pr_1: bool = False):
#     """Calculate preliminary asymmetric structure factor by preliminary defined parameters.
#     For more details see documentation module "Structure factor".
#     """
#     # dimension of pr_1 is [hkl, symmetry, a]
#     if len(scat_length_neutron.shape) == 1:
#         scat_length = scat_length_neutron[na, na, :] # neutron diffraction [atoms]
#     elif len(scat_length_neutron.shape) == 2:
#         scat_length = scat_length_neutron[:, na, :] # X-ray diffraction [hkl, atoms]
# 
#     res = (pr_2*(pr_1*atom_multiplicity[na, na, :]*scat_length*atom_occupancy[na, na, :]*debye_waller).sum(axis=2)).sum(axis=1)/pr_2.shape[-1]
#     dder = {}
#     if flag_scat_length_neutron:
#         dder["scat_length_neutron_real"] = (numpy.expand_dims(pr_2, axis=2)*\
#             (pr_1*atom_multiplicity*atom_occupancy*debye_waller)).sum(axis=1)/pr_2.shape[-1] # FIXME: only for neutron diffraction
#         dder["scat_length_neutron_imag"] = 1j*(numpy.expand_dims(pr_2, axis=2)*\
#             (pr_1*atom_multiplicity*atom_occupancy*debye_waller)).sum(axis=1)/pr_2.shape[-1] # FIXME: only for neutron diffraction
#     if flag_debye_waller:
#         dder["debye_waller"] = (numpy.expand_dims(pr_2, axis=2)*\
#             (pr_1*atom_multiplicity*atom_occupancy)*scat_length).sum(axis=1)/pr_2.shape[-1]
#     if flag_atom_occupancy:
#         dder["atom_occupancy"] = (numpy.expand_dims(pr_2, axis=2)*\
#             (pr_1*atom_multiplicity*debye_waller)*scat_length).sum(axis=1)/pr_2.shape[-1]
#     if flag_pr_1:
#         dder["pr_1_real"] = (pr_2[:,:,na]*atom_multiplicity[na, na, :]*atom_occupancy[na, na, :]*scat_length*debye_waller)/pr_2.shape[-1]
#         dder["pr_1_imag"] = 1j*(pr_2[:,:,na]*atom_multiplicity[na, na, :]*atom_occupancy[na, na, :]*scat_length*debye_waller)/pr_2.shape[-1]
#     return res, dder



def calc_f_by_f_asym_a_pr(f_asym_a, scattering_length, pr_3, centrosymmetry, pr_4, flag_f_asym_a: bool = False, flag_scattering_length: bool = False):
    """Calculate structure factor by preliminary defined parameters.
    For more details see documentation module "Structure factor".

    Dimensions:
    f_asym_a = [9, hkl, a] or [hkl, a]
    scattering length = [hkl, a] or [a]
    pr_3 = [hkl]
    """
    if len(scattering_length.shape) == 1:
        scat_length_2d = scattering_length[na, :] # neutron diffraction [atoms]
    elif len(scattering_length.shape) == 2:
        scat_length_2d = scattering_length[:, :] # X-ray diffraction [hkl, atoms]
    
    if len(f_asym_a.shape) == 2: # for structure factor [hkl, a]
        sum_axis = 1
        hh = scat_length_2d
        pr_3_ext = pr_3
    elif len(f_asym_a.shape) == 3: # for tensor structure factor [9, hkl, a]
        sum_axis = 2
        hh = scat_length_2d[na, :, :]
        pr_3_ext = numpy.expand_dims(pr_3, axis=0) # [9, hkl]
        

    f_asym = (hh * f_asym_a).sum(axis=sum_axis)
    f_asym_conj = (hh * f_asym_a.conjugate()).sum(axis=sum_axis)
    f_h = pr_3_ext * f_asym
    f_h_conj = pr_3_ext.conjugate() * f_asym_conj
    if centrosymmetry:
        res= 0.5*(f_h+pr_4*f_h_conj)
    else:
        res= f_h
    dder = {}
    if flag_f_asym_a:
        ofh = numpy.ones(f_h.shape, dtype=float)
        if centrosymmetry:
            hhh_real = pr_3 + pr_4*pr_3.conjugate()
            hhh_imag = pr_3 - pr_4*pr_3.conjugate()
            dder["f_asym_a_real"] = 0.5*(numpy.expand_dims(hhh_real, axis=-1))*hh
            dder["f_asym_a_imag"] = 0.5*1j*(numpy.expand_dims(hhh_imag, axis=-1))*hh
        else:
            dder["f_asym_a_real"] = numpy.expand_dims(pr_3_ext, axis=-1)*hh
            dder["f_asym_a_imag"] = numpy.expand_dims(pr_3_ext, axis=-1)*hh*1j

        # ofh = numpy.ones(f_h.shape, dtype=float)
        # if centrosymmetry:
        #     dder["f_asym_real"] = 0.5*(pr_3+pr_4*pr_3.conjugate())*ofh
        #     dder["f_asym_imag"] = 0.5*(pr_3-1j*pr_4*pr_3.conjugate())*ofh
        # else:
        #     dder["f_asym_real"] = pr_3*ofh
        #     dder["f_asym_imag"] = pr_3*1j*ofh
    if flag_scattering_length:
        pass
    return res, dder

# DELETE iT
# def calc_f_by_f_asym_pr(f_asym, pr_3, centrosymmetry, pr_4, flag_f_asym: bool = False):
#     """Calculate structure factor by preliminary defined parameters.
#     For more details see documentation module "Structure factor".
#     """
#     f_h = pr_3 * f_asym
#     if centrosymmetry:
#         res= 0.5*(f_h+pr_4*f_h.conjugate())
#     else:
#         res= f_h
#     dder = {}
#     if flag_f_asym:
#         ofh = numpy.ones(f_h.shape, dtype=float)
#         if centrosymmetry:
#             dder["f_asym_real"] = 0.5*(pr_3+pr_4*pr_3.conjugate())*ofh
#             dder["f_asym_imag"] = 0.5*(pr_3-1j*pr_4*pr_3.conjugate())*ofh
#         else:
#             dder["f_asym_real"] = pr_3*ofh
#             dder["f_asym_imag"] = pr_3*1j*ofh
#     return res, dder


def calc_sft_ccs_asym_a_by_pr(
        atom_para_multiplicity, debye_waller_factor, atom_para_occupancy,
        atom_para_susceptibility, atom_para_sc_chi,
        pr_1, pr_2, pr_5, theta=None,
        flag_debye_waller: bool = False,
        flag_atom_para_occupancy: bool = False, flag_atom_para_susceptibility: bool = False,
        flag_pr_1: bool = False, flag_pr_5: bool = False):
    """Calculate preliminary asymmetric structure factor tensor by preliminary defined parameters in 10**-12 cm.
    For more details see documentation module "Structure factor".

    The susceptibility parameters are give in mu_B
    """
    mas_constr = (0.2695*atom_para_sc_chi * atom_para_susceptibility[na, :, :]).sum(axis=1)

    hh, dder_hh = calc_m_q_inv_m(pr_5[:, :, na], mas_constr[:, na, :], flag_m=False, flag_q=flag_atom_para_susceptibility)
    if theta is not None:
        hh *= numpy.expand_dims(theta, axis=(0,2))
        if flag_atom_para_susceptibility:
            dder_hh["q"] *= numpy.expand_dims(theta, axis=(0, 1, 3))

    hh_1 = atom_para_multiplicity * atom_para_occupancy
    hh_3 = pr_1*debye_waller_factor*hh_1[na, na, :]
    res = (pr_2[na, :, :, na] * hh_3[na, :, :, :] * hh[:, na, :, :]).sum(axis=2)/pr_2.shape[-1]
    dder = {}
    if flag_atom_para_susceptibility:
        dder_hh_2 = 0.2695*(dder_hh["q"][:,:, na,:, :]* atom_para_sc_chi[na, :, :, na,:]).sum(axis=1)
        dder["atom_para_susceptibility"] =  (pr_2[na, na, :, :, na] * hh_3[na, na, :, :, :] * dder_hh_2[:, :, na, :, :]).sum(axis=3)/pr_2.shape[-1]
    return res, dder

# DELETE IT
# def calc_sft_ccs_asym_by_pr(
#         atom_para_multiplicity, atom_para_form_factor, debye_waller_factor, atom_para_occupancy,
#         atom_para_susceptibility, atom_para_sc_chi,
#         pr_1, pr_2, pr_5,
#         flag_atom_para_form_factor: bool = False, flag_debye_waller: bool = False,
#         flag_atom_para_occupancy: bool = False, flag_atom_para_susceptibility: bool = False,
#         flag_pr_1: bool = False, flag_pr_5: bool = False):
#     """Calculate preliminary asymmetric structure factor tensor by preliminary defined parameters in 10**-12 cm.
#     For more details see documentation module "Structure factor".
# 
#     The susceptibility parameters are give in mu_B
#     """
#     mas_constr = (0.2695*atom_para_sc_chi * atom_para_susceptibility[na, :, :]).sum(axis=1)
# 
#     hh, dder_hh = calc_m_q_inv_m(pr_5[:, :, na], mas_constr[:, na, :], flag_m=False, flag_q=flag_atom_para_susceptibility)
#     hh_1 = atom_para_multiplicity * atom_para_occupancy
#     hh_2 = atom_para_form_factor * hh_1[na, :]
#     hh_3 = pr_1*debye_waller_factor*hh_2[:, na, :]
#     res = (pr_2[na, :, :] * (hh_3[na, :, :, :] * hh[:, na, :, :]).sum(axis=3)).sum(axis=2)/pr_2.shape[-1]
#     dder = {}
#     if flag_atom_para_susceptibility:
#         dder_hh_2 = 0.2695*(dder_hh["q"][:,:, na,:, :]* atom_para_sc_chi[na, :, :, na,:]).sum(axis=1)
#         dder["atom_para_susceptibility"] =  (pr_2[na, na, :, :, na] * hh_3[na, na, :, :, :] * dder_hh_2[:, :, na, :, :]).sum(axis=3)/pr_2.shape[-1]
#     return res, dder


def calc_f_nucl_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False):
    """Calculate nuclear structure factor based on the information given in dictionary.
    Output information is written in the same dictionary. The following keys have to be defined.
    """
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(["atom_fract_xyz", "atom_occupancy",
        "atom_scat_length_neutron", "atom_b_iso", "atom_beta", "unit_cell_parameters"])
    diff_set_crystal = necessary_crystal_keys.difference(set(dict_crystal_keys))
    if len(diff_set_crystal) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_crystal:}")

    flag_reduced_symm_elems = len(set(["reduced_symm_elems", "centrosymmetry", "translation_elems"]).difference(set(dict_crystal_keys))) == 0
    flag_full_symm_elems = len(set(["full_symm_elems", ]).difference(set(dict_crystal_keys))) == 0
    flag_full_mcif_elems = len(set(["full_mcif_elems", ]).difference(set(dict_crystal_keys))) == 0

    if not(flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems):
        raise AttributeError("The symmetry elements have to be defined.")

    necessary_in_out_keys = set(["index_hkl", ])
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_in_out:}")

    index_hkl = dict_in_out["index_hkl"]
    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems[:13]
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]

    atom_fract_xyz = dict_crystal["atom_fract_xyz"]
    atom_site_sc_fract = dict_crystal["atom_site_sc_fract"] 
    atom_site_sc_b = dict_crystal["atom_site_sc_b"] 
    atom_fract_xyz = calc_m_v(atom_site_sc_fract, numpy.mod(atom_fract_xyz, 1), flag_m=False, flag_v=False)[0] + atom_site_sc_b

    atom_occupancy = dict_crystal["atom_occupancy"]
    scat_length_neutron = dict_crystal["atom_scat_length_neutron"]
    atom_b_iso = dict_crystal["atom_b_iso"]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,)+atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta*numpy.expand_dims(atom_beta, axis=0)).sum(axis=1)

    flag_unit_cell_parameters = numpy.any(dict_crystal["flags_unit_cell_parameters"])
    flag_atom_fract_xyz = numpy.any(dict_crystal["flags_atom_fract_xyz"])
    flag_atom_occupancy = numpy.any(dict_crystal["flags_atom_occupancy"])
    flag_atom_b_iso = numpy.any(dict_crystal["flags_atom_b_iso"])
    flag_atom_beta = numpy.any(dict_crystal["flags_atom_beta"])

    
    f_nucl, dder = calc_f_nucl(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_fract_xyz, atom_occupancy, scat_length_neutron, atom_b_iso, atom_beta,
        dict_in_out, 
        flag_unit_cell_parameters=flag_unit_cell_parameters, flag_atom_fract_xyz=flag_atom_fract_xyz,
        flag_atom_occupancy=flag_atom_occupancy, flag_atom_b_iso=flag_atom_b_iso, flag_atom_beta=flag_atom_beta,
        flag_use_precalculated_data=flag_use_precalculated_data)
    
    if "atom_multiplicity" in dict_in_out.keys():
        dict_crystal["atom_multiplicity"] = dict_in_out["atom_multiplicity"]
    return f_nucl, dder


def calc_f_nucl(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_fract_xyz, atom_occupancy, scat_length_neutron, atom_b_iso, atom_beta,
        dict_in_out: dict = None,
        flag_unit_cell_parameters: bool = False, flag_atom_fract_xyz: bool = False,
        flag_atom_occupancy: bool = False, flag_atom_b_iso: bool = False, flag_atom_beta: bool = False,
        flag_use_precalculated_data: bool = False):
    """Calculate nuclear structure factor based on the information given in dictionary.
    Output information is written in the same dictionary. The following keys have to be defined.
    """
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if (flag_use_precalculated_data and ('index_hkl' in dict_in_out_keys)):
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if (flag_use_precalculated_data and ("atom_multiplicity" in dict_in_out_keys)):
        atom_multiplicity = dict_in_out["atom_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(atom_fract_xyz)

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems)
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        atom_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, atom_symm_elems)
        if flag_dict:
            dict_in_out["atom_multiplicity"] = atom_multiplicity

    flag_pr_1 = flag_atom_fract_xyz
    if (flag_use_precalculated_data and ("pr_1" in dict_in_out_keys) and not(flag_atom_fract_xyz)):
        pr_1 = dict_in_out["pr_1"] 
    else:
        pr_1, dder_pr_1 = calc_pr1(index_hkl, reduced_symm_elems, atom_fract_xyz, flag_fract_xyz=flag_atom_fract_xyz)
        if flag_dict:
            dict_in_out["pr_1"] = pr_1

    if (flag_use_precalculated_data and ("pr_2" in dict_in_out_keys)):
        pr_2 = dict_in_out["pr_2"] 
    else:
        pr_2 = calc_pr2(index_hkl, reduced_symm_elems)
        if flag_dict:
            dict_in_out["pr_2"] = pr_2 

    if (flag_use_precalculated_data and ("pr_3" in dict_in_out_keys)):
        pr_3 = dict_in_out["pr_3"] 
    else:
        pr_3 = calc_pr3(index_hkl, translation_elems)
        if flag_dict:
            dict_in_out["pr_3"] = pr_3 

    if (flag_use_precalculated_data and ("pr_4" in dict_in_out_keys)):
        pr_4 = dict_in_out["pr_4"] 
    else:
        pr_4 = calc_pr4(index_hkl, centrosymmetry_position)
        if flag_dict:
            dict_in_out["pr_4"] = pr_4 


    flag_sthovl = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("sthovl" in dict_in_out_keys) and not(flag_sthovl)):
        sthovl = dict_in_out["sthovl"] 
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = flag_sthovl or flag_atom_b_iso or flag_atom_beta
    if (flag_use_precalculated_data and ("debye_waller_factor" in dict_in_out_keys) and not(flag_debye_waller_factor)):
        debye_waller_factor = dict_in_out["debye_waller_factor"] 
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na], sthovl[:, na, na], atom_b_iso[na, na, :],
            atom_beta[:, na, na, :], reduced_symm_elems[:, na, :, na],
            flag_sthovl=flag_sthovl, flag_b_iso=flag_atom_b_iso, flag_beta=flag_atom_beta)
        if flag_dict:
            dict_in_out["debye_waller_factor"] = debye_waller_factor

    flag_scat_length_neutron = False
    flag_debye_waller = flag_atom_b_iso or flag_atom_beta
    flag_f_asym = flag_scat_length_neutron or flag_debye_waller or flag_pr_1 

    if (flag_use_precalculated_data and ("f_asym" in dict_in_out_keys) and
            not(flag_f_asym)):
        f_asym = dict_in_out["f_asym"]
    else: 
        f_asym, dder_f_asym = calc_f_asym_a_by_pr(
            atom_multiplicity, debye_waller_factor, atom_occupancy,
            pr_1, pr_2, 
            flag_debye_waller=flag_debye_waller, flag_atom_occupancy=flag_atom_occupancy,
            flag_pr_1=flag_pr_1)
        if flag_dict:
            dict_in_out["f_asym"] = f_asym

    flag_f_nucl = flag_f_asym
    
    if (flag_use_precalculated_data and ("f_nucl" in dict_in_out_keys) and
            not(flag_f_nucl)):
        f_nucl = dict_in_out["f_nucl"]
    else:
        f_nucl, dder_f_nucl = calc_f_by_f_asym_a_pr(f_asym, scat_length_neutron, pr_3, centrosymmetry, pr_4, flag_f_asym_a=flag_f_asym, flag_scattering_length=flag_scat_length_neutron)
        if flag_dict:
            dict_in_out["f_nucl"] = f_nucl

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = None

    if flag_atom_fract_xyz:
        dder["atom_fract_xyz"] = None

    if flag_atom_occupancy:
        dder["atom_occupancy"] = None

    if flag_atom_b_iso:
        dder["atom_b_iso"] = None

    if flag_atom_beta:
        dder["atom_beta"] = None

    return f_nucl, dder



def calc_f_charge_by_dictionary(dict_crystal, wavelength:float, dict_in_out, flag_use_precalculated_data: bool = False):
    """Calculate nuclear structure factor based on the information given in dictionary.
    Output information is written in the same dictionary. The following keys have to be defined.
    """
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(["atom_fract_xyz", "atom_occupancy",
        "atom_scat_length_neutron", "atom_b_iso", "atom_beta", "unit_cell_parameters"])
    diff_set_crystal = necessary_crystal_keys.difference(set(dict_crystal_keys))
    if len(diff_set_crystal) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_crystal:}")

    flag_reduced_symm_elems = len(set(["reduced_symm_elems", "centrosymmetry", "translation_elems"]).difference(set(dict_crystal_keys))) == 0
    flag_full_symm_elems = len(set(["full_symm_elems", ]).difference(set(dict_crystal_keys))) == 0
    flag_full_mcif_elems = len(set(["full_mcif_elems", ]).difference(set(dict_crystal_keys))) == 0

    if not(flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems):
        raise AttributeError("The symmetry elements have to be defined.")

    necessary_in_out_keys = set(["index_hkl", ])
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_in_out:}")

    index_hkl = dict_in_out["index_hkl"]
    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems[:13]
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]

    atom_fract_xyz = dict_crystal["atom_fract_xyz"]
    atom_site_sc_fract = dict_crystal["atom_site_sc_fract"] 
    atom_site_sc_b = dict_crystal["atom_site_sc_b"] 
    atom_fract_xyz = calc_m_v(atom_site_sc_fract, numpy.mod(atom_fract_xyz, 1), flag_m=False, flag_v=False)[0] + atom_site_sc_b

    atom_occupancy = dict_crystal["atom_occupancy"]
    
    table_sthovl = dict_crystal["table_sthovl"]
    table_atom_scattering_amplitude = dict_crystal["table_atom_scattering_amplitude"]

    table_wavelength = dict_crystal["table_wavelength"]
    table_atom_dispersion = dict_crystal["table_atom_dispersion"]
    atom_dispersion = numpy.array([numpy.interp(float(wavelength), table_wavelength, hh) for hh in table_atom_dispersion], dtype=complex)
    dict_in_out["atom_dispersion"] = atom_dispersion

    atom_b_iso = dict_crystal["atom_b_iso"]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,)+atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta*numpy.expand_dims(atom_beta, axis=0)).sum(axis=1)

    flag_unit_cell_parameters = numpy.any(dict_crystal["flags_unit_cell_parameters"])
    flag_atom_fract_xyz = numpy.any(dict_crystal["flags_atom_fract_xyz"])
    flag_atom_occupancy = numpy.any(dict_crystal["flags_atom_occupancy"])
    flag_atom_b_iso = numpy.any(dict_crystal["flags_atom_b_iso"])
    flag_atom_beta = numpy.any(dict_crystal["flags_atom_beta"])

    
    f_charge, dder = calc_f_charge(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_fract_xyz, atom_occupancy, table_sthovl, table_atom_scattering_amplitude, atom_dispersion, atom_b_iso, atom_beta,
        dict_in_out, 
        flag_unit_cell_parameters=flag_unit_cell_parameters, flag_atom_fract_xyz=flag_atom_fract_xyz,
        flag_atom_occupancy=flag_atom_occupancy, flag_atom_b_iso=flag_atom_b_iso, flag_atom_beta=flag_atom_beta,
        flag_use_precalculated_data=flag_use_precalculated_data)
    
    if "atom_multiplicity" in dict_in_out.keys():
        dict_crystal["atom_multiplicity"] = dict_in_out["atom_multiplicity"]
    return f_charge, dder



def calc_f_charge(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_fract_xyz, atom_occupancy, table_sthovl, table_atom_scattering_amplitude, atom_dispersion, atom_b_iso, atom_beta,
        dict_in_out: dict = None,
        flag_unit_cell_parameters: bool = False, flag_atom_fract_xyz: bool = False,
        flag_atom_occupancy: bool = False, flag_atom_b_iso: bool = False, flag_atom_beta: bool = False,
        flag_use_precalculated_data: bool = False):
    """Calculate nuclear structure factor based on the information given in dictionary.
    Output information is written in the same dictionary. The following keys have to be defined.
    """
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if (flag_use_precalculated_data and ('index_hkl' in dict_in_out_keys)):
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if (flag_use_precalculated_data and ("atom_multiplicity" in dict_in_out_keys)):
        atom_multiplicity = dict_in_out["atom_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(atom_fract_xyz)

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems)
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        atom_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, atom_symm_elems)
        if flag_dict:
            dict_in_out["atom_multiplicity"] = atom_multiplicity

    flag_pr_1 = flag_atom_fract_xyz
    if (flag_use_precalculated_data and ("pr_1" in dict_in_out_keys) and not(flag_atom_fract_xyz)):
        pr_1 = dict_in_out["pr_1"] 
    else:
        pr_1, dder_pr_1 = calc_pr1(index_hkl, reduced_symm_elems, atom_fract_xyz, flag_fract_xyz=flag_atom_fract_xyz)
        if flag_dict:
            dict_in_out["pr_1"] = pr_1

    if (flag_use_precalculated_data and ("pr_2" in dict_in_out_keys)):
        pr_2 = dict_in_out["pr_2"] 
    else:
        pr_2 = calc_pr2(index_hkl, reduced_symm_elems)
        if flag_dict:
            dict_in_out["pr_2"] = pr_2 

    if (flag_use_precalculated_data and ("pr_3" in dict_in_out_keys)):
        pr_3 = dict_in_out["pr_3"] 
    else:
        pr_3 = calc_pr3(index_hkl, translation_elems)
        if flag_dict:
            dict_in_out["pr_3"] = pr_3 

    if (flag_use_precalculated_data and ("pr_4" in dict_in_out_keys)):
        pr_4 = dict_in_out["pr_4"] 
    else:
        pr_4 = calc_pr4(index_hkl, centrosymmetry_position)
        if flag_dict:
            dict_in_out["pr_4"] = pr_4 


    flag_sthovl = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("sthovl" in dict_in_out_keys) and not(flag_sthovl)):
        sthovl = dict_in_out["sthovl"] 
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = flag_sthovl or flag_atom_b_iso or flag_atom_beta
    if (flag_use_precalculated_data and ("debye_waller_factor" in dict_in_out_keys) and not(flag_debye_waller_factor)):
        debye_waller_factor = dict_in_out["debye_waller_factor"] 
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na], sthovl[:, na, na], atom_b_iso[na, na, :],
            atom_beta[:, na, na, :], reduced_symm_elems[:, na, :, na],
            flag_sthovl=flag_sthovl, flag_b_iso=flag_atom_b_iso, flag_beta=flag_atom_beta)
        if flag_dict:
            dict_in_out["debye_waller_factor"] = debye_waller_factor

    flag_scat_length_neutron = False
    flag_debye_waller = flag_atom_b_iso or flag_atom_beta
    flag_f_asym = flag_scat_length_neutron or flag_debye_waller or flag_pr_1 

    if (flag_use_precalculated_data and ("f_asym" in dict_in_out_keys) and
            not(flag_f_asym)):
        f_asym = dict_in_out["f_asym"]
    else: 
        l_scat_length_xray = [
            numpy.interp(sthovl, table_sthovl, table_sc_ampl) for table_sc_ampl in table_atom_scattering_amplitude]
        hh = numpy.stack(l_scat_length_xray, axis=1)
        scat_length_xray = (
            numpy.stack(l_scat_length_xray, axis=1) + 
            numpy.expand_dims(atom_dispersion, axis=0)
        )

        f_asym, dder_f_asym = calc_f_asym_a_by_pr(
            atom_multiplicity, debye_waller_factor, atom_occupancy,
            pr_1, pr_2, 
            flag_debye_waller=flag_debye_waller, flag_atom_occupancy=flag_atom_occupancy,
            flag_pr_1=flag_pr_1)
        if flag_dict:
            dict_in_out["f_asym"] = f_asym
    
    flag_f_charge = flag_f_asym
    if (flag_use_precalculated_data and ("f_charge" in dict_in_out_keys) and
            not(flag_f_charge)):
        f_charge = dict_in_out["f_charge"]
    else:
        f_charge, dder_f_charge = calc_f_by_f_asym_a_pr(f_asym, scat_length_xray, pr_3, centrosymmetry, pr_4, flag_f_asym_a=flag_f_asym, flag_scattering_length=flag_scat_length_neutron)
        if flag_dict:
            dict_in_out["f_charge"] = f_charge

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = None

    if flag_atom_fract_xyz:
        dder["atom_fract_xyz"] = None

    if flag_atom_occupancy:
        dder["atom_occupancy"] = None

    if flag_atom_b_iso:
        dder["atom_b_iso"] = None

    if flag_atom_beta:
        dder["atom_beta"] = None

    return f_charge, dder


def calc_sft_ccs_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False):
    """Calculate structure factor tensor in CCS (X||a*, Z||c) based on the information given in dictionary.
    Output information is written in the same dictionary. 
    """
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(["unit_cell_parameters", ])
    diff_set_crystal = necessary_crystal_keys.difference(set(dict_crystal_keys))
    if len(diff_set_crystal) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_crystal:}")

    flag_reduced_symm_elems = len(set(["reduced_symm_elems", "centrosymmetry", "translation_elems"]).difference(set(dict_crystal_keys))) == 0
    flag_full_symm_elems = len(set(["full_symm_elems", ]).difference(set(dict_crystal_keys))) == 0
    flag_full_mcif_elems = len(set(["full_mcif_elems", ]).difference(set(dict_crystal_keys))) == 0

    if not(flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems):
        raise AttributeError("The symmetry elements have to be defined.")

    necessary_in_out_keys = set(["index_hkl", ])
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_in_out:}")

    index_hkl = dict_in_out["index_hkl"]

    non_zero_keys = set(["mag_atom_lande_factor", "mag_atom_kappa",
        "mag_atom_j0_parameters", "mag_atom_j2_parameters"])
    diff_set_crystal = non_zero_keys.difference(set(dict_crystal_keys))
    
    if len(diff_set_crystal) != 0:
        sft_ccs = numpy.zeros((9, index_hkl.shape[-1]), dtype=complex)
        dder = {}
        return sft_ccs, dder

    if "flag_only_orbital" in dict_in_out_keys:
        flag_only_orbital = dict_in_out["flag_only_orbital"]
    else:
        flag_only_orbital = False

    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems 
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_para_index = dict_crystal["atom_para_index"]
    atom_para_fract_xyz = dict_crystal["atom_fract_xyz"][:, atom_para_index]
    atom_para_sc_fract = dict_crystal["atom_site_sc_fract"][:, atom_para_index]
    atom_para_sc_b = dict_crystal["atom_site_sc_b"][:, atom_para_index]
    atom_para_fract_xyz = calc_m_v(
         atom_para_sc_fract, numpy.mod(atom_para_fract_xyz, 1), flag_m=False, flag_v=False)[0] + atom_para_sc_b

    atom_para_occupancy = dict_crystal["atom_occupancy"][atom_para_index]
    atom_para_b_iso = dict_crystal["atom_b_iso"][atom_para_index]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,)+atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta*numpy.expand_dims(atom_beta, axis=0)).sum(axis=1)
    atom_para_beta = atom_beta[:, atom_para_index]

    mag_atom_para_index = dict_crystal["mag_atom_para_index"]
    atom_para_lande_factor = dict_crystal["mag_atom_lande_factor"][mag_atom_para_index]
    atom_para_kappa = dict_crystal["mag_atom_kappa"][mag_atom_para_index]
    atom_para_j0_parameters = dict_crystal["mag_atom_j0_parameters"][:, mag_atom_para_index]
    atom_para_j2_parameters = dict_crystal["mag_atom_j2_parameters"][:, mag_atom_para_index]

    atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
    atom_para_sc_chi = dict_crystal["atom_para_sc_chi"] 
    flag_unit_cell_parameters = numpy.any(dict_crystal["flags_unit_cell_parameters"])
    flag_atom_para_fract_xyz = numpy.any(dict_crystal["flags_atom_fract_xyz"][:, atom_para_index])
    flag_atom_para_occupancy = numpy.any(dict_crystal["flags_atom_occupancy"][atom_para_index])
    flag_atom_para_b_iso = numpy.any(dict_crystal["flags_atom_b_iso"][atom_para_index])
    flag_atom_para_beta = numpy.any(dict_crystal["flags_atom_beta"][:, atom_para_index])
    flag_atom_para_susceptibility = numpy.any(dict_crystal["flags_atom_para_susceptibility"])
    flag_atom_para_lande_factor = numpy.any(dict_crystal["flags_mag_atom_lande_factor"][mag_atom_para_index])
    flag_atom_para_kappa = numpy.any(dict_crystal["flags_mag_atom_kappa"][mag_atom_para_index])

    
    sft_ccs, dder = calc_sft_ccs(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_para_fract_xyz, atom_para_occupancy, atom_para_susceptibility, atom_para_b_iso, atom_para_beta,
        atom_para_lande_factor, atom_para_kappa, atom_para_j0_parameters, atom_para_j2_parameters, atom_para_sc_chi,
        dict_in_out=dict_in_out, flag_only_orbital=flag_only_orbital,
        flag_unit_cell_parameters=flag_unit_cell_parameters, flag_atom_para_fract_xyz=flag_atom_para_fract_xyz,
        flag_atom_para_occupancy=flag_atom_para_occupancy, flag_atom_para_susceptibility=flag_atom_para_susceptibility,
        flag_atom_para_b_iso=flag_atom_para_b_iso, flag_atom_para_beta=flag_atom_para_beta,
        flag_atom_para_lande_factor=flag_atom_para_lande_factor, flag_atom_para_kappa=flag_atom_para_kappa, 
        flag_use_precalculated_data=flag_use_precalculated_data)
    return sft_ccs, dder


def calc_sft_ccs(index_hkl,
        reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_para_fract_xyz, atom_para_occupancy, atom_para_susceptibility, atom_para_b_iso, atom_para_beta,
        atom_para_lande_factor, atom_para_kappa, atom_para_j0_parameters, atom_para_j2_parameters, atom_para_sc_chi,
        dict_in_out: dict = None, flag_only_orbital: bool = False,
        flag_unit_cell_parameters: bool = False, flag_atom_para_fract_xyz: bool = False,
        flag_atom_para_occupancy: bool = False, flag_atom_para_susceptibility: bool = False,
        flag_atom_para_b_iso: bool = False, flag_atom_para_beta: bool = False,
        flag_atom_para_lande_factor: bool = False, flag_atom_para_kappa: bool = False, 
        flag_use_precalculated_data: bool = False):
    """Calculate structure factor tensor in Cartesian coordinate system with X||a*, Z||c in 10**-12 cm.
    Note, that the susceptibility parameters are given in mu_B. 
    """
    
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if 'index_hkl' in dict_in_out_keys:
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if (flag_use_precalculated_data and ("atom_para_multiplicity" in dict_in_out_keys)):
        mag_atom_multiplicity = dict_in_out["atom_para_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(atom_para_fract_xyz)

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems[:13], centrosymmetry, centrosymmetry_position, translation_elems)
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        mag_atom_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, atom_symm_elems)
        if flag_dict:
            dict_in_out["atom_para_multiplicity"] = mag_atom_multiplicity

    flag_pr_1 = flag_atom_para_fract_xyz
    if (flag_use_precalculated_data and ("pr_1_atom_para" in dict_in_out_keys) and not(flag_atom_para_fract_xyz)):
        pr_1 = dict_in_out["pr_1_atom_para"] 
    else:
        pr_1, dder_pr_1 = calc_pr1(index_hkl, reduced_symm_elems[:13], atom_para_fract_xyz, flag_fract_xyz=flag_atom_para_fract_xyz)
        if flag_dict:
            dict_in_out["pr_1_atom_para"] = pr_1

    if (flag_use_precalculated_data and ("pr_2" in dict_in_out_keys)):
        pr_2 = dict_in_out["pr_2"] 
    else:
        pr_2 = calc_pr2(index_hkl, reduced_symm_elems[:13])
        if flag_dict:
            dict_in_out["pr_2"] = pr_2 

    if (flag_use_precalculated_data and ("pr_3" in dict_in_out_keys)):
        pr_3 = dict_in_out["pr_3"] 
    else:
        pr_3 = calc_pr3(index_hkl, translation_elems)
        if flag_dict:
            dict_in_out["pr_3"] = pr_3 

    if (flag_use_precalculated_data and ("pr_4" in dict_in_out_keys)):
        pr_4 = dict_in_out["pr_4"] 
    else:
        pr_4 = calc_pr4(index_hkl, centrosymmetry_position)
        if flag_dict:
            dict_in_out["pr_4"] = pr_4 

    flag_sthovl = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("sthovl" in dict_in_out_keys) and not(flag_sthovl)):
        sthovl = dict_in_out["sthovl"] 
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    flag_pr_5 = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("pr_5" in dict_in_out_keys) and not(flag_pr_5)):
        pr_5 = dict_in_out["pr_5"] 
    else:
        pr_5, dder_pr_5 = calc_pr5(reduced_symm_elems[:13], unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["pr_5"] = pr_5

    flag_atom_para_form_factor = (flag_sthovl or flag_atom_para_lande_factor or flag_atom_para_kappa)
    flag_hh = True
    if "flag_only_orbital" in dict_in_out_keys:
        flag_hh = flag_only_orbital == dict_in_out["flag_only_orbital"]
    dict_in_out["flag_only_orbital"] = flag_only_orbital
    if (flag_use_precalculated_data and ("atom_para_form_factor" in dict_in_out_keys) and not(flag_atom_para_form_factor) and flag_hh):
        atom_para_form_factor = dict_in_out["atom_para_form_factor"] 
    else:
        atom_para_form_factor, dder_ff = calc_form_factor(
            sthovl[:, na], atom_para_lande_factor[na, :], atom_para_kappa[na, :], atom_para_j0_parameters[:, na, :], atom_para_j2_parameters[:, na, :],
            flag_lande_factor=flag_atom_para_lande_factor,
            flag_only_orbital=flag_only_orbital,
            flag_sthovl=flag_sthovl, 
            flag_kappa=flag_atom_para_kappa)
        if flag_dict:
            dict_in_out["atom_para_form_factor"] = atom_para_form_factor

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = flag_sthovl or flag_atom_para_b_iso or flag_atom_para_beta
    if (flag_use_precalculated_data and ("atom_para_debye_waller_factor" in dict_in_out_keys) and not(flag_debye_waller_factor)):
        debye_waller_factor = dict_in_out["atom_para_debye_waller_factor"] 
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na], sthovl[:, na, na], atom_para_b_iso[na, na, :],
            atom_para_beta[:, na, na, :], reduced_symm_elems[:13, na, :, na],
            flag_sthovl=flag_sthovl, flag_b_iso=flag_atom_para_b_iso, flag_beta=flag_atom_para_beta)
        if flag_dict:
            dict_in_out["atom_para_debye_waller_factor"] = debye_waller_factor

    flag_scat_length_neutron = False
    flag_debye_waller = flag_atom_para_b_iso or flag_atom_para_beta

    flag_sft_ccs_asym = flag_atom_para_form_factor or flag_debye_waller or flag_atom_para_occupancy or flag_atom_para_susceptibility or flag_pr_1 or flag_pr_5
    if (flag_use_precalculated_data and ("sft_ccs_asym" in dict_in_out_keys) and
            not(flag_sft_ccs_asym)):
        sft_ccs_asym = dict_in_out["sft_ccs_asym"]
    else: 
        theta = None
        # if reduced_symm_elems.shape[0] == 14:
        #     theta = reduced_symm_elems[13] # * calc_det_m(reduced_symm_elems[4:13], flag_m=False)[0]
        #     # print("theta: ", theta)
        sft_ccs_asym, dder_sft_ccs_asym = calc_sft_ccs_asym_a_by_pr(
            mag_atom_multiplicity, debye_waller_factor, atom_para_occupancy, atom_para_susceptibility, atom_para_sc_chi,
            pr_1, pr_2, pr_5, theta=theta,
            flag_debye_waller=flag_debye_waller, flag_atom_para_occupancy=flag_atom_para_occupancy,
            flag_atom_para_susceptibility = flag_atom_para_susceptibility,
            flag_pr_1=flag_pr_1, flag_pr_5=flag_pr_5)
        if flag_dict:
            dict_in_out["sft_ccs_asym"] = sft_ccs_asym

    flag_sft_ccs = flag_sft_ccs_asym
    if (flag_use_precalculated_data and ("sft_ccs" in dict_in_out_keys) and
            not(flag_sft_ccs)):
        sft_ccs = dict_in_out["sft_ccs"]
    else:
        sft_ccs, dder_sft_ccs = calc_f_by_f_asym_a_pr(sft_ccs_asym, atom_para_form_factor, pr_3, centrosymmetry, pr_4, flag_f_asym_a=flag_sft_ccs_asym, flag_scattering_length=flag_atom_para_form_factor)
        if flag_dict:
            dict_in_out["sft_ccs"] = sft_ccs

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = None

    if flag_atom_para_fract_xyz:
        dder["atom_para_fract_xyz"] = None

    if flag_atom_para_occupancy:
        dder["atom_para_occupancy"] = None

    if flag_atom_para_b_iso:
        dder["atom_para_b_iso"] = None

    if flag_atom_para_beta:
        dder["atom_para_beta"] = None
    if flag_atom_para_susceptibility:
        dder["atom_para_susceptibility"] = (
            dder_sft_ccs["f_asym_a_real"][:, na, :, :]*dder_sft_ccs_asym["atom_para_susceptibility"]+
            dder_sft_ccs["f_asym_a_imag"][:, na, :, :]*dder_sft_ccs_asym["atom_para_susceptibility"])
    return sft_ccs, dder


def calc_index_hkl_multiplicity_in_range(sthovl_min, sthovl_max, unit_cell_parameters, reduced_symm_elems, translation_elems, centrosymmetry: bool):
    a, b, c = unit_cell_parameters[0], unit_cell_parameters[1], unit_cell_parameters[2]
    h_max = int(2.*a*sthovl_max)
    k_max = int(2.*b*sthovl_max)
    l_max = int(2.*c*sthovl_max)
    
    index_h = numpy.arange(-h_max, h_max+1, 1, dtype=int)
    index_k = numpy.arange(-k_max, k_max+1, 1, dtype=int)
    index_l = numpy.arange(-l_max, l_max+1, 1, dtype=int)

    index_h, index_k, index_l = numpy.meshgrid(index_h, index_k, index_l, indexing="ij")
    index_h, index_k, index_l = index_h.flatten(), index_k.flatten(), index_l.flatten()
    index_hkl_full = numpy.stack([index_h, index_k, index_l], axis=0)
    index_hkl_equivalent = calc_equivalent_reflections(index_hkl_full, reduced_symm_elems, centrosymmetry=centrosymmetry)
    
    label_hkl_equivalent = 1000000*index_hkl_equivalent[0] + 1000*index_hkl_equivalent[1] + index_hkl_equivalent[2]
    index_max = numpy.argsort(label_hkl_equivalent, axis=1)[:,-1]
    index_hkl_sort = index_hkl_equivalent[:, numpy.arange(index_max.size),index_max]
    index_hkl_unique, counts_unique = numpy.unique(index_hkl_sort, axis=1, return_counts=True)

    pr_3 = calc_pr3(index_hkl_unique, translation_elems)
    flag = numpy.logical_not(numpy.isclose(pr_3, 0.))

    index_hkl = index_hkl_unique[:, flag]
    counts = counts_unique[flag]
    sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)

    arg_sort_sthovl = numpy.argsort(sthovl)
    index_hkl_sort = index_hkl[:, arg_sort_sthovl]
    counts_sort = counts[arg_sort_sthovl]
    sthovl_sort = sthovl[arg_sort_sthovl]
    
    flag = numpy.logical_and(sthovl_sort>= sthovl_min, sthovl_sort <= sthovl_max)
    index_hkl_out = index_hkl_sort[:, flag]
    counts_out = counts_sort[flag]
    return index_hkl_out, counts_out


def calc_f_m_perp_ordered_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False):
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(["unit_cell_parameters", "full_mcif_elems"])
    diff_set_crystal = necessary_crystal_keys.difference(set(dict_crystal_keys))
    if len(diff_set_crystal) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_crystal:}")

    necessary_in_out_keys = set(["index_hkl", ])
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_in_out:}")

    index_hkl = dict_in_out["index_hkl"]

    non_zero_keys = set(["mag_atom_lande_factor", "mag_atom_kappa",
        "mag_atom_j0_parameters", "mag_atom_j2_parameters"])
    diff_set_crystal = non_zero_keys.difference(set(dict_crystal_keys))
    
    if len(diff_set_crystal) != 0:
        sft_ccs = numpy.zeros((3, index_hkl.shape[-1]), dtype=complex)
        dder = {}
        return sft_ccs, dder

    if "flag_only_orbital" in dict_in_out_keys:
        flag_only_orbital = dict_in_out["flag_only_orbital"]
    else:
        flag_only_orbital = False

    full_mcif_elems = dict_crystal["full_mcif_elems"]

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_ordered_index = dict_crystal["atom_ordered_index"]
    atom_ordered_fract_xyz = dict_crystal["atom_fract_xyz"][:, atom_ordered_index]
    atom_ordered_sc_fract = dict_crystal["atom_site_sc_fract"][:, atom_ordered_index]
    atom_ordered_sc_b = dict_crystal["atom_site_sc_b"][:, atom_ordered_index]
    atom_ordered_fract_xyz = calc_m_v(
        atom_ordered_sc_fract, numpy.mod(atom_ordered_fract_xyz, 1), flag_m=False, flag_v=False)[0] + atom_ordered_sc_b


    atom_ordered_occupancy = dict_crystal["atom_occupancy"][atom_ordered_index]
    atom_ordered_b_iso = dict_crystal["atom_b_iso"][atom_ordered_index]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,)+atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta*numpy.expand_dims(atom_beta, axis=0)).sum(axis=1)
    atom_ordered_beta = atom_beta[:, atom_ordered_index]


    mag_atom_ordered_index = dict_crystal["mag_atom_ordered_index"]
    atom_ordered_lande_factor = dict_crystal["mag_atom_lande_factor"][mag_atom_ordered_index]
    atom_ordered_kappa = dict_crystal["mag_atom_kappa"][mag_atom_ordered_index]
    atom_ordered_j0_parameters = dict_crystal["mag_atom_j0_parameters"][:, mag_atom_ordered_index]
    atom_ordered_j2_parameters = dict_crystal["mag_atom_j2_parameters"][:, mag_atom_ordered_index]

    atom_ordered_moment_crystalaxis_xyz = dict_crystal["atom_ordered_moment_crystalaxis_xyz"]
    flags_atom_ordered_moment_crystalaxis_xyz = dict_crystal["flags_atom_ordered_moment_crystalaxis_xyz"]
    flag_unit_cell_parameters = numpy.any(dict_crystal["flags_unit_cell_parameters"])
    flag_atom_ordered_fract_xyz = numpy.any(dict_crystal["flags_atom_fract_xyz"][:, atom_ordered_index])
    flag_atom_ordered_occupancy = numpy.any(dict_crystal["flags_atom_occupancy"][atom_ordered_index])
    flag_atom_ordered_b_iso = numpy.any(dict_crystal["flags_atom_b_iso"][atom_ordered_index])
    flag_atom_ordered_beta = numpy.any(dict_crystal["flags_atom_beta"][:, atom_ordered_index])
    flag_atom_ordered_lande_factor = numpy.any(dict_crystal["flags_mag_atom_lande_factor"][mag_atom_ordered_index])
    flag_atom_ordered_kappa = numpy.any(dict_crystal["flags_mag_atom_kappa"][mag_atom_ordered_index])

    

    flag_atom_ordered_moment_crystalaxis_xyz = numpy.any(flags_atom_ordered_moment_crystalaxis_xyz)
    f_m_perp_o, dder = calc_f_m_perp_ordered(index_hkl,
        full_mcif_elems,
        unit_cell_parameters, atom_ordered_fract_xyz, atom_ordered_occupancy, atom_ordered_moment_crystalaxis_xyz, atom_ordered_b_iso, atom_ordered_beta,
        atom_ordered_lande_factor, atom_ordered_kappa, atom_ordered_j0_parameters, atom_ordered_j2_parameters, 
        dict_in_out=dict_in_out, flag_only_orbital=flag_only_orbital,
        flag_unit_cell_parameters=flag_unit_cell_parameters, flag_atom_ordered_fract_xyz=flag_atom_ordered_fract_xyz,
        flag_atom_ordered_occupancy=flag_atom_ordered_occupancy, flag_atom_ordered_moment_crystalaxis_xyz=flag_atom_ordered_moment_crystalaxis_xyz,
        flag_atom_ordered_b_iso=flag_atom_ordered_b_iso, flag_atom_ordered_beta=flag_atom_ordered_beta,
        flag_atom_ordered_lande_factor=flag_atom_ordered_lande_factor, flag_atom_ordered_kappa=flag_atom_ordered_kappa, 
        flag_use_precalculated_data=flag_use_precalculated_data)
    return f_m_perp_o, dder


def calc_f_m_perp_ordered(index_hkl,
        full_mcif_elems,
        unit_cell_parameters, atom_ordered_fract_xyz, atom_ordered_occupancy, atom_ordered_moment_crystalaxis_xyz, atom_ordered_b_iso, atom_ordered_beta,
        atom_ordered_lande_factor, atom_ordered_kappa, atom_ordered_j0_parameters, atom_ordered_j2_parameters, 
        dict_in_out: dict = None, flag_only_orbital: bool = False,
        flag_unit_cell_parameters: bool = False, flag_atom_ordered_fract_xyz: bool = False,
        flag_atom_ordered_occupancy: bool = False, flag_atom_ordered_moment_crystalaxis_xyz: bool = False,
        flag_atom_ordered_b_iso: bool = False, flag_atom_ordered_beta: bool = False,
        flag_atom_ordered_lande_factor: bool = False, flag_atom_ordered_kappa: bool = False, 
        flag_use_precalculated_data: bool = False):

    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if 'index_hkl' in dict_in_out_keys:
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if (flag_use_precalculated_data and ("atom_ordered_multiplicity" in dict_in_out_keys)):
        atom_ordered_multiplicity = dict_in_out["atom_ordered_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(atom_ordered_fract_xyz)

        atom_ordered_multiplicity = calc_multiplicity_by_atom_symm_elems(full_mcif_elems[:13], atom_symm_elems)
        if flag_dict:
            dict_in_out["atom_ordered_multiplicity"] = atom_ordered_multiplicity

    flag_pr_1 = flag_atom_ordered_fract_xyz
    if (flag_use_precalculated_data and ("pr_1_atom_ordered" in dict_in_out_keys) and not(flag_atom_ordered_fract_xyz)):
        pr_1 = dict_in_out["pr_1_atom_ordered"] 
    else:
        pr_1, dder_pr_1 = calc_pr1(index_hkl, full_mcif_elems[:13], atom_ordered_fract_xyz, flag_fract_xyz=flag_atom_ordered_fract_xyz)
        if flag_dict:
            dict_in_out["pr_1_atom_ordered"] = pr_1

    if (flag_use_precalculated_data and ("pr_2" in dict_in_out_keys)):
        pr_2 = dict_in_out["pr_2"] 
    else:
        pr_2 = calc_pr2(index_hkl, full_mcif_elems[:13])
        if flag_dict:
            dict_in_out["pr_2"] = pr_2 


    flag_sthovl = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("sthovl" in dict_in_out_keys) and not(flag_sthovl)):
        sthovl = dict_in_out["sthovl"] 
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    flag_mag_atom_form_factor = (flag_sthovl or flag_atom_ordered_lande_factor or flag_atom_ordered_kappa)
    flag_hh = True
    if "flag_only_orbital" in dict_in_out_keys:
        flag_hh = flag_only_orbital == dict_in_out["flag_only_orbital"]
    dict_in_out["flag_only_orbital"] = flag_only_orbital

    flag_atom_ordered_form_factor = flag_atom_ordered_lande_factor or flag_sthovl or flag_atom_ordered_kappa
    if (flag_use_precalculated_data and ("atom_ordered_form_factor" in dict_in_out_keys) and not(flag_atom_ordered_form_factor) and flag_hh):
        atom_ordered_form_factor = dict_in_out["atom_ordered_form_factor"] 
    else:
        atom_ordered_form_factor, dder_ff = calc_form_factor(
            sthovl[:, na], atom_ordered_lande_factor[na, :], atom_ordered_kappa[na, :], atom_ordered_j0_parameters[:, na, :], atom_ordered_j2_parameters[:, na, :],
            flag_lande_factor=flag_atom_ordered_lande_factor,
            flag_only_orbital=flag_only_orbital,
            flag_sthovl=flag_sthovl, 
            flag_kappa=flag_atom_ordered_kappa)
        if flag_dict:
            dict_in_out["atom_ordered_form_factor"] = atom_ordered_form_factor

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = flag_sthovl or flag_atom_ordered_b_iso or flag_atom_ordered_beta
    if (flag_use_precalculated_data and ("atom_ordered_debye_waller_factor" in dict_in_out_keys) and not(flag_debye_waller_factor)):
        debye_waller_factor = dict_in_out["atom_ordered_debye_waller_factor"] 
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na], sthovl[:, na, na], atom_ordered_b_iso[na, na, :],
            atom_ordered_beta[:, na, na, :], full_mcif_elems[:13, na, :, na],
            flag_sthovl=flag_sthovl, flag_b_iso=flag_atom_ordered_b_iso, flag_beta=flag_atom_ordered_beta)
        if flag_dict:
            dict_in_out["atom_ordered_debye_waller_factor"] = debye_waller_factor

    flag_debye_waller = flag_atom_ordered_b_iso or flag_atom_ordered_beta

    m_norm, der_m_norm = calc_m_m_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    r_direct = full_mcif_elems[4:13, :] 
    rm, der_rm = calc_m_v(
        r_direct[:, :, na], atom_ordered_moment_crystalaxis_xyz[:, na, :],
        flag_m=False, flag_v=flag_atom_ordered_moment_crystalaxis_xyz)
    rm_ccs, der_rm_ccs = calc_m_v(
        m_norm, rm,
        flag_m=flag_unit_cell_parameters, flag_v=flag_atom_ordered_moment_crystalaxis_xyz)
    det_r, der_det_r = calc_det_m(r_direct)
    theta_s = full_mcif_elems[13,:] 

    moment_ccs = 0.2695*rm_ccs*(theta_s*det_r)[na, :, na]
    
    hh_1 = atom_ordered_multiplicity*atom_ordered_occupancy
    hh_2 = atom_ordered_form_factor*hh_1[na, :]
    hh_3 = pr_1*debye_waller_factor*hh_2[:, na, :]
    f_m = (pr_2[na, :, :] * (hh_3[na, :, :, :] * moment_ccs[:, na, :, :]).sum(axis=3)).sum(axis=2)/pr_2.shape[-1]
    eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
    dict_in_out["eq_ccs"] = eq_ccs

    dict_in_out["f_m_o"] = f_m

    flag_f_m = False
    f_m_perp_o, dder_f_m_perp_o = calc_vector_product_v1_v2_v1(eq_ccs, f_m, flag_v1=flag_unit_cell_parameters, flag_v2=flag_f_m)
    dict_in_out["f_m_perp_o"] = f_m_perp_o
    dder = {}
    if flag_atom_ordered_moment_crystalaxis_xyz:
        dder["atom_ordered_moment_crystalaxis_xyz"] = None

    dder = {}
    return f_m_perp_o, dder


def calc_bulk_susceptibility_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False):
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()

    necessary_crystal_keys = set(["unit_cell_parameters", ])
    diff_set_crystal = necessary_crystal_keys.difference(set(dict_crystal_keys))
    if len(diff_set_crystal) != 0:
        raise AttributeError(f"The following attributes have to be defined {diff_set_crystal:}")

    flag_reduced_symm_elems = len(set(["reduced_symm_elems", "centrosymmetry", "translation_elems"]).difference(set(dict_crystal_keys))) == 0
    flag_full_symm_elems = len(set(["full_symm_elems", ]).difference(set(dict_crystal_keys))) == 0
    flag_full_mcif_elems = len(set(["full_mcif_elems", ]).difference(set(dict_crystal_keys))) == 0

    if not(flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems):
        raise AttributeError("The symmetry elements have to be defined.")

    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems[:13]
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_para_index = dict_crystal["atom_para_index"]
    atom_para_fract_xyz = dict_crystal["atom_fract_xyz"][:, atom_para_index]
    atom_para_sc_fract = dict_crystal["atom_site_sc_fract"][:, atom_para_index]
    atom_para_sc_b = dict_crystal["atom_site_sc_b"][:, atom_para_index]
    atom_para_fract_xyz = calc_m_v(
         atom_para_sc_fract, numpy.mod(atom_para_fract_xyz, 1), flag_m=False, flag_v=False)[0] + atom_para_sc_b
    atom_para_occupancy = dict_crystal["atom_occupancy"][atom_para_index]

    atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
    atom_para_sc_chi = dict_crystal["atom_para_sc_chi"] #FIXME not sure that it should be in crystal dictionary but not dict_in_out
    flag_unit_cell_parameters = numpy.any(dict_crystal["flags_unit_cell_parameters"])
    flag_atom_para_fract_xyz = numpy.any(dict_crystal["flags_atom_fract_xyz"][:, atom_para_index])
    flag_atom_para_occupancy = numpy.any(dict_crystal["flags_atom_occupancy"][atom_para_index])
    flag_atom_para_susceptibility = numpy.any(dict_crystal["flags_atom_para_susceptibility"])

    bulk_susceptibility, dder = calc_bulk_susceptibility(reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_para_fract_xyz, atom_para_occupancy, atom_para_susceptibility, atom_para_sc_chi,
        dict_in_out = dict_in_out, flag_unit_cell_parameters = flag_unit_cell_parameters, flag_atom_para_fract_xyz = flag_atom_para_fract_xyz,
        flag_atom_para_occupancy = flag_atom_para_occupancy, flag_atom_para_susceptibility = flag_atom_para_susceptibility,
        flag_use_precalculated_data = flag_use_precalculated_data)
    error_bars = numpy.zeros_like(bulk_susceptibility)
    if 'atom_para_susceptibility' in dder.keys():
        dder_aps = dder['atom_para_susceptibility']
        sigma_aps = dict_crystal["atom_para_susceptibility_sigma"]
        error_bars = numpy.sum(dder_aps*numpy.expand_dims(sigma_aps, axis=0), axis=1)


    return bulk_susceptibility, error_bars


def calc_bulk_susceptibility(reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems,
        unit_cell_parameters, atom_para_fract_xyz, atom_para_occupancy, atom_para_susceptibility, atom_para_sc_chi,
        dict_in_out: dict = None, flag_unit_cell_parameters: bool = False, flag_atom_para_fract_xyz: bool = False,
        flag_atom_para_occupancy: bool = False, flag_atom_para_susceptibility: bool = False, 
        flag_use_precalculated_data: bool = False):

    if dict_in_out is None:
        flag_dict = True
        dict_in_out_keys = []
    else:       
        dict_in_out_keys = dict_in_out.keys()
        flag_dict = True

    if (flag_use_precalculated_data and ("atom_para_multiplicity" in dict_in_out_keys)):
        atom_para_multiplicity = dict_in_out["atom_para_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(atom_para_fract_xyz)

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems)
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        atom_para_multiplicity = calc_multiplicity_by_atom_symm_elems(full_symm_elems, atom_symm_elems)
        if flag_dict:
            dict_in_out["atom_para_multiplicity"] = atom_para_multiplicity

    flag_pr_5 = flag_unit_cell_parameters
    if (flag_use_precalculated_data and ("pr_5" in dict_in_out_keys) and not(flag_pr_5)):
        pr_5 = dict_in_out["pr_5"] 
    else:
        pr_5, dder_pr_5 = calc_pr5(reduced_symm_elems, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
        if flag_dict:
            dict_in_out["pr_5"] = pr_5


    mas_constr = (atom_para_sc_chi * atom_para_susceptibility[na, :, :]).sum(axis=1)

    hh, dder_hh = calc_m_q_inv_m(pr_5[:, :, na], mas_constr[:, na, :], flag_m=False, flag_q=flag_atom_para_susceptibility)

    chi_bulk = (((atom_para_multiplicity * atom_para_occupancy)[na, na, :]*hh).sum(axis=2)).sum(axis=1)/pr_5.shape[-1]
    dder = {}
    return chi_bulk, dder