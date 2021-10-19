"""Calculation for MagCrystal.

Functions
---------
    - calc_b_iso_beta
    - calc_f_nucl
    - calc_f_mag
"""
import numpy

from cryspy.A_functions_base.function_3_mcif import \
    calc_full_sym_elems, calc_multiplicity, calc_moment_by_sym_elem

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_phase_by_hkl_xyz_rb, calc_dwf

from cryspy.A_functions_base.symmetry_elements import calc_full_mag_elems

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAnisoL

from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import \
    SpaceGroupSymopMagnCenteringL

from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import \
    SpaceGroupSymopMagnOperationL

from cryspy.C_item_loop_classes.cl_1_atom_site_moment import AtomSiteMomentL

from cryspy.C_item_loop_classes.cl_2_atom_site_scat import AtomSiteScatL


def calc_b_iso_beta(cell: Cell, atom_site: AtomSiteL,
                    atom_site_aniso: AtomSiteAnisoL):
    """Calculate b_iso and beta_ij based on atom_site and atom_sites.

    For each atom defined in atom_site.
    """
    a_s = atom_site
    a_s_a = atom_site_aniso
    l_b_iso, l_beta = [], []
    coeff = float(8.*numpy.pi**2)
    for item_a_s in a_s.items:
        label_atom = item_a_s.label
        try:
            adp_type = item_a_s.adp_type
        except AttributeError:
            adp_type = None
        b_iso = 0.
        beta = (0., 0., 0., 0., 0., 0.)
        if adp_type == "Uiso":
            u_iso = float(item_a_s.u_iso_or_equiv)
            b_iso = float(8.*numpy.pi**2*u_iso)
        elif adp_type == "Biso":
            b_iso = float(item_a_s.b_iso_or_equiv)
        elif adp_type == "Uovl":
            # FIXME: correct it
            u_iso = float(item_a_s.u_iso_or_equiv)
            b_iso = coeff*u_iso
        elif adp_type == "Umpe":
            # FIXME: correct it
            u_iso = float(item_a_s.u_iso_or_equiv)
            b_iso = float(8.*numpy.pi**2*u_iso)
        elif adp_type == "Uani":
            item_a_s_a = a_s_a[label_atom]
            beta = item_a_s_a.calc_beta(cell)
        elif adp_type == "Bovl":
            # FIXME: correct it
            b_iso = float(item_a_s.b_iso_or_equiv)
        elif adp_type == "Bani":
            item_a_s_a = a_s_a[label_atom]
            beta = (float(item_a_s_a.b_11), float(item_a_s_a.b_22),
                    float(item_a_s_a.b_33), float(item_a_s_a.b_12),
                    float(item_a_s_a.b_13), float(item_a_s_a.b_23))
        l_b_iso.append(b_iso)
        l_beta.append(beta)
    np_b_iso = numpy.array(l_b_iso, dtype=float)
    np_beta = numpy.array(l_beta, dtype=float)
    return np_b_iso, np_beta


# FIXME: full_space_group_symop is temporary slow solution.
def calc_f_nucl(
        index_hkl,
        space_group_symop_magn_operation: SpaceGroupSymopMagnOperationL,
        space_group_symop_magn_centering: SpaceGroupSymopMagnCenteringL,
        cell: Cell, atom_site: AtomSiteL, atom_site_aniso: AtomSiteAnisoL,
        flag_derivatives: bool = False):
    """
    Calculate nuclear structure factor. TEST.

    Keyword Arguments
    -----------------
        index_h, index_k, index_l: 1D numpy array of Miller indexes

    Output
    ------
        f_nucl: 1D numpy array of Nuclear structure factor

    Example
    -------
        >>> import numpy as np
        >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                      np.array([1,0],dtype=int)
        >>> f_nucl, der = crystal.calc_f_nucl(h, k, l)

    """
    dder = {}

    se_sgsmo = space_group_symop_magn_operation.get_sym_elems()
    se_sgsmc = space_group_symop_magn_centering.get_sym_elems()

    elem_symm_fs = calc_full_mag_elems(se_sgsmo, se_sgsmc)


    occupancy = numpy.array(atom_site.occupancy, dtype=float)
    x = numpy.array(atom_site.fract_x, dtype=float)
    y = numpy.array(atom_site.fract_y, dtype=float)
    z = numpy.array(atom_site.fract_z, dtype=float)

    fract_xyz = numpy.array([x, y, z], dtype=float)
    # FIXME: temporary solution
    try:
        atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
    except AttributeError:
        atom_multiplicity = calc_multiplicity(elem_symm_fs, fract_xyz)

    scat_length_neutron = numpy.array(atom_site.scat_length_neutron,
                                      dtype=complex)

    occ_mult = occupancy*atom_multiplicity

    r_11 = elem_symm_fs[4].astype(float)
    r_12 = elem_symm_fs[5].astype(float)
    r_13 = elem_symm_fs[6].astype(float)
    r_21 = elem_symm_fs[7].astype(float)
    r_22 = elem_symm_fs[8].astype(float)
    r_23 = elem_symm_fs[9].astype(float)
    r_31 = elem_symm_fs[10].astype(float)
    r_32 = elem_symm_fs[11].astype(float)
    r_33 = elem_symm_fs[12].astype(float)
    b_1 = elem_symm_fs[0].astype(float)/elem_symm_fs[3].astype(float)
    b_2 = elem_symm_fs[1].astype(float)/elem_symm_fs[3].astype(float)
    b_3 = elem_symm_fs[2].astype(float)/elem_symm_fs[3].astype(float)

    phase_3d = calc_phase_by_hkl_xyz_rb(
        index_hkl, x, y, z, r_11, r_12, r_13, r_21, r_22,
        r_23, r_31, r_32, r_33, b_1, b_2, b_3)

    b_iso, beta = calc_b_iso_beta(cell, atom_site, atom_site_aniso)

    dwf_3d = calc_dwf(cell, index_hkl, b_iso, beta, r_11,
                      r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

    hh = phase_3d*dwf_3d
    phase_2d = hh.sum(axis=2)  # sum over symmetry

    b_scat_2d = numpy.meshgrid(index_hkl[0], scat_length_neutron,
                               indexing="ij")[1]
    occ_mult_2d = numpy.meshgrid(index_hkl[0], occ_mult, indexing="ij")[1]

    hh = phase_2d * b_scat_2d * occ_mult_2d
    # nuclear structure factor in assymetric unit cell
    f_hkl_as = hh.sum(axis=1)*1./r_11.size

    # f_nucl = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
    #                                             f_hkl_as)

    return f_hkl_as, dder


def calc_f_mag(
        index_hkl,
        space_group_symop_magn_operation: SpaceGroupSymopMagnOperationL,
        space_group_symop_magn_centering: SpaceGroupSymopMagnCenteringL,
        cell: Cell, atom_site: AtomSiteL, atom_site_aniso: AtomSiteAnisoL,
        atom_site_scat: AtomSiteScatL, atom_site_moment: AtomSiteMomentL,
        flag_derivatives: bool = False, flag_only_orbital: bool = False):
    """Calculate magnetic structure factor.

    It's given in Cartesian coordianate system x||a* z||c
    (!!!!!!it should be checked!!!!!)

    Unity is 10**-12 cm
    """
    dder = {}
    index_h, index_k, index_l = index_hkl[0], index_hkl[1], index_hkl[2]

    m_x = numpy.array(atom_site_moment.crystalaxis_x, dtype=float)
    m_y = numpy.array(atom_site_moment.crystalaxis_y, dtype=float)
    m_z = numpy.array(atom_site_moment.crystalaxis_z, dtype=float)
    moment = numpy.array([m_x, m_y, m_z], dtype=float)

    l_it_a_s_mag = [atom_site[item.label] for item in atom_site_moment.items]
    atom_site_mag = AtomSiteL()
    atom_site_mag.items = l_it_a_s_mag

    l_x, l_y, l_z, l_occ, l_mult = [], [], [], [], []
    flag_calc_multiplicity = False
    for atom_site_item in atom_site_mag.items:
        l_x.append(atom_site_item.fract_x)
        l_y.append(atom_site_item.fract_y)
        l_z.append(atom_site_item.fract_z)
        l_occ.append(atom_site_item.occupancy)
        try:
            l_mult.append(atom_site_item.multiplicity)
        except AttributeError:
            flag_calc_multiplicity = True

    occupancy = numpy.array(l_occ, dtype=float)
    x = numpy.array(l_x, dtype=float)
    y = numpy.array(l_y, dtype=float)
    z = numpy.array(l_z, dtype=float)

    fract_xyz = numpy.array([x, y, z], dtype=float)

    b_iso, beta = calc_b_iso_beta(cell, atom_site_mag, atom_site_aniso)

    sthovl = cell.calc_sthovl(index_h, index_k, index_l)
    form_factor = atom_site_scat.calc_form_factor(
        sthovl, flag_only_orbital=flag_only_orbital)

    se_sgsmo = space_group_symop_magn_operation.get_sym_elems()
    se_sgsmc = space_group_symop_magn_centering.get_sym_elems()

    elem_symm_fs = calc_full_mag_elems(se_sgsmo, se_sgsmc)

    if flag_calc_multiplicity:
        atom_multiplicity = calc_multiplicity(elem_symm_fs[:13], fract_xyz)
    else:
        atom_multiplicity = numpy.array(l_mult, dtype=int)
    

    occ_mult = occupancy*atom_multiplicity

    r_11 = elem_symm_fs[4].astype(float)
    r_12 = elem_symm_fs[5].astype(float)
    r_13 = elem_symm_fs[6].astype(float)
    r_21 = elem_symm_fs[7].astype(float)
    r_22 = elem_symm_fs[8].astype(float)
    r_23 = elem_symm_fs[9].astype(float)
    r_31 = elem_symm_fs[10].astype(float)
    r_32 = elem_symm_fs[11].astype(float)
    r_33 = elem_symm_fs[12].astype(float)
    b_1 = elem_symm_fs[0].astype(float)/elem_symm_fs[3].astype(float)
    b_2 = elem_symm_fs[1].astype(float)/elem_symm_fs[3].astype(float)
    b_3 = elem_symm_fs[2].astype(float)/elem_symm_fs[3].astype(float)

    # [xyz, symm, mag_at]
    moments_3d = calc_moment_by_sym_elem(elem_symm_fs, moment)
    # print("moments_3d: \n", moments_3d)
    # print(get_str_for_sym_elem(sym_elems), end=2*"\n")

    phase_3d = calc_phase_by_hkl_xyz_rb(
        index_hkl, x, y, z, r_11, r_12, r_13, r_21, r_22,
        r_23, r_31, r_32, r_33, b_1, b_2, b_3)

    dwf_3d = calc_dwf(cell, index_hkl, b_iso, beta, r_11,
                      r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

    n_a = numpy.newaxis

    hh_1 = phase_3d*dwf_3d*(moments_3d[0, :, :]).transpose()[n_a, :, :]
    hh_2 = phase_3d*dwf_3d*(moments_3d[1, :, :]).transpose()[n_a, :, :]
    hh_3 = phase_3d*dwf_3d*(moments_3d[2, :, :]).transpose()[n_a, :, :]

    phase_2d_1 = hh_1.sum(axis=2)  # sum over symmetry
    phase_2d_2 = hh_2.sum(axis=2)  # sum over symmetry
    phase_2d_3 = hh_3.sum(axis=2)  # sum over symmetry

    occ_mult_2d = numpy.meshgrid(index_h, occ_mult, indexing="ij")[1]

    # print("phase_2d_1:", phase_2d_1)

    hh_1 = phase_2d_1 * form_factor * occ_mult_2d
    hh_2 = phase_2d_2 * form_factor * occ_mult_2d
    hh_3 = phase_2d_3 * form_factor * occ_mult_2d
    # nuclear structure factor in assymetric unit cell
    f_hkl_1 = hh_1.sum(axis=1)*1./r_11.size
    f_hkl_2 = hh_2.sum(axis=1)*1./r_11.size
    f_hkl_3 = hh_3.sum(axis=1)*1./r_11.size

    # f_nucl = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
    #                                             f_hkl_as)
    # f_mag should be given in Cartezian coordinate system (x||a*, z||c)
    m_m = cell.m_m_norm
    f_hkl_cart_1 = (m_m[0, 0]*f_hkl_1 + m_m[0, 1]*f_hkl_2 + m_m[0, 2]*f_hkl_3
                    )*0.2695
    f_hkl_cart_2 = (m_m[1, 0]*f_hkl_1 + m_m[1, 1]*f_hkl_2 + m_m[1, 2]*f_hkl_3
                    )*0.2695
    f_hkl_cart_3 = (m_m[2, 0]*f_hkl_1 + m_m[2, 1]*f_hkl_2 + m_m[2, 2]*f_hkl_3
                    )*0.2695
    f_mag = numpy.array([f_hkl_cart_1, f_hkl_cart_2, f_hkl_cart_3],
                        dtype=complex)
    return f_mag, dder
