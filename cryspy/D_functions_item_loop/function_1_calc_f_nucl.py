"""Calculation of flip ratio."""
import numpy

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_phase_by_hkl_xyz_rb, calc_dwf

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAnisoL

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
        adp_type = item_a_s.adp_type
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
        index_h, index_k, index_l, full_space_group_symop,
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

    r_s_g_s = full_space_group_symop

    occupancy = numpy.array(atom_site.occupancy, dtype=float)
    x = numpy.array(atom_site.fract_x, dtype=float)
    y = numpy.array(atom_site.fract_y, dtype=float)
    z = numpy.array(atom_site.fract_z, dtype=float)

    atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
    scat_length_neutron = numpy.array(atom_site.scat_length_neutron,
                                      dtype=complex)

    occ_mult = occupancy*atom_multiplicity

    r_11 = r_s_g_s.numpy_r_11.astype(float)
    r_12 = r_s_g_s.numpy_r_12.astype(float)
    r_13 = r_s_g_s.numpy_r_13.astype(float)
    r_21 = r_s_g_s.numpy_r_21.astype(float)
    r_22 = r_s_g_s.numpy_r_22.astype(float)
    r_23 = r_s_g_s.numpy_r_23.astype(float)
    r_31 = r_s_g_s.numpy_r_31.astype(float)
    r_32 = r_s_g_s.numpy_r_32.astype(float)
    r_33 = r_s_g_s.numpy_r_33.astype(float)
    b_1 = r_s_g_s.numpy_b_1.astype(float)
    b_2 = r_s_g_s.numpy_b_2.astype(float)
    b_3 = r_s_g_s.numpy_b_3.astype(float)

    phase_3d = calc_phase_by_hkl_xyz_rb(
        index_h, index_k, index_l, x, y, z, r_11, r_12, r_13, r_21, r_22,
        r_23, r_31, r_32, r_33, b_1, b_2, b_3)

    b_iso, beta = calc_b_iso_beta(cell, atom_site, atom_site_aniso)

    dwf_3d = calc_dwf(cell, index_h, index_k, index_l, b_iso, beta, r_11,
                      r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

    hh = phase_3d*dwf_3d
    phase_2d = hh.sum(axis=2)  # sum over symmetry

    b_scat_2d = numpy.meshgrid(index_h, scat_length_neutron,
                               indexing="ij")[1]
    occ_mult_2d = numpy.meshgrid(index_h, occ_mult, indexing="ij")[1]

    hh = phase_2d * b_scat_2d * occ_mult_2d
    # nuclear structure factor in assymetric unit cell
    f_hkl_as = hh.sum(axis=1)*1./r_11.size

    # f_nucl = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
    #                                             f_hkl_as)

    return f_hkl_as, dder


# FIXME: full_space_group_symop is temporary slow solution.
def calc_f_mag(
        index_h, index_k, index_l, full_space_group_symop,
        cell: Cell, atom_site: AtomSiteL, atom_site_aniso: AtomSiteAnisoL,
        atom_site_moment: AtomSiteMomentL, atom_site_scat: AtomSiteScatL,
        flag_derivatives: bool = False, flag_only_orbital: bool = False):
    """Calculate magnetic structure factor. TEST."""
    dder = {}

    r_s_g_s = full_space_group_symop

    occupancy = numpy.array(atom_site.occupancy, dtype=float)
    x = numpy.array(atom_site.fract_x, dtype=float)
    y = numpy.array(atom_site.fract_y, dtype=float)
    z = numpy.array(atom_site.fract_z, dtype=float)

    atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
    scat_length_neutron = numpy.array(atom_site.scat_length_neutron,
                                      dtype=complex)

    sthovl = cell.calc_sthovl(index_h, index_k, index_l)
    form_factor = atom_site_scat.calc_form_factor(
        sthovl, flag_only_orbital=flag_only_orbital)


    m_x = numpy.array(atom_site_moment.crystalaxis_x, dtype=float)
    m_y = numpy.array(atom_site_moment.crystalaxis_y, dtype=float)
    m_z = numpy.array(atom_site_moment.crystalaxis_z, dtype=float)

    occ_mult = occupancy*atom_multiplicity

    r_11 = r_s_g_s.numpy_r_11.astype(float)
    r_12 = r_s_g_s.numpy_r_12.astype(float)
    r_13 = r_s_g_s.numpy_r_13.astype(float)
    r_21 = r_s_g_s.numpy_r_21.astype(float)
    r_22 = r_s_g_s.numpy_r_22.astype(float)
    r_23 = r_s_g_s.numpy_r_23.astype(float)
    r_31 = r_s_g_s.numpy_r_31.astype(float)
    r_32 = r_s_g_s.numpy_r_32.astype(float)
    r_33 = r_s_g_s.numpy_r_33.astype(float)
    b_1 = r_s_g_s.numpy_b_1.astype(float)
    b_2 = r_s_g_s.numpy_b_2.astype(float)
    b_3 = r_s_g_s.numpy_b_3.astype(float)

    phase_3d = calc_phase_by_hkl_xyz_rb(
        index_h, index_k, index_l, x, y, z, r_11, r_12, r_13, r_21, r_22,
        r_23, r_31, r_32, r_33, b_1, b_2, b_3)

    b_iso, beta = calc_b_iso_beta(cell, atom_site, atom_site_aniso)

    dwf_3d = calc_dwf(cell, index_h, index_k, index_l, b_iso, beta, r_11,
                      r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

    hh = phase_3d*dwf_3d
    phase_2d = hh.sum(axis=2)  # sum over symmetry

    b_scat_2d = numpy.meshgrid(index_h, scat_length_neutron,
                               indexing="ij")[1]
    occ_mult_2d = numpy.meshgrid(index_h, occ_mult, indexing="ij")[1]

    hh = phase_2d * b_scat_2d * occ_mult_2d
    # nuclear structure factor in assymetric unit cell
    f_hkl_as = hh.sum(axis=1)*1./r_11.size

    # f_nucl = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
    #                                             f_hkl_as)

    return f_hkl_as, dder
