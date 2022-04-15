# -*- coding: utf-8 -*-
"""
Calculations structure factor based on items and loops
"""
import numpy

from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import AtomSiteAnisoL

def calc_b_iso_beta(atom_site: AtomSiteL, atom_site_aniso: AtomSiteAnisoL, cell: Cell):
    """
    Calculate b_iso and beta_ij based on atom_site and atom_sites.
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
            if a_s_a is not None:
                item_a_s_a = a_s_a[label_atom]
                beta = item_a_s_a.calc_beta(cell)
        elif adp_type == "Bovl":
            # FIXME: correct it
            b_iso = float(item_a_s.b_iso_or_equiv)
        elif adp_type == "Bani":
            if a_s_a is not None:
                item_a_s_a = a_s_a[label_atom]
                beta = (float(item_a_s_a.b_11), float(item_a_s_a.b_22),
                        float(item_a_s_a.b_33), float(item_a_s_a.b_12),
                        float(item_a_s_a.b_13), float(item_a_s_a.b_23))
        l_b_iso.append(b_iso)
        l_beta.append(beta)
    np_b_iso = numpy.array(l_b_iso, dtype=float)
    np_beta = numpy.array(l_beta, dtype=float).transpose()
    return np_b_iso, np_beta

