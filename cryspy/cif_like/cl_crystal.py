__author__ = 'ikibalin'
__version__ = "2019_12_04"
import os
import math
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_data_constr import DataConstr
from cryspy.common.cl_fitable import Fitable

from cryspy.symcif.cl_space_group import SpaceGroup
import cryspy.cif_like.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS
from cryspy.corecif.cl_cell import Cell
from cryspy.corecif.cl_atom_site import AtomSite, AtomSiteL
from cryspy.corecif.cl_atom_type import AtomTypeL
from cryspy.corecif.cl_atom_site_aniso import AtomSiteAnisoL
from cryspy.corecif.cl_refln import Refln, ReflnL
from cryspy.magneticcif.cl_atom_site_susceptibility import AtomSiteSusceptibilityL
from cryspy.magneticcif.cl_atom_type_scat import AtomTypeScatL
from cryspy.magneticcif.cl_atom_site_scat import AtomSiteScatL
from cryspy.magneticcif.cl_refln_susceptibility import ReflnSusceptibility, ReflnSusceptibilityL

class Crystal(DataConstr):
    """
Data items in the CRYSTAL category record details about
crystal structure.

Description in cif file::

 data_Fe3O4             
 _cell_angle_alpha 90.0                    
 _cell_angle_beta 90.0
 _cell_angle_gamma 90.0
 _cell_length_a 8.56212()
 _cell_length_b 8.56212
 _cell_length_c 8.56212
 _space_group_it_coordinate_system_code 2  
 _space_group_IT_number    227
 
 loop_                                     
 _atom_site_adp_type
 _atom_site_B_iso_or_equiv
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_label
 _atom_site_occupancy
 _atom_site_type_symbol
  Uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
  Uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
  Uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-
 
 loop_                                     
 _atom_type_scat_length_neutron
 _atom_type_symbol
   0.945 Fe3+
  0.5803 O2-
 
 loop_
 _atom_site_aniso_U_11
 _atom_site_aniso_U_12
 _atom_site_aniso_U_13
 _atom_site_aniso_U_22
 _atom_site_aniso_U_23
 _atom_site_aniso_U_33
 _atom_site_aniso_label
  0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
  0.0 0.0 0.0 0.0 0.0 0.0 Fe3B
 
 loop_
 _atom_site_scat_label
 _atom_site_scat_lande
 Fe3A 2.0 
 Fe3B 2.0 
 
 loop_     
 _atom_site_susceptibility_label
 _atom_site_susceptibility_chi_type
 _atom_site_susceptibility_chi_11
 _atom_site_susceptibility_chi_12
 _atom_site_susceptibility_chi_13
 _atom_site_susceptibility_chi_22
 _atom_site_susceptibility_chi_23
 _atom_site_susceptibility_chi_33
  Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
  Fe3B Cani 3.041      0.0 0.0  3.041 0.0  3.041
    """
    MANDATORY_CLASSES = (SpaceGroup, Cell, AtomSiteL)
    OPTIONAL_CLASSES = (AtomTypeL, AtomSiteAnisoL, AtomSiteSusceptibilityL, AtomSiteScatL, AtomTypeScatL)
    INTERNAL_CLASSES = ()
    def __init__(self, cell=None, atom_site=None, 
                 space_group=None, atom_type=None, atom_site_aniso=None, atom_site_susceptibility=None, 
                 atom_site_scat=None, atom_type_scat=None,  
                 data_name=""):
        super(Crystal, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
                                      optional_classes=self.OPTIONAL_CLASSES,
                                      internal_classes=self.INTERNAL_CLASSES)
        self.cell = cell
        self.atom_site = atom_site
        self.data_name = data_name

        self.space_group = space_group
        self.atom_type = atom_type
        self.atom_site_aniso = atom_site_aniso
        self.atom_site_susceptibility = atom_site_susceptibility
        self.atom_site_scat = atom_site_scat
        self.atom_type_scat = atom_type_scat

        if self.is_defined:
            self.form_object

    @property
    def space_group(self):
        l_res = self[SpaceGroup]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @space_group.setter
    def space_group(self, x):
        if x is None:
            pass
        elif isinstance(x, SpaceGroup):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, SpaceGroup):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)


    @property
    def cell(self):
        l_res = self[Cell]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @cell.setter
    def cell(self, x):
        if x is None:
            pass
        elif isinstance(x, Cell):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Cell):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def atom_site(self):
        """
        """
        l_res = self[AtomSiteL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_site.setter
    def atom_site(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomSiteL):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, AtomSiteL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def atom_type(self):
        l_res = self[AtomTypeL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_type.setter
    def atom_type(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomTypeL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomTypeL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def atom_site_aniso(self):
        l_res = self[AtomSiteAnisoL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_site_aniso.setter
    def atom_site_aniso(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomSiteAnisoL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomSiteAnisoL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def atom_site_susceptibility(self):
        l_res = self[AtomSiteSusceptibilityL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_site_susceptibility.setter
    def atom_site_susceptibility(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomSiteSusceptibilityL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomSiteSusceptibilityL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def atom_type_scat(self):
        l_res = self[AtomTypeScatL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_type_scat.setter
    def atom_type_scat(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomTypeScatL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomTypeScatL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def atom_site_scat(self):
        l_res = self[AtomSiteScatL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_site_scat.setter
    def atom_site_scat(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomSiteScatL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomSiteScatL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    def __repr__(self):
        ls_out = ["Crystal:"]
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***", s_out, UserWarning, stacklevel=2)


    def apply_constraint(self)->bool:
        space_group = self.space_group
        space_group_wyckoff = space_group.space_group_wyckoff
        cell = self.cell
        flag_cell = cell.apply_constraint(space_group.bravais_type)
        atom_site = self.atom_site
        flag_atom_site = atom_site.apply_constraint(space_group_wyckoff)
        atom_site_aniso = self.atom_site_aniso
        flag_adp = True
        if atom_site_aniso is not None:
            flag_adp = atom_site_aniso.apply_space_group_constraint(atom_site, space_group)
        flag_sucs_1, flag_sucs_2, flag_sucs_3 = True, True, True
        atom_site_susceptibility = self.atom_site_susceptibility
        if atom_site_susceptibility is not None:
            flag_sucs_1 = atom_site_susceptibility.apply_chi_iso_constraint(cell)
            flag_sucs_2 = atom_site_susceptibility.apply_moment_iso_constraint(cell)
            flag_sucs_3 = atom_site_susceptibility.apply_space_group_constraint(atom_site, space_group)
            
        flag = all([flag_cell, flag_atom_site, flag_adp, flag_sucs_1, flag_sucs_2, flag_sucs_3])
        return flag


    def calc_refln(self, h, k, l):
        """
calculate nuclear structure factor

FIXME: introduce Debye-Waller factor
        """
        space_group = self.space_group
        r_s_g_s = space_group.reduced_space_group_symop

        cell = self.cell
        atom_site = self.atom_site
        atom_site_aniso = self.atom_site_aniso

        occupancy = numpy.array(atom_site.occupancy, dtype=float)
        x = numpy.array(atom_site.fract_x, dtype=float)
        y = numpy.array(atom_site.fract_y, dtype=float)
        z = numpy.array(atom_site.fract_z, dtype=float)

        atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
        scat_length_neutron = numpy.array(atom_site.scat_length_neutron, dtype=complex)
        

        occ_mult = occupancy*atom_multiplicity 

        r_11, r_12, r_13 = r_s_g_s.r_11, r_s_g_s.r_12, r_s_g_s.r_13
        r_21, r_22, r_23 = r_s_g_s.r_21, r_s_g_s.r_22, r_s_g_s.r_23
        r_31, r_32, r_33 = r_s_g_s.r_31, r_s_g_s.r_32, r_s_g_s.r_33
        b_1, b_2, b_3 = r_s_g_s.b_1, r_s_g_s.b_2, r_s_g_s.b_3

        phase_3d = CONSTANTS_AND_FUNCTIONS.calc_phase_by_hkl_xyz_rb(h, k, l, x, y, z, 
        r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3)

        flag_adp = atom_site_aniso is not None
        flag_adp = False
        if flag_adp:
            if atom_site.is_defined_attribute("b_iso_or_equiv"):
                b_iso = numpy.array(atom_site.b_iso_or_equiv, dtype=float) 
            else:
                b_iso = numpy.zeros(shape=(len(atom_site.label)), dtype=float)
            beta = atom_site_aniso.calc_beta(cell)
            dwf_3d = CONSTANTS_AND_FUNCTIONS.calc_dwf(cell, h, k, l, b_iso, beta, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        else:
            dwf_3d = numpy.ones(phase_3d.shape, dtype=float)

        hh = phase_3d*dwf_3d
        phase_2d = hh.sum(axis=2)#sum over symmetry

        b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]

        hh = phase_2d * b_scat_2d * occ_mult_2d
        f_hkl_as = hh.sum(axis=1)*1./len(r_11)#nuclear structure factor in assymetric unit cell

        f_nucl = space_group.calc_f_hkl_by_f_hkl_as(h, k, l, f_hkl_as)
   
        

        item=[Refln(index_h=_1, index_k=_2, index_l=_3, f_calc=_4) for _1, _2, _3, _4 in zip(h, k, l, f_nucl)]
        res = ReflnL(item=item)
        return res


    def calc_refln_susceptibility(self, h, k, l):
        """
calculate susceptibility structure factor tensor

FIXME: introduce Debye-Waller factor
        """
        space_group = self.space_group
        r_s_g_s = space_group.reduced_space_group_symop

        cell = self.cell
        atom_site_scat = self.atom_site_scat
        atom_site_susceptibility = self.atom_site_susceptibility
        sthovl = cell.calc_sthovl(h, k, l)

        if atom_site_susceptibility is None:
            item = [ReflnSusceptibility(index_h=_h, index_k=_k, index_l=_l, sintlambda=_s,
                                    chi_11_calc = 0., chi_12_calc = 0., chi_13_calc = 0., 
                                    chi_21_calc = 0., chi_22_calc = 0., chi_23_calc = 0., 
                                    chi_31_calc = 0., chi_32_calc = 0., chi_33_calc = 0., 
                                    moment_11_calc = 0., moment_12_calc = 0., moment_13_calc = 0., 
                                    moment_21_calc = 0., moment_22_calc = 0., moment_23_calc = 0., 
                                    moment_31_calc = 0., moment_32_calc = 0., moment_33_calc = 0.)
                for _h, _k, _l, _s in zip(h, k, l, sthovl)]
            res = ReflnSusceptibilityL(item=item)
            return res

        chi_11 = numpy.array(atom_site_susceptibility.chi_11, float)
        chi_22 = numpy.array(atom_site_susceptibility.chi_22, float)
        chi_33 = numpy.array(atom_site_susceptibility.chi_33, float)
        chi_12 = numpy.array(atom_site_susceptibility.chi_12, float) 
        chi_13 = numpy.array(atom_site_susceptibility.chi_13, float) 
        chi_23 = numpy.array(atom_site_susceptibility.chi_23, float) 

        chi_11[numpy.isnan(chi_11)] = 0.
        chi_22[numpy.isnan(chi_22)] = 0.
        chi_33[numpy.isnan(chi_33)] = 0.
        chi_12[numpy.isnan(chi_12)] = 0.
        chi_13[numpy.isnan(chi_13)] = 0.
        chi_23[numpy.isnan(chi_23)] = 0.

        moment_11 = numpy.array(atom_site_susceptibility.moment_11, float)
        moment_22 = numpy.array(atom_site_susceptibility.moment_22, float)
        moment_33 = numpy.array(atom_site_susceptibility.moment_33, float)
        moment_12 = numpy.array(atom_site_susceptibility.moment_12, float) 
        moment_13 = numpy.array(atom_site_susceptibility.moment_13, float) 
        moment_23 = numpy.array(atom_site_susceptibility.moment_23, float) 

        moment_11[numpy.isnan(moment_11)] = 0.
        moment_22[numpy.isnan(moment_22)] = 0.
        moment_33[numpy.isnan(moment_33)] = 0.
        moment_12[numpy.isnan(moment_12)] = 0.
        moment_13[numpy.isnan(moment_13)] = 0.
        moment_23[numpy.isnan(moment_23)] = 0.

        atom_site = self.atom_site
        atom_site_scat.load_atom_type_scat_by_atom_site(atom_site)


        np_x_y_z_occ_mult = numpy.array([(atom_site[_item.label].fract_x, atom_site[_item.label].fract_y,
                               atom_site[_item.label].fract_z, atom_site[_item.label].occupancy, 
                               atom_site[_item.label].multiplicity)
                               for _item in atom_site_susceptibility.item], 
                               dtype = float)

        atom_site_aniso = self.atom_site_aniso

        x = np_x_y_z_occ_mult[:, 0]
        y = np_x_y_z_occ_mult[:, 1]
        z = np_x_y_z_occ_mult[:, 2]
        occupancy = np_x_y_z_occ_mult[:, 3]
        atom_multiplicity = np_x_y_z_occ_mult[:, 4]

        occ_mult = occupancy*atom_multiplicity 

        r_11, r_12, r_13 = r_s_g_s.r_11, r_s_g_s.r_12, r_s_g_s.r_13
        r_21, r_22, r_23 = r_s_g_s.r_21, r_s_g_s.r_22, r_s_g_s.r_23
        r_31, r_32, r_33 = r_s_g_s.r_31, r_s_g_s.r_32, r_s_g_s.r_33
        b_1, b_2, b_3 = r_s_g_s.b_1, r_s_g_s.b_2, r_s_g_s.b_3

        phase_3d = CONSTANTS_AND_FUNCTIONS.calc_phase_by_hkl_xyz_rb(h, k, l, x, y, z, 
        r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3)

        flag_adp = atom_site_aniso is not None
        flag_adp = False
        if flag_adp:
            if atom_site.is_defined_attribute("b_iso_or_equiv"):
                b_iso = numpy.array(atom_site.b_iso_or_equiv, dtype=float) 
            else:
                b_iso = numpy.zeros(shape=(len(atom_site.label)), dtype=float)
            beta = atom_site_aniso.calc_beta(cell)
            dwf_3d = CONSTANTS_AND_FUNCTIONS.calc_dwf(cell, h, k, l, b_iso, beta, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        else:
            dwf_3d = numpy.ones(phase_3d.shape, dtype=float)

        hh = phase_3d*dwf_3d
        
        #phase_2d = hh.sum(axis=2)#sum over symmetry

        #b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        form_factor = atom_site_scat.calc_form_factor(sthovl)

        #dimensions: hkl, magnetic atoms, reduced symmetry operators
        ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 = \
        CONSTANTS_AND_FUNCTIONS.calc_form_factor_tensor_susceptibility(
            chi_11, chi_22, chi_33, chi_12, chi_13, chi_23, 
            r_s_g_s, form_factor, cell, h, k, l)

        ffm_11, ffm_12, ffm_13, ffm_21, ffm_22, ffm_23, ffm_31, ffm_32, ffm_33 = \
        CONSTANTS_AND_FUNCTIONS.calc_form_factor_tensor_susceptibility(
            moment_11, moment_22, moment_33, moment_12, moment_13, moment_23, 
            r_s_g_s, form_factor, cell, h, k, l)


        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]

        #dimensions: hkl, number of atoms, reduced symmetry operators, 18 elements of susceptribility tensor
        hh = numpy.stack([ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33, 
                          ffm_11, ffm_12, ffm_13, ffm_21, ffm_22, ffm_23, ffm_31, ffm_32, ffm_33], axis=-1)

        b_scat_4d = hh

        hh = (phase_3d * dwf_3d * occ_mult_2d[:, :, numpy.newaxis])[:, :, :, numpy.newaxis] * b_scat_4d
        f_hkl_as_2d = (hh.sum(axis=2)*1./len(r_11)).sum(axis=1)#nuclear structure factor in assymetric unit cell

        f_2d = space_group.calc_f_hkl_by_f_hkl_as(h, k, l, f_hkl_as_2d)
        hh = 0.2695*f_2d
        #sft_ij form the structure factor tensor in local coordinate system (ia, ib, ic)
        #chi in 10-12 cm; chim in muB (it is why here 0.2695)
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = self._orto_matrix(
                cell,
                hh[:, 0], hh[:, 1], hh[:, 2],
                hh[:, 3], hh[:, 4], hh[:, 5],
                hh[:, 6], hh[:, 7], hh[:, 8])

        sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33 = self._orto_matrix(
                cell,
                hh[:, 9], hh[:, 10], hh[:, 11], 
                hh[:, 12], hh[:, 13], hh[:, 14], 
                hh[:, 15], hh[:, 16], hh[:, 17])

        item = [ReflnSusceptibility(index_h=_h, index_k=_k, index_l=_l, sintlambda=_s,
                                    chi_11_calc = _s_11, chi_12_calc = _s_12, chi_13_calc = _s_13, 
                                    chi_21_calc = _s_21, chi_22_calc = _s_22, chi_23_calc = _s_23, 
                                    chi_31_calc = _s_31, chi_32_calc = _s_32, chi_33_calc = _s_33, 
                                    moment_11_calc = _sm_11, moment_12_calc = _sm_12, moment_13_calc = _sm_13, 
                                    moment_21_calc = _sm_21, moment_22_calc = _sm_22, moment_23_calc = _sm_23, 
                                    moment_31_calc = _sm_31, moment_32_calc = _sm_32, moment_33_calc = _sm_33)
        for _h, _k, _l, _s, _s_11, _s_12, _s_13, _s_21, _s_22, _s_23, _s_31, _s_32, _s_33, 
            _sm_11, _sm_12, _sm_13, _sm_21, _sm_22, _sm_23, _sm_31, _sm_32, _sm_33 in zip(
            h, k, l, sthovl, s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33,
            sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33)    
        ]
        res = ReflnSusceptibilityL(item=item)
        return res


    def _orto_matrix(self, cell, l_11, l_12, l_13, l_21, l_22, l_23, l_31, 
                     l_32, l_33):
        """
rewrite matrix l_ij defined in coordinate (ia, ib, ic) to matrix s_ij, 
which is denined in Chartesian coordinate system, such as:
x||ia, y in blane (ia, ib), z perpendicular to that plane.
...

...
representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
expressions are taken from international tables
matrix_ib is inversed matrix B

ia, ib, ic is inversed unit cell parameters (it can be estimated from matrix matrix_ib)

.. math::

    X = B x, x = iB X
    xT*CHI*x = XT iBT CHI iB X

output chiLOC = iBT CHI iB
        """
        m_ib_norm = cell.m_ib_norm
        m_ibt_norm = m_ib_norm.transpose()
        
        r11, r12, r13 = m_ibt_norm[0, 0], m_ibt_norm[0, 1], m_ibt_norm[0, 2]
        r21, r22, r23 = m_ibt_norm[1, 0], m_ibt_norm[1, 1], m_ibt_norm[1, 2]
        r31, r32, r33 = m_ibt_norm[2, 0], m_ibt_norm[2, 1], m_ibt_norm[2, 2]
        
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = CONSTANTS_AND_FUNCTIONS.calc_mRmCmRT(
                r11, r12, r13, r21, r22, r23, r31, r32, r33,
                l_11, l_12, l_13, l_21, l_22, l_23, l_31, l_32, l_33)        

        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33

    
    def calc_hkl(self, sthol_min, sthovl_max):
        cell = self.cell
        space_group = self.space_group
        res = cell.calc_hkl(space_group, sthol_min, sthovl_max)
        return res

    def calc_hkl_in_range(self, sthol_min, sthovl_max):
        cell = self.cell
        res = cell.calc_hkl_in_range(sthol_min, sthovl_max)
        return res
    
    def calc_magnetization_ellipsoid(self):
        """
Magnetization ellipsoids are given in the same coordinate system as U_ij (anisotropic Debye-Waller factor)
        """
        cell = self.cell
        a_s_m_a = self.atom_site_susceptibility
        m_ib_norm = cell.m_ib_norm
        m_ibt_norm = numpy.transpose(m_ib_norm)
        l_res = []
        for _l, _11, _22, _33, _12, _13, _23 in zip(a_s_m_a.label,
                a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33, 
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):
            m_chi = numpy.array([[_11, _12, _13],
                                 [_12, _22, _23],
                                 [_13, _23, _33]], dtype=float)
            _m1 = numpy.matmul(m_chi, m_ib_norm)
            _m2 = numpy.matmul(m_ibt_norm, m_chi)
            _m_u = numpy.matmul(_m1, _m2)
            l_res.append(_m_u)
        return l_res
    
    def calc_main_axes_of_magnetization_ellipsoids(self):
        cell = self.cell
        a_s_m_a = self.atom_site_susceptibility
        m_ib_norm = cell.m_ib_norm
        m_ibt_norm = numpy.transpose(m_ib_norm)
        ll_moments = []
        ll_directions = []
        for _l, _11, _22, _33, _12, _13, _23 in zip(a_s_m_a.label,
                a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33, 
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):

            s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = self._orto_matrix(cell, _11, _12, _13, _12, _22, _23, _13, _23, _33)
            m_chi_norm = numpy.array([[s_11, s_12, s_13],
                                      [s_12, s_22, s_23],
                                      [s_13, s_23, s_33]], dtype=float)
            eig, mat = numpy.linalg.eig(m_chi_norm)
            l_moments = list(eig)
            l_directions = [mat[:,0], mat[:,1], mat[:,2]]
            ll_moments.append(l_moments)
            ll_directions.append(l_directions)
        return ll_moments, ll_directions
    
    def calc_magnetic_moments_with_field_loc(self, field_abc):
        """
:FIXME: IN GENERAL CASE NOT CORRECT
        """
        np_field = numpy.array(field_abc, dtype=float)
        spgr = self.space_group
        
        a_s = self.atom_site
        a_s_m_a = self.atom_site_susceptibility
        l_lab_out, l_xyz_out, l_moment_out = [], [], []
        for _l, _11, _22, _33, _12, _13, _23 in zip(a_s_m_a.label,
                a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33, 
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):
            m_chi = numpy.array([[_11, _12, _13],
                                 [_12, _22, _23],
                                 [_13, _23, _33]], dtype=float)
            x, y, z = float(a_s[_l].fract_x), float(a_s[_l].fract_y), float(a_s[_l].fract_z)
            l_out = spgr.calc_rotated_matrix_for_position(m_chi, x, y, z)
            for _i_out, _out in enumerate(l_out):
                _xyz = _out[0]
                _chi = _out[1]
                _moment = numpy.matmul(_chi, np_field)
                l_lab_out.append(f"{_l:}_{_i_out+1:}")
                l_xyz_out.append(_xyz)
                l_moment_out.append(_moment)
        return l_lab_out, l_xyz_out, l_moment_out