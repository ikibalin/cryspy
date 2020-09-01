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
from cryspy.rhocif.cl_atom_local_axes import AtomLocalAxes, AtomLocalAxesL
from cryspy.rhocif.cl_atom_electron_configuration import AtomElectronConfiguration, AtomElectronConfigurationL



def calc_m_sigma(np_M, np_sigma):
    """
Estimate errors for equation:

    np_M * np_val * np_M^T
    """
    np_sig_1 = numpy.zeros((3, 3), dtype=float)
    for _i in range(3):
        for _j in range(3):
            res = []
            for _k in range(3):
                for _l in range(3):
                    res.append((np_sigma[_k, _l]*np_M[_i, _k]*np_M[_j, _l])**2)
            np_sig_1[_i, _j] = sum(res)
    return numpy.sqrt(np_sig_1)

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
    OPTIONAL_CLASSES = (AtomTypeL, AtomSiteAnisoL, AtomSiteSusceptibilityL, AtomSiteScatL, AtomTypeScatL, 
                        AtomLocalAxesL, AtomElectronConfigurationL)
    INTERNAL_CLASSES = ()
    def __init__(self, cell=None, atom_site=None, 
                 space_group=None, atom_type=None, atom_site_aniso=None, atom_site_susceptibility=None, 
                 atom_site_scat=None, atom_type_scat=None,  atom_local_axes=None, atom_electron_configuration=None,
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
        self.atom_local_axes = atom_local_axes

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


    @property
    def atom_local_axes(self):
        l_res = self[AtomLocalAxesL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_local_axes.setter
    def atom_local_axes(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomLocalAxesL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomLocalAxesL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    @property
    def atom_electron_configuration(self):
        l_res = self[AtomElectronConfigurationL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @atom_electron_configuration.setter
    def atom_electron_configuration(self, x):
        if x is None:
            pass
        elif isinstance(x, AtomElectronConfigurationL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, AtomElectronConfigurationL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    def apply_constraint(self)->bool:
        space_group = self.space_group
        space_group_wyckoff = space_group.space_group_wyckoff
        cell = self.cell
        flag_cell = cell.apply_constraint(space_group.bravais_type, space_group.it_coordinate_system_code)
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

    def calc_b_iso_beta(self):
        """
Calculate b_iso and beta_ij based on atom_site and atom_sites 
for each atom defined in atom_site
        """
        a_s = self.atom_site
        a_s_a = self.atom_site_aniso
        l_b_iso, l_beta = [], []
        coeff = float(8.*numpy.pi**2)
        cell = self.cell
        for item_a_s in a_s.item:
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
                #FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = coeff*u_iso
            elif adp_type == "Umpe":
                #FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = float(8.*numpy.pi**2*u_iso)
            elif adp_type == "Uani":
                item_a_s_a = a_s_a[label_atom]
                beta = item_a_s_a.calc_beta(cell)
            elif adp_type == "Bovl":
                #FIXME: correct it
                b_iso = float(item_a_s.b_iso_or_equiv)
            elif adp_type == "Bani":
                item_a_s_a = a_s_a[label_atom]
                beta  = (float(item_a_s_a.b_11), float(item_a_s_a.b_22), float(item_a_s_a.b_33),
                         float(item_a_s_a.b_12), float(item_a_s_a.b_13), float(item_a_s_a.b_23))
            l_b_iso.append(b_iso) 
            l_beta.append(beta)
        np_b_iso = numpy.array(l_b_iso, dtype=float)
        np_beta = numpy.array(l_beta, dtype=float)
        return np_b_iso, np_beta

    def calc_f_nucl(self, h, k, l):
        """
Calculate nuclear structure factor

Keyword arguments: 

    h, k, l: 1D numpy array of Miller indexes

Output: 

    f_nucl: 1D numpy array of Nuclear structure factor

Example:

>>> import numpy as np
>>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int), np.array([1,0],dtype=int)
>>> f_nucl = crystal.calc_f_nucl(h, k, l)
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
        
        r_11 = numpy.array(r_s_g_s.r_11, dtype=float)
        r_12 = numpy.array(r_s_g_s.r_12, dtype=float)
        r_13 = numpy.array(r_s_g_s.r_13, dtype=float)
        r_21 = numpy.array(r_s_g_s.r_21, dtype=float)
        r_22 = numpy.array(r_s_g_s.r_22, dtype=float)
        r_23 = numpy.array(r_s_g_s.r_23, dtype=float)
        r_31 = numpy.array(r_s_g_s.r_31, dtype=float)
        r_32 = numpy.array(r_s_g_s.r_32, dtype=float)
        r_33 = numpy.array(r_s_g_s.r_33, dtype=float)
        b_1 = numpy.array(r_s_g_s.b_1, dtype=float)
        b_2 = numpy.array(r_s_g_s.b_2, dtype=float)
        b_3 = numpy.array(r_s_g_s.b_3, dtype=float)


        phase_3d = CONSTANTS_AND_FUNCTIONS.calc_phase_by_hkl_xyz_rb(h, k, l, x, y, z, 
        r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3)

        b_iso, beta = self.calc_b_iso_beta()
        dwf_3d = CONSTANTS_AND_FUNCTIONS.calc_dwf(cell, h, k, l, b_iso, beta, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)


        hh = phase_3d*dwf_3d
        phase_2d = hh.sum(axis=2)#sum over symmetry

        b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]

        hh = phase_2d * b_scat_2d * occ_mult_2d
        f_hkl_as = hh.sum(axis=1)*1./len(r_11)#nuclear structure factor in assymetric unit cell

        f_nucl = space_group.calc_f_hkl_by_f_hkl_as(h, k, l, f_hkl_as)
        return f_nucl
 
        
    def calc_refln(self, h, k, l, flag_internal=True):
        """
Calculate Refln cryspy object where nuclear structure factor is stored.

Keyword arguments: 

    h, k, l: 1D numpy array of Miller indexes
    flag_internal: a flag to calculate or to use internal objects. 
                   It should be True if user call the function.
                   It's True by default.

Output: 

    refln: object cryspy.Refln

Example:

>>> import numpy as np
>>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int), np.array([1,0],dtype=int)
>>> refln = crystal.calc_refln(h, k, l)
>>> print(refln.to_cif())
        """
        f_nucl = self.calc_f_nucl(h, k, l)
        res = ReflnL()
        res.set_numpy_index_h(h)
        res.set_numpy_index_k(k)
        res.set_numpy_index_l(l)
        res.set_numpy_f_calc(f_nucl)
        if flag_internal:
            res.transform_numpy_arrays_to_items()
        return res


    def calc_susceptibility_moment_tensor(self, h, k, l, flag_only_orbital=False):
        """
Calculate susceptibility tensor and moment tensor in Cartesian orthogonal system (x||a*, z||c)

Keyword arguments: 

    h, k, l: 1D numpy array of Miller indexes
    flag_only_orbital default is False. When only orbital form-factor should be used put True

Output: 

    CHI_11, CHI_12, CHI_13, CHI_21, CHI_22, CHI_23, CHI_31, CHI_32, CHI_33: 1D numpy array of susceptibility tensor
    M_11, M_12, M_13, M_21, M_22, M_23, M_31, M_32, M_33: 1D numpy array of moment tensor

Example:

>>> import numpy as np
>>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int), np.array([1,0],dtype=int)
>>> CHI_M = crystal.calc_susceptibility_moment_tensor(h, k, l)
>>> CHI_11, CHI_12, CHI_13, CHI_21, CHI_22, CHI_23, CHI_31, CHI_32, CHI_33 = CHI_M[:9]
>>> M_11, M_12, M_13, M_21, M_22, M_23, M_31, M_32, M_33 = CHI_M[9:]
        """
        space_group = self.space_group
        r_s_g_s = space_group.reduced_space_group_symop

        cell = self.cell
        atom_site_scat = self.atom_site_scat
        atom_site_susceptibility = self.atom_site_susceptibility
        sthovl = cell.calc_sthovl(h, k, l)

        if atom_site_susceptibility is None:
            s_11, s_12, s_13 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            s_21, s_22, s_23 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            s_31, s_32, s_33 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            sm_11, sm_12, sm_13 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            sm_21, sm_22, sm_23 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            sm_31, sm_32, sm_33 = numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float), numpy.zeros(len(h), dtype=float)
            return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33, sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33

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

        r_11 = numpy.array(r_s_g_s.r_11, dtype=float)
        r_12 = numpy.array(r_s_g_s.r_12, dtype=float)
        r_13 = numpy.array(r_s_g_s.r_13, dtype=float)
        r_21 = numpy.array(r_s_g_s.r_21, dtype=float)
        r_22 = numpy.array(r_s_g_s.r_22, dtype=float)
        r_23 = numpy.array(r_s_g_s.r_23, dtype=float)
        r_31 = numpy.array(r_s_g_s.r_31, dtype=float)
        r_32 = numpy.array(r_s_g_s.r_32, dtype=float)
        r_33 = numpy.array(r_s_g_s.r_33, dtype=float)
        b_1 = numpy.array(r_s_g_s.b_1, dtype=float)
        b_2 = numpy.array(r_s_g_s.b_2, dtype=float)
        b_3 = numpy.array(r_s_g_s.b_3, dtype=float)


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
        form_factor = atom_site_scat.calc_form_factor(sthovl, flag_only_orbital=flag_only_orbital)

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
                (hh[:, 0], hh[:, 1], hh[:, 2],
                hh[:, 3], hh[:, 4], hh[:, 5],
                hh[:, 6], hh[:, 7], hh[:, 8]))

        sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33 = self._orto_matrix(
                (hh[:, 9], hh[:, 10], hh[:, 11],
                hh[:, 12], hh[:, 13], hh[:, 14], 
                hh[:, 15], hh[:, 16], hh[:, 17]))

        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33, sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33



    def calc_refln_susceptibility(self, h, k, l, flag_internal=True, flag_only_orbital=False):
        """
Calculate susceptibility tensor and moment tensor in Cartesian orthogonal system (x||a*, z||c)

Keyword arguments: 

    h, k, l: 1D numpy array of Miller indexes

Output: 

    ReflnSusceptibilityL object of cryspy
    
Example:

>>> import numpy as np
>>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int), np.array([1,0],dtype=int)
>>> refln_suscept = crystal.calc_refln_susceptibility(h, k, l)
        """
        CHI_M = self.calc_susceptibility_moment_tensor(h, k, l, flag_only_orbital=flag_only_orbital)

        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = CHI_M[:9]
        sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33 = CHI_M[9:]

        res = ReflnSusceptibilityL()
        res.set_numpy_index_h(h)
        res.set_numpy_index_k(k)
        res.set_numpy_index_l(l)
        res.set_numpy_chi_11_calc(s_11)
        res.set_numpy_chi_12_calc(s_12)
        res.set_numpy_chi_13_calc(s_13)
        res.set_numpy_chi_21_calc(s_21)
        res.set_numpy_chi_22_calc(s_22)
        res.set_numpy_chi_23_calc(s_23)
        res.set_numpy_chi_31_calc(s_31)
        res.set_numpy_chi_32_calc(s_32)
        res.set_numpy_chi_33_calc(s_33)
        res.set_numpy_moment_11_calc(sm_11)
        res.set_numpy_moment_12_calc(sm_12)
        res.set_numpy_moment_13_calc(sm_13)
        res.set_numpy_moment_21_calc(sm_21)
        res.set_numpy_moment_22_calc(sm_22)
        res.set_numpy_moment_23_calc(sm_23)
        res.set_numpy_moment_31_calc(sm_31)
        res.set_numpy_moment_32_calc(sm_32)
        res.set_numpy_moment_33_calc(sm_33)
        if flag_internal:
            res.transform_numpy_arrays_to_items()
        return res

    def _orto_matrix(self, l_ij):
        """
matrix l_ij is defined in coordinate system (a, b, c)
l_ij: = l_11, l_12, l_13, l_21, l_22, l_23, l_31, l_32, l_33
output matrix s_ij is defined in Cartezian coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        """
        cell = self.cell
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = cell.ortogonalize_matrix(l_ij)
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

Negtive eigenvalues of ellipsoid are replaced by positive.
        """
        cell = self.cell
        a_s_m_a = self.atom_site_susceptibility
        m_m_norm = cell.m_m_norm
        m_im_norm = numpy.linalg.inv(m_m_norm)
        l_res = []
        for _l, _11, _22, _33, _12, _13, _23 in zip(a_s_m_a.label,
                a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33, 
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):
            m_chi_loc = numpy.array([[_11, _12, _13],
                                     [_12, _22, _23],
                                     [_13, _23, _33]], dtype=float)
            m_chi_orto = numpy.matmul(numpy.matmul(m_m_norm, m_chi_loc), m_m_norm.transpose())
            val, vec = numpy.linalg.eig(m_chi_orto)
            D = numpy.diag(numpy.abs(val))#only positive eigenvalues
            R = numpy.array(vec).transpose()
            chi_as_u_orto = numpy.matmul(R.transpose(), numpy.matmul(D, R))
            chi_as_u_loc = numpy.matmul(numpy.matmul(m_im_norm, chi_as_u_orto), m_im_norm.transpose())
            l_res.append(chi_as_u_loc)
        return l_res
    
    def calc_main_axes_of_magnetization_ellipsoids(self):
        """
Calculates the susceptibility along the main axes of magnetization ellipsoid 
in Bohr magneton for each magnetic atom.

The main axes are given in Cartezian coordinate system (x||a*, z||c).

Output arguments:

    arg1: the list of the suscetibility along main axes for each atom
    arg2: the list of directions for each atom

Example:

>>> # crystal is defined object of cryspy library; 
>>> # type(crystal) is Crystal
>>> l_susceptibilities, l_directions = crystal.calc_main_axes_of_magnetization_ellipsoids()
>>> #cycle over magnetic atoms
>>> for susceptibilities, directions in zip(l_susceptibilities, l_directions):
>>>     #cycle over three main axes
>>>     for _val1, _direction in zip(susceptibilities, directions):
>>>         print(f"{_val1: 9.5f} mu_B/T along: {_direction[0]: 9.5f} {_direction[1]: 9.5f} {_direction[2]: 9.5f}")
>>>     print("")
>>> print("Cartezian coordinate system is x||a*, z||c.")
        """
        cell = self.cell
        a_s_m_a = self.atom_site_susceptibility
        if a_s_m_a is None:
            return None
        m_m_norm = cell.m_m_norm
        ll_moments = []
        ll_directions = []
        for _l, _11, _22, _33, _12, _13, _23 in zip(a_s_m_a.label,
                a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33, 
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):

            m_chi_loc = numpy.array([[_11, _12, _13],
                                     [_12, _22, _23],
                                     [_13, _23, _33]], dtype=float)
            m_chi_orto = numpy.matmul(numpy.matmul(m_m_norm, m_chi_loc), m_m_norm.transpose())

            val, vec = numpy.linalg.eig(m_chi_orto)
            flag_error = (math.isclose(_11.sigma, 0.) & math.isclose(_22.sigma, 0.) & math.isclose(_33.sigma, 0.) &
                          math.isclose(_12.sigma, 0.) & math.isclose(_13.sigma, 0.) & math.isclose(_23.sigma, 0.))
            if not(flag_error):
                np_sigma = numpy.array([[_11.sigma, _12.sigma, _13.sigma],
                                        [_12.sigma, _22.sigma, _23.sigma],
                                        [_13.sigma, _23.sigma, _33.sigma]], dtype=float)
                M1 = numpy.matmul(vec.transpose(), m_m_norm)
                M2 = calc_m_sigma(M1, np_sigma)
                l_sig = [sum(M2[:, 0]**2)**0.5, sum(M2[:, 1]**2)**0.5, sum(M2[:, 2]**2)**0.5] 
                val = [Fitable(value=_1, sigma=_2) for _1, _2 in zip(val, l_sig)]
            l_moments = list(val)
            l_directions = [vec[:,0], vec[:,1], vec[:,2]]
            ll_moments.append(l_moments)
            ll_directions.append(l_directions)
        return ll_moments, ll_directions

    def report_main_axes_of_magnetization_ellipsoids(self):
        """
Make a report about magnetization ellipsoids in string format.
Calculations are performed by get_main_axes_of_magnetization_ellipsoids method 
and calc_magnetization_ellipsoid method of the crystal object.
        """
        ls_out = []
        # crystal is defined object of cryspy library; 
        # type(crystal) is Crystal
        a_s_m_a = self.atom_site_susceptibility
        if a_s_m_a is None:
            return None
        l_susceptibilities, l_directions = self.calc_main_axes_of_magnetization_ellipsoids()
        l_chi_as_u = self.calc_magnetization_ellipsoid()
        #cycle over magnetic atoms
        for label, susceptibilities, directions, chi_as_u in zip(a_s_m_a.label, l_susceptibilities, l_directions, l_chi_as_u):
            ls_out.append(f"For `{label:}` the susceptibility is:")
            #cycle over three main axes
            for _val1, _direction in zip(susceptibilities, directions):
                s_val = f"{_val1: 9.5f}".rjust(9)
                ls_out.append(f"{s_val:} mu_B/T along: {_direction[0]: 9.5f} {_direction[1]: 9.5f} {_direction[2]: 9.5f}")
            ls_out.append(f"To plot magn. ellispoid as a thermal one:")
            ls_out.append(f"         U11, U22, U33: {chi_as_u[0, 0]:9.2f} {chi_as_u[1, 1]:9.2f} {chi_as_u[2, 2]:9.2f}")
            ls_out.append(f"         U12, U13, U23: {chi_as_u[0, 1]:9.2f} {chi_as_u[0, 2]:9.2f} {chi_as_u[1, 2]:9.2f}")
            ls_out.append("")
        ls_out.append("Cartezian coordinate system is x||a*, z||c.")

        return "\n".join(ls_out)

    def calc_magnetic_moments_with_field_loc(self, field_abc):
        """
Calculates the orientation of magetic moment for atoms at applied magnetic field.

The input magnetic field should be given in normalized unit cell (a/|a|, b/|b|, c/|c|)
The output magnetic moment are given in normalized unit cell (a/|a|, b/|b|, c/|c|)
        """
        np_field = numpy.array(field_abc, dtype=float)
        spgr = self.space_group

        cell = self.cell
        m_m_norm = cell.m_m_norm
        m_mt_norm_m_norm_field = numpy.matmul(numpy.matmul(m_m_norm.transpose(), m_m_norm), np_field)
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
                _moment = numpy.matmul(_chi, m_mt_norm_m_norm_field)

                l_lab_out.append(f"{_l:}_{_i_out+1:}")
                l_xyz_out.append(_xyz)
                l_moment_out.append(_moment)
        return l_lab_out, l_xyz_out, l_moment_out