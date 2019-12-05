"""
define classes to describe Crystal
"""
__author__ = 'ikibalin'
__version__ = "2019_12_04"
import os
import math
import numpy
from pycifstar import Global

import cryspy.f_crystal.corecif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS

import warnings
from typing import List, Tuple
from cryspy.f_common.cl_item_constr import ItemConstr
from cryspy.f_common.cl_loop_constr import LoopConstr
from cryspy.f_common.cl_data_constr import DataConstr

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_crystal.symcif.cl_space_group import SpaceGroup
from cryspy.f_crystal.corecif.cl_cell import Cell
from cryspy.f_crystal.corecif.cl_atom_site import AtomSiteL

class Crystal(DataConstr):
    """
Crystal:
===============
Data items in the CRYSTAL category record details about
crystal structure.

Description in cif file:
---------------------------
data_Fe3O4                                
_cell_angle_alpha 90.0                    
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_cell_length_a 8.56212()
_cell_length_b 8.56212
_cell_length_c 8.56212
_space_group_it_coordinate_system_code 2  
_space_group_name_H-M_alt "F d -3 m"
_space_group_IT_number    232

loop_                                     
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_label
_atom_site_occupancy
_atom_site_type_symbol
 uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
 uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
 uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-

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
_atom_site_magnetism_label
_atom_site_magnetism_lande
_atom_site_magnetism_kappa
Fe3A 2.0 1.0()
Fe3B 2.0() 1.0

loop_     
_atom_site_magnetism_aniso_label
_atom_site_magnetism_aniso_chi_type
_atom_site_magnetism_aniso_chi_11
_atom_site_magnetism_aniso_chi_12
_atom_site_magnetism_aniso_chi_13
_atom_site_magnetism_aniso_chi_22
_atom_site_magnetism_aniso_chi_23
_atom_site_magnetism_aniso_chi_33
 Fe3A cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
 Fe3B cani 3.041      0.0 0.0  3.041 0.0  3.041
    """
    MANDATORY_CLASSES = (SpaceGroup, Cell, AtomSiteL)
    OPTIONAL_CLASSES = ()
    INTERNAL_CLASSES = ()
    def __init__(self, cell=None, atom_site=None, data_name=""):
        super(Crystal, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
                                      optional_classes=self.OPTIONAL_CLASSES,
                                      internal_classes=self.INTERNAL_CLASSES)
        self.cell = cell
        self.atom_site = atom_site
        self.data_name = data_name
        if self.is_defined:
            self.form_object

    @property
    def space_group(self):
        """
        """
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
    def cell(self):
        """
        """
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


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***", s_out, UserWarning, stacklevel=2)

    
    @property
    def is_variable(self):
        """
        Output: True if there is any refined parameter
        """
        res = any([_obj.is_variable for _obj in self.mandatory_objs] +
                  [_obj.is_variable for _obj in self.optional_objs])
        return res
        
    def get_variables(self):
        """
        Output: the list of the refined parameters
        """
        l_variable = []
        for _obj in self.mandatory_objs:
            l_variable.extend(_obj.get_variables())
        for _obj in self.optional_objs:
            l_variable.extend(_obj.get_variables())
        return l_variable

    def apply_constraint(self)->bool:
        type_cell = self.space_group.type_cell
        space_group_wyckoff = self.space_group.space_group_wyckoff
        flag = self.cell.apply_constraint(type_cell)
        return flag

