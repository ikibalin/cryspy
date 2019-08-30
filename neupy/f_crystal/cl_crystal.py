"""
define classes to describe Crystal
"""
__author__ = 'ikibalin'
__version__ = "2019_08_30"
import os
import numpy


from pystar import CIFdata
from neupy.f_common.cl_fitable import Fitable
from neupy.f_crystal.cl_space_group import SpaceGroup
from neupy.f_crystal.cl_cell import Cell
from neupy.f_crystal.cl_atom_site import AtomSite
from neupy.f_crystal.cl_atom_site_aniso import AtomSiteAniso
from neupy.f_crystal.cl_atom_site_magnetism import AtomSiteMagnetism
from neupy.f_crystal.cl_atom_site_magnetism_aniso import AtomSiteMagnetismAniso



class Crystal(object):
    """
    Data items in the CRYSTAL category record details about
    crystal structure.
    
    Description in cif file:

    data_Fe3O4                                # object Crystal with label 'Fe3O4'
    _cell_angle_alpha 90.0                    # object Cell
    _cell_angle_beta 90.0
    _cell_angle_gamma 90.0
    _cell_length_a 8.56212()
    _cell_length_b 8.56212
    _cell_length_c 8.56212

    _space_group_it_coordinate_system_code 2  # object SpaceGroup
    _space_group_name_H-M_alt "F d -3 m"
    _space_group_IT_number    232

    loop_                                     # object AtomSite
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

    loop_                                     # object AtomType (optional)
    _atom_type_scat_length_neutron
    _atom_type_symbol
      0.945 Fe3+
     0.5803 O2-

    loop_                                     # object AtomSiteAniso (optional)
    _atom_site_aniso_U_11
    _atom_site_aniso_U_12
    _atom_site_aniso_U_13
    _atom_site_aniso_U_22
    _atom_site_aniso_U_23
    _atom_site_aniso_U_33
    _atom_site_aniso_label
     0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
     0.0 0.0 0.0 0.0 0.0 0.0 Fe3B

    loop_                                     # object AtomSiteMagnetism (optional)
    _atom_site_magnetism_label
    _atom_site_magnetism_lande
    _atom_site_magnetism_kappa
    Fe3A 2.0 1.0()
    Fe3B 2.0() 1.0

    loop_                                     # object AtomSiteMagnetismAniso (optional)
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
    def __init__(self, label='',  cell=None, space_group=None, atom_type=None,
                       atom_site=None, atom_site_aniso=None, 
                       atom_site_magnetism=None, atom_site_magnetism_aniso=None):
        super(Crystal, self).__init__()

        self.__label = ""
        self.__cell = None
        self.__space_group = None
        self.__atom_type = None
        self.__atom_site = None
        self.__atom_site_aniso = None
        self.__atom_site_magnetism = None
        self.__atom_site_magnetism_aniso = None

        self.label = label
        self.cell = cell
        self.space_group = space_group
        self.atom_type = atom_type
        self.atom_site = atom_site
        self.atom_site_aniso = atom_site_aniso
        self.atom_site_magnetism = atom_site_magnetism
        self.atom_site_magnetism_aniso = atom_site_magnetism_aniso
        
    def __repr__(self):
        ls_out = ["Crystal:"]
        ls_out.append(self.label)
        if self.cell is not None:
            ls_out.append("\n")
            ls_out.append(str(self.cell))
        if self.space_group is not None:
            ls_out.append("\n")
            ls_out.append(str(self.space_group))
        if self.atom_type is not None:
            ls_out.append("\n")
            ls_out.append(str(self.atom_type))
        if self.atom_site is not None:
            ls_out.append("\n")
            ls_out.append(str(self.atom_site))
        if self.atom_site_aniso is not None:
            ls_out.append("\n")
            ls_out.append(str(self.atom_site_aniso))
        if self.atom_site_magnetism is not None:
            ls_out.append("\n")
            ls_out.append(str(self.atom_site_magnetism))
        if self.atom_site_magnetism_aniso is not None:
            ls_out.append("\n")
            ls_out.append(str(self.atom_site_magnetism_aniso))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        The label is a unique identifier for a particular crystal. 

        Type: char
        """
        return self.__label
    @label.setter
    def label(self, x):
        if x is None:
            x_in = ""
        else:
            x_in = str(x).strip()
        self.__label = x_in

    @property
    def cell(self):
        """
        Data items in the CELL category record details about the
        crystallographic cell parameters and their measurement.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ccell.html
        """
        return self.__cell
    @cell.setter
    def cell(self, x):
        if isinstance(x, Cell):
            x_in = x
        elif isinstance(x, str):
            x_in = Cell()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to Cell")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into Cell")
        self.__cell = x_in

    @property
    def space_group(self):
        """
        Contains all the data items that refer to the space group as a
        whole, such as its name or crystal system. They may be looped,
        for example, in a list of space groups and their properties.

        Only a subset of the SPACE_GROUP category items appear in the
        core dictionary.  The remainder are found in the symmetry CIF
        dictionary.

        Space-group types are identified by their number as given in
        International Tables for Crystallography Vol. A. Specific
        settings of the space groups can be identified either by their
        Hall symbol or by specifying their symmetry operations.

        The commonly used Hermann-Mauguin symbol determines the
        space-group type uniquely but several different Hermann-Mauguin
        symbols may refer to the same space-group type. A
        Hermann-Mauguin symbol contains information on the choice of
        the basis, but not on the choice of origin.  Different formats
        for the Hermann-Mauguin symbol are found in the symmetry CIF
        dictionary.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Cspace_group.html
        """
        return self.__space_group
    @space_group.setter
    def space_group(self, x):
        if isinstance(x, SpaceGroup):
            x_in = x
        elif isinstance(x, str):
            x_in = SpaceGroup()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to SpaceGroup")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into SpaceGroup")
        self.__space_group = x_in

    @property
    def atom_type(self):
        """
        Data items in the ATOM_TYPE category record details about
        properties of the atoms that occupy the atom sites, such as the
        atomic scattering factors.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_type.html
        """
        return self.__atom_type
    @atom_type.setter
    def atom_type(self, x):
        """
        if isinstance(x, AtomType):
            x_in = x
        elif isinstance(x, str):
            x_in = AtomType()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to AtomType")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into AtomType")
        self.__atom_type = x_in
        """
        self.__atom_type = None #class AtomType is still not introduced


    @property
    def atom_site(self):
        """
        Data items in the ATOM_SITE category record details about
        the atom sites in a crystal structure, such as the positional
        coordinates.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html
        """
        return self.__atom_site
    @atom_site.setter
    def atom_site(self, x):
        if isinstance(x, AtomSite):
            x_in = x
        elif isinstance(x, str):
            x_in = AtomSite()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to AtomSite")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into AtomSite")
        self.__atom_site = x_in

    @property
    def atom_site_aniso(self):
        """
        Data items in the ATOM_SITE_ANISO category record details about
        the atom sites in a crystal structure, such as atomic displacement 
        parameters.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html
        """
        return self.__atom_site_aniso
    @atom_site_aniso.setter
    def atom_site_aniso(self, x):
        if isinstance(x, AtomSiteAniso):
            x_in = x
        elif isinstance(x, str):
            x_in = AtomSiteAniso()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to AtomSiteAniso")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into AtomSiteAniso")
        self.__atom_site_aniso = x_in


    @property
    def atom_site_magnetism(self):
        """
        Data items in the ATOM_SITE_MAGNETISM category record details about
        the magnetic parameters.

        """
        return self.__atom_site_magnetism
    @atom_site_magnetism.setter
    def atom_site_magnetism(self, x):
        if isinstance(x, AtomSiteMagnetism):
            x_in = x
        elif isinstance(x, str):
            x_in = AtomSiteMagnetism()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to AtomSiteMagnetism")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into AtomSiteMagnetism")
        self.__atom_site_magnetism = x_in


    @property
    def atom_site_magnetism_aniso(self):
        """
        Data items in the ATOM_SITE_MAGNETISM_ANISO category record details about
        the susceptibility tensor.
        """
        return self.__atom_site_magnetism_aniso
    @atom_site_magnetism_aniso.setter
    def atom_site_magnetism_aniso(self, x):
        if isinstance(x, AtomSiteMagnetismAniso):
            x_in = x
        elif isinstance(x, str):
            x_in = AtomSiteMagnetismAniso()
            flag = x_in.from_cif(x)
            if not(flag):
                self._show_message("A induced string can not be converted to AtomSiteMagnetismAniso")
                x_in = None
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("A type of induced element is not recognized to convert it into AtomSiteMagnetismAniso")
        self.__atom_site_magnetism_aniso = x_in

    @property
    def magnetism_aniso(self):
        return self.atom_site_magnetism_aniso
    @magnetism_aniso.setter
    def magnetism_aniso(self, x):
        self.atom_site_magnetism_aniso = x


    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given
        """
        cond = any([self.cell is not None, self.space_group is not None, self.atom_site is not None])
        return cond
    @property
    def is_variable(self):
        if self.cell is not None: flag_1 = self.cell.is_variable
        if self.space_group is not None: flag_2 = self.space_group.is_variable
        if self.atom_type is not None: flag_3 = self.atom_type.is_variable
        if self.atom_site is not None: flag_4 = self.atom_site.is_variable
        if self.atom_site_aniso is not None: flag_5 = self.atom_site_aniso.is_variable
        if self.atom_site_magnetism is not None: flag_6 = self.atom_site_magnetism.is_variable
        if self.atom_site_magnetism_aniso is not None: flag_7 = self.atom_site_magnetism_aniso.is_variable
        res = any([flag_1, flag_2, flag_3, flag_4, flag_5, flag_6, flag_7])
        return res


    def get_variables(self):
        if self.cell is not None: l_val_1 = self.cell.get_variables()
        if self.space_group is not None: l_val_2 = self.space_group.get_variables()
        if self.atom_type is not None: l_val_3 = self.atom_type.get_variables()
        if self.atom_site is not None: l_val_4 = self.atom_site.get_variables()
        if self.atom_site_aniso is not None: l_val_5 = self.atom_site_aniso.get_variables()
        if self.atom_site_magnetism is not None: l_val_6 = self.atom_site_magnetism.get_variables()
        if self.atom_site_magnetism_aniso is not None: l_val_7 = self.atom_site_magnetism_aniso.get_variables()

        l_variable = []
        l_variable.extend(l_val_1)
        l_variable.extend(l_val_2)
        l_variable.extend(l_val_3)
        l_variable.extend(l_val_4)
        l_variable.extend(l_val_5)
        l_variable.extend(l_val_6)
        l_variable.extend(l_val_7)
        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            str_1, str_2, str_3, str_4, str_5, str_6, str_7 = "", "", "", "", "", "", ""
            ls_out.append("data_{:}".format(self.label))
            if self.cell is not None: str_1 = self.cell.to_cif + "\n"
            if self.space_group is not None: str_2 = self.space_group.to_cif + "\n"
            if self.atom_type is not None: str_3 = self.atom_type.to_cif + "\n"
            if self.atom_site is not None: str_4 = self.atom_site.to_cif + "\n"
            if self.atom_site_aniso is not None: str_5 = self.atom_site_aniso.to_cif + "\n"
            if self.atom_site_magnetism is not None: str_6 = self.atom_site_magnetism.to_cif + "\n"
            if self.atom_site_magnetism_aniso is not None: str_7 = self.atom_site_magnetism_aniso.to_cif + "\n"
            ls_out.extend([str_1, str_2, str_3, str_4, str_5, str_6, str_7])
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_data = CIFdata()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        self.name = cif_data.name
        cif_values = cif_data.values
        if cif_values is not None:
            if cif_values.is_prefix("cell"):
                self.cell = str(cif_values)
            if cif_values.is_prefix("space_group"):
                self.space_group = str(cif_values)
        if cif_data.is_prefix("atom_type"): self.atom_type = str(cif_data["atom_type"])
        if cif_data.is_prefix("atom_site"): self.atom_site = str(cif_data["atom_site"])
        if cif_data.is_prefix("atom_site_aniso"): self.atom_site_aniso = str(cif_data["atom_site_aniso"])
        if cif_data.is_prefix("atom_site_magnetism"): self.atom_site_magnetism = str(cif_data["atom_site_magnetism"])
        if cif_data.is_prefix("atom_site_magnetism_aniso"): self.atom_site_magnetism_aniso = str(cif_data["atom_site_magnetism_aniso"])
        return True

