"""
define classes to describe Crystal
"""
__author__ = 'ikibalin'
__version__ = "2019_09_09"
import os
import numpy
import copy

from pycifstar import Data
from cryspy.f_common.cl_fitable import Fitable
from .cl_space_group import SpaceGroup
from .cl_cell import Cell
from .cl_atom_site import AtomSite
from .cl_atom_site_aniso import AtomSiteAniso
from .cl_atom_site_magnetism import AtomSiteMagnetism
from .cl_atom_site_magnetism_aniso import AtomSiteMagnetismAniso

from .cl_magnetism import calc_mRmCmRT

from .cl_refln import Refln


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
        ls_out.append("label: "+self.label+"\n")
        if self.cell is not None:
            ls_out.append(str(self.cell)+"\n")
        if self.space_group is not None:
            ls_out.append(str(self.space_group)+"\n")
        if self.atom_type is not None:
            ls_out.append(str(self.atom_type)+"\n")
        if self.atom_site is not None:
            ls_out.append(str(self.atom_site)+"\n")
        if self.atom_site_aniso is not None:
            ls_out.append(str(self.atom_site_aniso)+"\n")
        if self.atom_site_magnetism is not None:
            ls_out.append(str(self.atom_site_magnetism)+"\n")
        if self.atom_site_magnetism_aniso is not None:
            ls_out.append(str(self.atom_site_magnetism_aniso)+"\n")
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
        flag_1, flag_2, flag_3, flag_4, flag_5, flag_6, flag_7 = False, False, False, False, False, False, False
        if self.cell is not None: flag_1 = self.cell.is_variable
        #if self.space_group is not None: flag_2 = self.space_group.is_variable
        if self.atom_type is not None: flag_3 = self.atom_type.is_variable
        if self.atom_site is not None: flag_4 = self.atom_site.is_variable
        if self.atom_site_aniso is not None: flag_5 = self.atom_site_aniso.is_variable
        if self.atom_site_magnetism is not None: flag_6 = self.atom_site_magnetism.is_variable
        if self.atom_site_magnetism_aniso is not None: flag_7 = self.atom_site_magnetism_aniso.is_variable
        res = any([flag_1, flag_2, flag_3, flag_4, flag_5, flag_6, flag_7])
        return res


    def get_variables(self):
        l_val_1, l_val_2, l_val_3, l_val_4, l_val_5, l_val_6, l_val_7 = [], [], [], [], [], [], []
        if self.cell is not None: l_val_1 = self.cell.get_variables()
        #if self.space_group is not None: l_val_2 = self.space_group.get_variables()
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
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        self.label = cif_data.name
        cif_values = cif_data.items
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

    def calc_sf(self, h, k, l):
        """
        calculate nuclear structure factor and components of structure factor tensor
        """

        space_group = self.space_group
        cell = self.cell
        atom_site = self.atom_site
        atom_site_aniso = self.atom_site_aniso
        atom_site_magnetism = self.atom_site_magnetism
        atom_site_magnetism_aniso = self.atom_site_magnetism_aniso


        occupancy = numpy.array(atom_site.occupancy, dtype=float)
        scat_length_neutron = numpy.array(atom_site.scat_length_neutron, dtype=complex)


        fract = atom_site._form_fract()
        flag_adp = atom_site_aniso is not None
        flag_magnetism = atom_site_magnetism_aniso is not None
        if flag_adp:
            adp = atom_site_aniso._form_adp(atom_site)
        if flag_magnetism:
            magnetism = atom_site_magnetism_aniso._form_magnetism(atom_site, atom_site_magnetism)


        x, y, z = fract.x, fract.y, fract.z

        atom_multiplicity = space_group.calc_atom_mult(x, y, z)
        occ_mult = occupancy*atom_multiplicity 
        

        phase_3d = fract.calc_phase(space_group, h, k, l)#3d object
        if flag_adp:
            dwf_3d = adp.calc_dwf(space_group, cell, h, k, l)
        else:
            dwf_3d = numpy.ones(phase_3d.shape, dtype=float)

        hh = phase_3d*dwf_3d

        if flag_magnetism:
            ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 = \
                   magnetism.calc_form_factor_tensor(space_group, cell, h, k, l)
        else:
            np_zeros = numpy.zeros(phase_3d.shape, dtype=float)
            ff_11, ff_12, ff_13 = np_zeros, np_zeros, np_zeros
            ff_21, ff_22, ff_23 = np_zeros, np_zeros, np_zeros 
            ff_31, ff_32, ff_33 = np_zeros, np_zeros, np_zeros


        phase_2d = hh.sum(axis=2)

        ft_11 = (ff_11*hh).sum(axis=2)
        ft_12 = (ff_12*hh).sum(axis=2)
        ft_13 = (ff_13*hh).sum(axis=2)
        ft_21 = (ff_21*hh).sum(axis=2)
        ft_22 = (ff_22*hh).sum(axis=2)
        ft_23 = (ff_23*hh).sum(axis=2)
        ft_31 = (ff_31*hh).sum(axis=2)
        ft_32 = (ff_32*hh).sum(axis=2)
        ft_33 = (ff_33*hh).sum(axis=2)
        
        b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]
        
        l_el_symm = space_group.el_symm
        l_orig = space_group.orig
        centr = space_group.centr

        #calculation of nuclear structure factor        
        hh = phase_2d * b_scat_2d * occ_mult_2d
        f_hkl_as = hh.sum(axis=1)*1./len(l_el_symm)
        
        orig_x = [hh[0] for hh in l_orig]
        orig_y = [hh[1] for hh in l_orig]
        orig_z = [hh[2] for hh in l_orig]
        
        np_h, np_orig_x = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k, np_orig_y = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l, np_orig_z = numpy.meshgrid(l, orig_z, indexing = "ij")
        
        np_orig_as = numpy.exp(2*numpy.pi*1j*(np_h*np_orig_x+np_k*np_orig_y+np_l*np_orig_z))
        f_hkl_as = f_hkl_as*np_orig_as.sum(axis=1)*1./len(l_orig)

        if (centr):
            orig = space_group.p_centr
            f_nucl = 0.5*(f_hkl_as+f_hkl_as.conjugate()*numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2])))
        else:
            f_nucl = f_hkl_as

        #calculation of structure factor tensor
        sft_as_11 = (ft_11 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_12 = (ft_12 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_13 = (ft_13 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_21 = (ft_21 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_22 = (ft_22 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_23 = (ft_23 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_31 = (ft_31 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_32 = (ft_32 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_33 = (ft_33 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)

        sft_as_11 = sft_as_11 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_12 = sft_as_12 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_13 = sft_as_13 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_21 = sft_as_21 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_22 = sft_as_22 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_23 = sft_as_23 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_31 = sft_as_31 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_32 = sft_as_32 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_33 = sft_as_33 * np_orig_as.sum(axis=1)*1./len(l_orig)
    
        if (centr):
            orig = space_group.p_centr
            hh = numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2]))
            sft_11 = 0.5*(sft_as_11+sft_as_11.conjugate()*hh)
            sft_12 = 0.5*(sft_as_12+sft_as_12.conjugate()*hh)
            sft_13 = 0.5*(sft_as_13+sft_as_13.conjugate()*hh)
            sft_21 = 0.5*(sft_as_21+sft_as_21.conjugate()*hh)
            sft_22 = 0.5*(sft_as_22+sft_as_22.conjugate()*hh)
            sft_23 = 0.5*(sft_as_23+sft_as_23.conjugate()*hh)
            sft_31 = 0.5*(sft_as_31+sft_as_31.conjugate()*hh)
            sft_32 = 0.5*(sft_as_32+sft_as_32.conjugate()*hh)
            sft_33 = 0.5*(sft_as_33+sft_as_33.conjugate()*hh)
        else:
            sft_11, sft_12, sft_13 = sft_as_11, sft_as_12, sft_as_13
            sft_21, sft_22, sft_23 = sft_as_21, sft_as_22, sft_as_23
            sft_31, sft_32, sft_33 = sft_as_31, sft_as_32, sft_as_33            

        #sft_ij form the structure factor tensor in local coordinate system (ia, ib, ic)
        #chi in 10-12 cm; chim in muB (it is why here 0.2695)
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = self._orto_matrix(
                cell,
                sft_11*0.2695, sft_12*0.2695, sft_13*0.2695, 
                sft_21*0.2695, sft_22*0.2695, sft_23*0.2695, 
                sft_31*0.2695, sft_32*0.2695, sft_33*0.2695)

        refln = Refln()
        refln.h, refln.k, refln.l = copy.deepcopy(h), copy.deepcopy(k), copy.deepcopy(l)
        refln.f_nucl = f_nucl
        refln.sft_11, refln.sft_12, refln.sft_13 = s_11, s_12, s_13
        refln.sft_21, refln.sft_22, refln.sft_23 = s_21, s_22, s_23
        refln.sft_31, refln.sft_32, refln.sft_33 = s_31, s_32, s_33
        return refln


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

        X = B x, x = iB X
        xT*CHI*x = XT iBT CHI iB X
    
        output chiLOC = iBT CHI iB
        """
        m_ib = cell.m_ib
        ia, ib, ic = cell.ia, cell.ib, cell.ic
        """
        matrix_chi = numpy.array(
                [[self["chi_11"], self["chi_12"], self["chi_13"]],
                 [self["chi_12"], self["chi_22"], self["chi_23"]],
                 [self["chi_13"], self["chi_23"], self["chi_33"]]], 
                 dtype = float)
        #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
        #[a,b,c,alpha,beta,gamma]=ucp
        y1 = m_ib[0,0]
        y2 = m_ib[1,1]
        y3 = m_ib[2,2]
        y4 = m_ib[0,1]
        y5 = m_ib[0,2]
        y6 = m_ib[1,2]
        #B=[[x1,x4,x5],
        #   [0.,x2,x6],
        #   [0.,0.,x3]]
        #it shuld be checked
        #iB=numpy.linalg.inv(B)
        y1 = 1./x1
        y2 = 1./x2
        y3 = 1./x3
        y4 = -1*x4*1./(x1*x2)
        y6 = -1*x6*1./(x2*x3)
        y5 = (x4*x6-x2*x5)*1./(x1*x2*x3)
        """
        m_ib_norm = numpy.copy(m_ib)
        m_ib_norm[:,0] *= ia
        m_ib_norm[:,1] *= ib
        m_ib_norm[:,2] *= ic
        
        m_ibt_norm = m_ib_norm.transpose()
        
        r11, r12, r13 = m_ibt_norm[0, 0], m_ibt_norm[0, 1], m_ibt_norm[0, 2]
        r21, r22, r23 = m_ibt_norm[1, 0], m_ibt_norm[1, 1], m_ibt_norm[1, 2]
        r31, r32, r33 = m_ibt_norm[2, 0], m_ibt_norm[2, 1], m_ibt_norm[2, 2]
        
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = calc_mRmCmRT(
                r11, r12, r13, r21, r22, r23, r31, r32, r33,
                l_11, l_12, l_13, l_21, l_22, l_23, l_31, l_32, l_33)        
        """
        ibt_chi = numpy.matmul(m_ibt_norm, matrix_chi)
        matrix_chi_loc = numpy.matmul(ibt_chi, m_ib_norm)
        d_out = dict(matrix_chi_loc = matrix_chi_loc)
        self.update(d_out)
        """
        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33

    def apply_constraint(self):
        space_group = self.space_group
        cell = self.cell
        cell.apply_constraint()
        atom_site = self.atom_site
        atom_site_aniso = self.atom_site_aniso
        if atom_site_aniso is not None:
            atom_site_aniso.apply_space_group_constraint(atom_site, space_group)
        atom_site_magnetism_aniso = self.atom_site_magnetism_aniso
        if atom_site_magnetism_aniso is not None:
            atom_site_magnetism_aniso.apply_chi_iso_constraint(cell)
            atom_site_magnetism_aniso.apply_space_group_constraint(atom_site, space_group)
    
    def calc_hkl(self, sthol_min, sthovl_max):
        cell = self.cell
        space_group = self.space_group
        res = cell.calc_hkl(space_group, sthol_min, sthovl_max)
        return res

    def calc_hkl_in_range(self, sthol_min, sthovl_max):
        cell = self.cell
        res = cell.calc_hkl_in_range(sthol_min, sthovl_max)
        return res
