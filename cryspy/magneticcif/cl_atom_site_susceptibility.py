"""
atom_site_susceptibility_moment was added according to suggestion of Henrik Thoma <h.thoma@fz-juelich.de>:

Since a large number of orientations and magnetic fields measured for one compound in the weak ferromagnetic state, 
it is introduced an additional weak ferromagnetic moment tensor to refine data. 

atom_moment_i = (moment_ij + abs(field) * chi_ij) * field_j
"""

__author__ = 'ikibalin'
__version__ = "2019_12_05"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


class AtomSiteSusceptibility(ItemConstr):
    """
Data items in the ATOM_SITE_MAGNETISM_ANISO category record details about
magnetic properties of the atoms that occupy the atom sites.

Description in cif file::
 
 _atom_site_susceptibility_label     Fe3A 
 _atom_site_susceptibility_chi_type  Cani 
 _atom_site_susceptibility_chi_11   -3.468(74) 
 _atom_site_susceptibility_chi_22   -3.468(74) 
 _atom_site_susceptibility_chi_33   -3.468(74) 
 _atom_site_susceptibility_chi_12    0.0
 _atom_site_susceptibility_chi_13    0.0
 _atom_site_susceptibility_chi_23    0.0 
    """    
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("chi_type", "moment_type", "chi_11", "chi_22", "chi_33", "chi_12", "chi_13", "chi_23", 
                          "moment_11", "moment_22", "moment_33", "moment_12", "moment_13", "moment_23")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_site_susceptibility"
    ACCESIBLE_CHI_TYPE = ("Ciso", "Cani")
    ACCESIBLE_MOMENT_TYPE = ("Miso", "Mani")
    def __init__(self, label=None, chi_type=None, 
                 chi_11=None, chi_22=None, chi_33=None,
                 chi_12=None, chi_13=None, chi_23=None,
                 moment_type=None, 
                 moment_11=None, moment_22=None, moment_33=None,
                 moment_12=None, moment_13=None, moment_23=None,
                 ):
        super(AtomSiteSusceptibility, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                     optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                     internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                     prefix=self.PREFIX)
        self.label = label
        self.chi_type = chi_type
        self.chi_11 = chi_11
        self.chi_12 = chi_12
        self.chi_13 = chi_13
        self.chi_22 = chi_22
        self.chi_23 = chi_23
        self.chi_33 = chi_33
        self.moment_type = moment_type
        self.moment_11 = moment_11
        self.moment_12 = moment_12
        self.moment_13 = moment_13
        self.moment_22 = moment_22
        self.moment_23 = moment_23
        self.moment_33 = moment_33
        if self.is_defined:
            self.form_object
        
    def __repr__(self):
        ls_out = ["AtomSiteSusceptibility:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)


    @property
    def label(self):
        """
The _atom_site_magnetism_aniso_label is a unique identifier for a particular site
in the crystal. 

Type: char
        """
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)



    @property
    def chi_type(self):
        """
Chi type. 

The data value must be one of the following::
 
 Ciso, Cani
        """
        return getattr(self, "__chi_type")
    @chi_type.setter
    def chi_type(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_CHI_TYPE):
                warnings.warn(f"chi_type '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__chi_type", x_in)

    @property
    def chi_11(self):
        """
chi_11
Default: 0.
Type: float
        """
        return getattr(self, "__chi_11")
    @chi_11.setter
    def chi_11(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_11", x_in)


    @property
    def chi_22(self):
        """
chi_22
Default: 0.
Type: float
        """
        return getattr(self, "__chi_22")
    @chi_22.setter
    def chi_22(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_22", x_in)

    @property
    def chi_33(self):
        """
chi_33
Default: 0.
Type: float
        """
        return getattr(self, "__chi_33")
    @chi_33.setter
    def chi_33(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_33", x_in)


    @property
    def chi_12(self):
        """
chi_12
Default: 0.
Type: float
        """
        return getattr(self, "__chi_12")
    @chi_12.setter
    def chi_12(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_12", x_in)

    @property
    def chi_13(self):
        """
chi_13
Default: 0.
Type: float
        """
        return getattr(self, "__chi_13")
    @chi_13.setter
    def chi_13(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_13", x_in)

    @property
    def chi_23(self):
        """
chi_23
Default: 0.
Type: float
        """
        return getattr(self, "__chi_23")
    @chi_23.setter
    def chi_23(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__chi_23", x_in)





    @property
    def moment_type(self):
        """
Chi type. 

The data value must be one of the following::
 
 Miso, Mani
        """
        return getattr(self, "__moment_type")
    @moment_type.setter
    def moment_type(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_MOMENT_TYPE):
                warnings.warn(f"moment_type '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__moment_type", x_in)

    @property
    def moment_11(self):
        """
moment_11
Default: 0.
Type: float
        """
        return getattr(self, "__moment_11")
    @moment_11.setter
    def moment_11(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_11", x_in)


    @property
    def moment_22(self):
        """
moment_22
Default: 0.
Type: float
        """
        return getattr(self, "__moment_22")
    @moment_22.setter
    def moment_22(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_22", x_in)

    @property
    def moment_33(self):
        """
moment_33
Default: 0.
Type: float
        """
        return getattr(self, "__moment_33")
    @moment_33.setter
    def moment_33(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_33", x_in)


    @property
    def moment_12(self):
        """
moment_12
Default: 0.
Type: float
        """
        return getattr(self, "__moment_12")
    @moment_12.setter
    def moment_12(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_12", x_in)

    @property
    def moment_13(self):
        """
moment_13
Default: 0.
Type: float
        """
        return getattr(self, "__moment_13")
    @moment_13.setter
    def moment_13(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_13", x_in)

    @property
    def moment_23(self):
        """
moment_23
Default: 0.
Type: float
        """
        return getattr(self, "__moment_23")
    @moment_23.setter
    def moment_23(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__moment_23", x_in)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)



    @property
    def is_variable(self):
        _l = []
        if self.chi_11 is not None: _l.append(self.chi_11.refinement)
        if self.chi_22 is not None: _l.append(self.chi_22.refinement)
        if self.chi_33 is not None: _l.append(self.chi_33.refinement)
        if self.chi_12 is not None: _l.append(self.chi_12.refinement)
        if self.chi_13 is not None: _l.append(self.chi_13.refinement)
        if self.chi_23 is not None: _l.append(self.chi_23.refinement)

        if self.moment_11 is not None: _l.append(self.moment_11.refinement)
        if self.moment_22 is not None: _l.append(self.moment_22.refinement)
        if self.moment_33 is not None: _l.append(self.moment_33.refinement)
        if self.moment_12 is not None: _l.append(self.moment_12.refinement)
        if self.moment_13 is not None: _l.append(self.moment_13.refinement)
        if self.moment_23 is not None: _l.append(self.moment_23.refinement)

        res = any(_l)
        return res


    def get_variables(self):
        l_variable = []
        if self.chi_11 is not None:
            if self.chi_11.refinement: l_variable.append(self.chi_11)
        if self.chi_22 is not None:
            if self.chi_22.refinement: l_variable.append(self.chi_22)
        if self.chi_33 is not None:
            if self.chi_33.refinement: l_variable.append(self.chi_33)
        if self.chi_12 is not None:
            if self.chi_12.refinement: l_variable.append(self.chi_12)
        if self.chi_13 is not None:
            if self.chi_13.refinement: l_variable.append(self.chi_13)
        if self.chi_23 is not None:
            if self.chi_23.refinement: l_variable.append(self.chi_23)

        if self.moment_11 is not None:
            if self.moment_11.refinement: l_variable.append(self.moment_11)
        if self.moment_22 is not None:
            if self.moment_22.refinement: l_variable.append(self.moment_22)
        if self.moment_33 is not None:
            if self.moment_33.refinement: l_variable.append(self.moment_33)
        if self.moment_12 is not None:
            if self.moment_12.refinement: l_variable.append(self.moment_12)
        if self.moment_13 is not None:
            if self.moment_13.refinement: l_variable.append(self.moment_13)
        if self.moment_23 is not None:
            if self.moment_23.refinement: l_variable.append(self.moment_23)
        return l_variable


    def _form_magnetism(self, atom_site, atom_magnetism):
        label = numpy.array(atom_site.label, dtype=str)
        label_aniso = numpy.array(self.label, dtype=str)
        if not(set(label_aniso).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False
        
        j0_A_in, j2_A_in = atom_site.j0_A, atom_site.j2_A
        j0_a_in, j2_a_in = atom_site.j0_a, atom_site.j2_a
        j0_B_in, j2_B_in = atom_site.j0_B, atom_site.j2_B
        j0_b_in, j2_b_in = atom_site.j0_b, atom_site.j2_b
        j0_C_in, j2_C_in = atom_site.j0_C, atom_site.j2_C
        j0_c_in, j2_c_in = atom_site.j0_c, atom_site.j2_c
        j0_D_in, j2_D_in = atom_site.j0_D, atom_site.j2_D

        lande_in, kappa_in = atom_magnetism._form_lande_kappa(atom_site)
        
        #chi_type = numpy.array(self.chi_type, dtype=str)
        chi_11 = numpy.array(self.chi_11, dtype=float)
        chi_22 = numpy.array(self.chi_22, dtype=float)
        chi_33 = numpy.array(self.chi_33, dtype=float)
        chi_12 = numpy.array(self.chi_12, dtype=float)
        chi_13 = numpy.array(self.chi_13, dtype=float)
        chi_23 = numpy.array(self.chi_23, dtype=float)

        moment_11 = numpy.array(self.moment_11, dtype=float)
        moment_22 = numpy.array(self.moment_22, dtype=float)
        moment_33 = numpy.array(self.moment_33, dtype=float)
        moment_12 = numpy.array(self.moment_12, dtype=float)
        moment_13 = numpy.array(self.moment_13, dtype=float)
        moment_23 = numpy.array(self.moment_23, dtype=float)


        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in label_aniso], dtype=int)

        chi_11_in = numpy.zeros(label.shape, dtype=float)
        chi_22_in = numpy.zeros(label.shape, dtype=float)
        chi_33_in = numpy.zeros(label.shape, dtype=float)
        chi_12_in = numpy.zeros(label.shape, dtype=float)
        chi_13_in = numpy.zeros(label.shape, dtype=float)
        chi_23_in = numpy.zeros(label.shape, dtype=float)

        moment_11_in = numpy.zeros(label.shape, dtype=float)
        moment_22_in = numpy.zeros(label.shape, dtype=float)
        moment_33_in = numpy.zeros(label.shape, dtype=float)
        moment_12_in = numpy.zeros(label.shape, dtype=float)
        moment_13_in = numpy.zeros(label.shape, dtype=float)
        moment_23_in = numpy.zeros(label.shape, dtype=float)

        chi_11_in[np_index], chi_22_in[np_index], chi_33_in[np_index] = chi_11, chi_22, chi_33
        chi_12_in[np_index], chi_13_in[np_index], chi_23_in[np_index] = chi_12, chi_13, chi_23

        moment_11_in[np_index], moment_22_in[np_index], moment_33_in[np_index] = moment_11, moment_22, moment_33
        moment_12_in[np_index], moment_13_in[np_index], moment_23_in[np_index] = moment_12, moment_13, moment_23

        magnetism = None
        #magnetism = Magnetism(factor_lande=lande_in, kappa=kappa_in, 
        #                      chi_11=chi_11_in, chi_22=chi_22_in, chi_33=chi_33_in,
        #                      chi_12=chi_12_in, chi_13=chi_13_in, chi_23=chi_23_in,
        #                      moment_11=moment_11_in, moment_22=moment_22_in, moment_33=moment_33_in,
        #                      moment_12=moment_12_in, moment_13=moment_13_in, moment_23=moment_23_in,
        #                      j0_A=j0_A_in, j0_a=j0_a_in, j0_B=j0_B_in, j0_b=j0_b_in,
        #                      j0_C=j0_C_in, j0_c=j0_c_in, j0_D=j0_D_in,
        #                      j2_A=j2_A_in, j2_a=j2_a_in, j2_B=j2_B_in, j2_b=j2_b_in,
        #                      j2_C=j2_C_in, j2_c=j2_c_in, j2_D=j2_D_in)

        return magnetism

    def apply_space_group_constraint(self, atom_site, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        flag = True
        l_numb = atom_site.calc_constr_number(space_group)
        label_aniso = self.label
        label = atom_site.label
        index = label.index(label_aniso)
        chi_11, chi_22, chi_33, chi_12, chi_13, chi_23 = self.chi_11, self.chi_22, self.chi_33, self.chi_12, self.chi_13, self.chi_23
        moment_11, moment_22, moment_33, moment_12, moment_13, moment_23 = self.moment_11, self.moment_22, self.moment_33, self.moment_12, self.moment_13, self.moment_23
        flag_chi = (self.chi_type is not None)
        if flag_chi:
            flag_chi = self.chi_type.lower().startswith("cani")
        flag_moment = (self.moment_type is not None)
        if flag_moment:
            flag_moment = self.moment_type.lower().startswith("mani")
        if flag_chi:
            chi_11.constraint_flag, chi_22.constraint_flag, chi_33.constraint_flag = False, False, False
            chi_12.constraint_flag, chi_13.constraint_flag, chi_23.constraint_flag = False, False, False
        if flag_moment:
            moment_11.constraint_flag, moment_22.constraint_flag, moment_33.constraint_flag = False, False, False
            moment_12.constraint_flag, moment_13.constraint_flag, moment_23.constraint_flag = False, False, False
        numb = l_numb[index]
        if numb == 1:
            if flag_chi:
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 2:
            if flag_chi:
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
            if flag_moment:
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
        elif numb == 3:
            if flag_chi:
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
        elif numb == 4:
            if flag_chi:
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 5:
            if flag_chi:
                chi_22.value = chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 6:
            if flag_chi:
                chi_22.value = chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_23.value = chi_13.value 
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_23.value = moment_13.value 
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 7:
            if flag_chi:
                chi_22.value = chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_23.value = -1.*chi_13.value 
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_23.value = -1.*moment_13.value 
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 8:
            if flag_chi:
                chi_22.value = chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 9:
            if flag_chi:
                chi_33.value = chi_22.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
            if flag_moment:
                moment_33.value = moment_22.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
        elif numb == 10:
            if flag_chi:
                chi_33.value = chi_22.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_13.value = 1.*chi_12.value 
                chi_13.refinement = False
                chi_13.constraint_flag = True
            if flag_moment:
                moment_33.value = moment_22.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_13.value = 1.*moment_12.value 
                moment_13.refinement = False
                moment_13.constraint_flag = True
        elif numb == 11:
            if flag_chi:
                chi_33.value = chi_22.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_13.value = -1.*chi_12.value 
                chi_13.refinement = False
                chi_13.constraint_flag = True
            if flag_moment:
                moment_33.value = moment_22.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_13.value = -1.*moment_12.value 
                moment_13.refinement = False
                moment_13.constraint_flag = True
        elif numb == 12:
            if flag_chi:
                chi_33.value = chi_22.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_33.value = moment_22.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 13:
            if flag_chi:
                chi_12.value = 0.5*chi_22.value
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.5*moment_22.value
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 14:
            if flag_chi:
                chi_12.value = 0.5*chi_22.value
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.5*moment_22.value
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 15:
            if flag_chi:
                chi_12.value = 0.5*chi_22.value
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_23.value = 2.*chi_13.value
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_12.value = 0.5*moment_22.value
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_23.value = 2.*moment_13.value
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 16:
            if flag_chi:
                chi_22.value = 1.0*chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_12.value = 0.5*chi_11.value
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = 1.0*moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_12.value = 0.5*moment_11.value
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 17:
            if flag_chi:
                chi_22.value = 1.0*chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_33.value = 1.0*chi_11.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_12.value = 0.
                chi_12.refinement = False
                chi_12.constraint_flag = True
                chi_13.value = 0.
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 0.
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = 1.0*moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_33.value = 1.0*moment_11.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_12.value = 0.
                moment_12.refinement = False
                moment_12.constraint_flag = True
                moment_13.value = 0.
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 0.
                moment_23.refinement = False
                moment_23.constraint_flag = True
        elif numb == 18:
            if flag_chi:
                chi_22.value = 1.0*chi_11.value
                chi_22.refinement = False
                chi_22.constraint_flag = True
                chi_33.value = 1.0*chi_11.value
                chi_33.refinement = False
                chi_33.constraint_flag = True
                chi_13.value = 1.0*chi_12.value
                chi_13.refinement = False
                chi_13.constraint_flag = True
                chi_23.value = 1.0*chi_12.value
                chi_23.refinement = False
                chi_23.constraint_flag = True
            if flag_moment:
                moment_22.value = 1.0*moment_11.value
                moment_22.refinement = False
                moment_22.constraint_flag = True
                moment_33.value = 1.0*moment_11.value
                moment_33.refinement = False
                moment_33.constraint_flag = True
                moment_13.value = 1.0*moment_12.value
                moment_13.refinement = False
                moment_13.constraint_flag = True
                moment_23.value = 1.0*moment_12.value
                moment_23.refinement = False
                moment_23.constraint_flag = True
        return flag


    def apply_chi_iso_constraint(self, cell):
        flag = True
        c_a = cell.cos_a
        s_ib = cell.sin_ib
        s_ig = cell.sin_ig
        c_ib = cell.cos_ib
        c_ig = cell.cos_ig
        #not sure, it is better to check
        chi_type = self.chi_type
        chi_11, chi_22, chi_33 = self.chi_11, self.chi_22, self.chi_33
        chi_12, chi_13, chi_23 = self.chi_12, self.chi_13, self.chi_23
        if chi_type is None:
            return flag
        if chi_type.lower().startswith("ciso"):
            chi_22.value = chi_11.value
            chi_33.value = chi_11.value
            chi_12.value = chi_11.value*c_ig
            chi_13.value = chi_11.value*c_ib
            chi_23.value = chi_11.value*(c_ib*c_ig-s_ib*s_ig*c_a)
            chi_22.refinement, chi_33.refinement, chi_12.refinement, chi_13.refinement, chi_23.refinement =False, False, False, False, False
        return  flag

    def apply_moment_iso_constraint(self, cell):
        flag = True
        c_a = cell.cos_a
        s_ib = cell.sin_ib
        s_ig = cell.sin_ig
        c_ib = cell.cos_ib
        c_ig = cell.cos_ig
        #not sure, it is better to check
        moment_type = self.moment_type
        moment_11, moment_22, moment_33 =  self.moment_11, self.moment_22, self.moment_33 
        moment_12, moment_13, moment_23 =  self.moment_12, self.moment_13, self.moment_23
        if moment_type is  None:
            return flag
        if moment_type.lower().startswith("miso"):
            moment_22.value = moment_11.value
            moment_33.value = moment_11.value
            moment_12.value = moment_11.value*c_ig
            moment_13.value = moment_11.value*c_ib
            moment_23.value = moment_11.value*(c_ib*c_ig-s_ib*s_ig*c_a)
            moment_22.refinement, moment_33.refinement, moment_12.refinement, moment_13.refinement, moment_23.refinement =False, False, False, False, False
        return  flag



class AtomSiteSusceptibilityL(LoopConstr):
    """
Data items in the ATOM_SITE_MAGNETISM_ANISO category record details about
magnetic properties of the atoms that occupy the atom sites.

Description in cif file::

 loop_  
 _atom_site_susceptibility_label
 _atom_site_susceptibility_chi_type
 _atom_site_susceptibility_chi_11
 _atom_site_susceptibility_chi_12
 _atom_site_susceptibility_chi_13
 _atom_site_susceptibility_chi_22
 _atom_site_susceptibility_chi_23
 _atom_site_susceptibility_chi_33
 _atom_site_susceptibility_moment_type
 _atom_site_susceptibility_moment_11
 _atom_site_susceptibility_moment_12
 _atom_site_susceptibility_moment_13
 _atom_site_susceptibility_moment_22
 _atom_site_susceptibility_moment_23
 _atom_site_susceptibility_moment_33
  Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468 Mani 0. 0. 0. 0. 0. 0.
  Fe3B Cani 3.041      0.0 0.0  3.041 0.0  3.041 Mani 0. 0. 0. 0. 0. 0.
    
    """    
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomSiteSusceptibility
    def __init__(self, item=[], loop_name=""):
        super(AtomSiteSusceptibilityL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomSiteSusceptibilityL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)


    def apply_space_group_constraint(self, atom_site, space_group):
        l_flag = [_item.apply_space_group_constraint(atom_site, space_group) for _item in self.item]
        return all(l_flag)
    def apply_chi_iso_constraint(self, cell):
        l_flag = [_item.apply_chi_iso_constraint(cell) for _item in self.item]
        return all(l_flag)

    def apply_moment_iso_constraint(self, cell):
        l_flag = [_item.apply_moment_iso_constraint(cell) for _item in self.item]
        return all(l_flag)
