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

class AtomSiteAniso(ItemConstr):
    """
AtomSiteAniso
=================

Data items in the ATOM_SITE_ANISO category record details about
the atom sites in a crystal structure, such as atomic displacement parameters.
Data items in the ATOM_site_ANISO category record details about
magnetic properties of the atoms that occupy the atom sites.

Description in cif file:
--------------------------
::

 _atom_site_aniso_label    O1
 _atom_site_aniso_U_11     0.071(1)
 _atom_site_aniso_U_22     0.076(1)
 _atom_site_aniso_U_33     0.0342(9)
 _atom_site_aniso_U_12     0.008(1)
 _atom_site_aniso_U_13     0.0051(9) 
 _atom_site_aniso_U_23    -0.0030(9) 

Reference: 
-------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html>`_
    """    
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("u_11", "u_22", "u_33", "u_12", "u_13", "u_23", 
                          "b_11", "b_22", "b_33", "b_12", "b_13", "b_23")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_site_aniso"
    def __init__(self, label=None, 
                 u_11=None, u_22=None, u_33=None, u_12=None, u_13=None, u_23=None, 
                 b_11=None, b_22=None, b_33=None, b_12=None, b_13=None, b_23=None):

        super(AtomSiteAniso, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                            optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                            internal_attribute=self.INTERNAL_ATTRIBUTE,
                                            prefix=self.PREFIX)
        self.label = label
        self.u_11 = u_11
        self.u_12 = u_12
        self.u_13 = u_13
        self.u_22 = u_22
        self.u_23 = u_23
        self.u_33 = u_33
        self.b_11 = b_11
        self.b_12 = b_12
        self.b_13 = b_13
        self.b_22 = b_22
        self.b_23 = b_23
        self.b_33 = b_33
        if self.is_defined:
            self.form_object
        
    def __repr__(self):
        ls_out = ["AtomSiteAniso:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)        

    @property
    def label(self):
        """
Anisotropic atomic displacement parameters are looped in
a separate list. This code must match the
_atom_site_label of the associated atom in the atom coordinate
list and conform with the same rules described in
_atom_site_label.

Type: char

Reference: 
------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_aniso_label.html>`_
        """
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)

    @property
    def u_11(self):
        """
These are the standard anisotropic atomic displacement
components in angstroms squared which appear in the
structure-factor term

.. math::
 T = \\exp \\left[-2\\pi^2 \\Sigma_{i}\\Sigma_{j} \\left(U^{ij} h_{i} h_{j} a^{*}_{i}  a^{*}_{j} \\right) \\right]

- h is the Miller indices
- a* is the reciprocal-space cell lengths

The unique elements of the real symmetric matrix are
entered by row.

Default: 0.

Type: float

Reference:
--------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_aniso_U_.html>`_
        """
        return getattr(self, "__u_11")
    @u_11.setter
    def u_11(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_11", x_in)


    @property
    def u_22(self):
        """
see help for u_11 parameter
        """
        return getattr(self, "__u_22")
    @u_22.setter
    def u_22(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_22", x_in)

    @property
    def u_33(self):
        """
see help for u_11 parameter
        """
        return getattr(self, "__u_33")
    @u_33.setter
    def u_33(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_33", x_in)

    @property
    def u_12(self):
        """
see help for u_11 parameter
        """
        return getattr(self, "__u_12")
    @u_12.setter
    def u_12(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_12", x_in)

    @property
    def u_13(self):
        """
see help for u_11 parameter
        """
        return getattr(self, "__u_13")
    @u_13.setter
    def u_13(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_13", x_in)

    @property
    def u_23(self):
        """
see help for u_23 parameter
        """
        return getattr(self, "__u_23")
    @u_23.setter
    def u_23(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_23", x_in)




    @property
    def b_11(self):
        """
These are the standard anisotropic atomic displacement compo-
nents in Angströms squared which appear in the structure-factor
term 

.. math::
 T = \\exp \\left[-\frac{1}{4} \\Sigma_{i}\\Sigma_{j} \\left(B^{ij} h_{i} h_{j} a^{*}_{i}  a^{*}_{j} \\right) \\right]

where h = the Miller indices and a∗ = the reciprocal-space cell
lengths.

The unique elements of the real symmetric matrix are entered by
row.The IUCr Commission on Nomenclature recommends against
the use of B for reporting atomic displacement parameters. U,
being directly proportional to B, is preferred.


- h is the Miller indices
- a* is the reciprocal-space cell lengths

The unique elements of the real symmetric matrix are
entered by row.

Default: 0.

Type: float

Reference:
--------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_aniso_B_.html>`_
        """
        return getattr(self, "__b_11")
    @b_11.setter
    def b_11(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_11", x_in)


    @property
    def b_22(self):
        """
see help for b_11 parameter
        """
        return getattr(self, "__b_22")
    @b_22.setter
    def b_22(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_22", x_in)

    @property
    def b_33(self):
        """
see help for b_11 parameter
        """
        return getattr(self, "__b_33")
    @b_33.setter
    def b_33(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_33", x_in)

    @property
    def b_12(self):
        """
see help for b_11 parameter
        """
        return getattr(self, "__b_12")
    @b_12.setter
    def b_12(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_12", x_in)

    @property
    def b_13(self):
        """
see help for b_11 parameter
        """
        return getattr(self, "__b_13")
    @b_13.setter
    def b_13(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_13", x_in)

    @property
    def b_23(self):
        """
see help for b_11 parameter
        """
        return getattr(self, "__b_23")
    @b_23.setter
    def b_23(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_23", x_in)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)



    @property
    def is_variable(self):
        res = any([self.u_11.refinement,
                   self.u_22.refinement,
                   self.u_33.refinement,
                   self.u_12.refinement,
                   self.u_13.refinement,
                   self.u_23.refinement])
        return res


    def get_variables(self):
        l_variable = []
        if self.u_11.refinement: l_variable.append(self.u_11)
        if self.u_22.refinement: l_variable.append(self.u_22)
        if self.u_33.refinement: l_variable.append(self.u_33)
        if self.u_12.refinement: l_variable.append(self.u_12)
        if self.u_13.refinement: l_variable.append(self.u_13)
        if self.u_23.refinement: l_variable.append(self.u_23)
        return l_variable

    def apply_space_group_constraint(self, atom_site, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        l_numb = atom_site.calc_constr_number(space_group)
        label_aniso = self.label
        label = atom_site.label
        l_ind = [label.index(_1) for _1 in label_aniso]
        for index, u_11, u_22, u_33, u_12, u_13, u_23 in zip(l_ind, self.u_11, self.u_22, self.u_33, 
                                                                    self.u_12, self.u_13, self.u_23):
            u_11.constraint_flag, u_22.constraint_flag, u_33.constraint_flag = False, False, False
            u_12.constraint_flag, u_13.constraint_flag, u_23.constraint_flag = False, False, False
            numb = l_numb[index]
            if numb == 1:
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 2:
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
            elif numb == 3:
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
            elif numb == 4:
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 5:
                u_22.value = u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 6:
                u_22.value = u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_23.value = u_13.value 
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 7:
                u_22.value = u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_23.value = -1.*u_13.value 
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 8:
                u_22.value = u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 9:
                u_33.value = u_22.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
            elif numb == 10:
                u_33.value = u_22.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_13.value = 1.*u_12.value 
                u_13.refinement = False
                u_13.constraint_flag = True
            elif numb == 11:
                u_33.value = u_22.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_13.value = -1.*u_12.value 
                u_13.refinement = False
                u_13.constraint_flag = True
            elif numb == 12:
                u_33.value = u_22.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 13:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_12.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 14:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 15:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_12.constraint_flag = True
                u_23.value = 2.*u_13.value
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 16:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_12.value = 0.5*u_11.value
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 17:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_33.value = 1.0*u_11.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_12.value = 0.
                u_12.refinement = False
                u_12.constraint_flag = True
                u_13.value = 0.
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 0.
                u_23.refinement = False
                u_23.constraint_flag = True
            elif numb == 18:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_22.constraint_flag = True
                u_33.value = 1.0*u_11.value
                u_33.refinement = False
                u_33.constraint_flag = True
                u_13.value = 1.0*u_12.value
                u_13.refinement = False
                u_13.constraint_flag = True
                u_23.value = 1.0*u_12.value
                u_23.refinement = False
                u_23.constraint_flag = True


    def calc_beta(self, cell):
        """
        calculate beta_ij from U_ij
        """
        ia, ib, ic = cell.reciprocal_length_a, cell.reciprocal_length_b, cell.reciprocal_length_c
        u_11, u_22, u_33 = float(self.u_11), float(self.u_22), float(self.u_33)
        u_12, u_13, u_23 = float(self.u_12), float(self.u_13), float(self.u_23)
        beta_11 = 2.*numpy.pi**2*u_11*ia**2
        beta_22 = 2.*numpy.pi**2*u_22*ib**2
        beta_33 = 2.*numpy.pi**2*u_33*ic**2
        beta_12 = 2.*numpy.pi**2*u_12*ia*ib
        beta_13 = 2.*numpy.pi**2*u_13*ia*ic
        beta_23 = 2.*numpy.pi**2*u_23*ib*ic
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23


    def _form_adp(self, atom_site):
        label = numpy.array(atom_site.label, dtype=str)
        b_iso = numpy.array(atom_site.b_iso, dtype=float)
        adp_type = numpy.array(atom_site.adp_type, dtype=str)

        b_iso_in = numpy.zeros(b_iso.shape, dtype=float)
        b_iso_in = numpy.where(adp_type=="uiso", b_iso, b_iso_in)
        label_aniso = numpy.array(self.label, dtype=str)
        if not(set(label_aniso).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False
        
        u_11 = numpy.array(self.u_11, dtype=float)
        u_22 = numpy.array(self.u_22, dtype=float)
        u_33 = numpy.array(self.u_33, dtype=float)
        u_12 = numpy.array(self.u_12, dtype=float)
        u_13 = numpy.array(self.u_13, dtype=float)
        u_23 = numpy.array(self.u_23, dtype=float)


        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in label_aniso], dtype=int)

        u_11_in = numpy.zeros(label.shape, dtype=float)
        u_22_in = numpy.zeros(label.shape, dtype=float)
        u_33_in = numpy.zeros(label.shape, dtype=float)
        u_12_in = numpy.zeros(label.shape, dtype=float)
        u_13_in = numpy.zeros(label.shape, dtype=float)
        u_23_in = numpy.zeros(label.shape, dtype=float)

        u_11_in[np_index], u_22_in[np_index], u_33_in[np_index] = u_11, u_22, u_33
        u_12_in[np_index], u_13_in[np_index], u_23_in[np_index] = u_12, u_13, u_23

        np_zeros = numpy.zeros(label.shape, dtype=float)
        u_11_in = numpy.where(adp_type=="uani", u_11_in, np_zeros)
        u_22_in = numpy.where(adp_type=="uani", u_22_in, np_zeros)
        u_33_in = numpy.where(adp_type=="uani", u_33_in, np_zeros)
        u_12_in = numpy.where(adp_type=="uani", u_12_in, np_zeros)
        u_13_in = numpy.where(adp_type=="uani", u_13_in, np_zeros)
        u_23_in = numpy.where(adp_type=="uani", u_23_in, np_zeros)

        adp = None
        #adp = ADP(u_11=u_11_in, u_22=u_22_in, u_33=u_33_in, 
        #          u_12=u_12_in, u_13=u_13_in, u_23=u_23_in, 
        #          b_iso=b_iso_in)
        return adp



class AtomSiteAnisoL(LoopConstr):
    """
AtomSiteAnisoL
=================

Data items in the ATOM_SITE_ANISO category record details about
the atom sites in a crystal structure, such as atomic displacement parameters.
Data items in the ATOM_site_ANISO category record details about
magnetic properties of the atoms that occupy the atom sites.

Description in cif file:
--------------------------
::

 loop_
 _atom_site_aniso_label
 _atom_site_aniso_U_11
 _atom_site_aniso_U_22
 _atom_site_aniso_U_33
 _atom_site_aniso_U_12
 _atom_site_aniso_U_13
 _atom_site_aniso_U_23
  O1   .071(1) .076(1) .0342(9) .008(1)   .0051(9) -.0030(9) 
  C2   .060(2) .072(2) .047(1)  .002(2)   .013(1)  -.009(1)  


Reference: 
-------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html>`_
    """    
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomSiteAniso
    def __init__(self, item=[], loop_name=""):
        super(AtomSiteAnisoL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomSiteAnisoL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
