__author__ = 'ikibalin'
__version__ = "2019_12_03"
import os
import math
import numpy
import warnings
from typing import List, Tuple

from pycifstar import Global

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


import cryspy.corecif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS
from cryspy.symcif.cl_space_group import SpaceGroup



class Cell(ItemConstr):
    """
Data items in the Cell class record details about the
crystallographic cell parameters and their measurement.

Description in cif file::

    _cell_length_a                     5.959(1)
    _cell_length_b                     14.956(1)
    _cell_length_c                     19.737(3)
    _cell_angle_alpha                  90
    _cell_angle_beta                   90
    _cell_angle_gamma                  90

:FIXME: the following attributes are not introduced: 
        measurement_pressure, measurement_radiation, measurement_reflns_used,
        measurement_temperature, measurement_theta_max, measurement_theta_min, 
        measurement_wavelength, special_details.
    """
    MANDATORY_ATTRIBUTE = ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ("m_m", "m_b", "m_ib", "m_ib_norm", "type_cell",
    "cos_a", "cos_b", "cos_g", "sin_a", "sin_b", "sin_g",
    "cos_a_sq", "cos_b_sq", "cos_g_sq", "sin_a_sq", "sin_b_sq", "sin_g_sq",
    "cos_ia", "cos_ib", "cos_ig", "sin_ia", "sin_ib", "sin_ig",
    "cos_ia_sq", "cos_ib_sq", "cos_ig_sq", "sin_ia_sq", "sin_ib_sq", "sin_ig_sq",
    "reciprocal_length_a", "reciprocal_length_b", "reciprocal_length_c", 
                          "reciprocal_angle_alpha", "reciprocal_angle_beta", "reciprocal_angle_gamma",
                          "volume", "formula_units_z", "it_coordinate_system_code")
    PREFIX = "cell"
    def __init__(self, length_a=None, length_b=None, length_c=None, angle_alpha=None, angle_beta=None, 
    angle_gamma=None, formula_units_z=None):
        super(Cell, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                   optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                   internal_attribute=self.INTERNAL_ATTRIBUTE,
                                   prefix=self.PREFIX)
        self.length_a = length_a
        self.length_b = length_b
        self.length_c = length_c
        self.angle_alpha = angle_alpha
        self.angle_beta = angle_beta
        self.angle_gamma = angle_gamma
        self.formula_units_z = formula_units_z

        if self.is_defined:
            self.form_object
        
    @property
    def length_a(self):
        """
Unit-cell lengths in angstroms corresponding to the structure
reported. 

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html>`_
        """
        return getattr(self, "__length_a")
    @length_a.setter
    def length_a(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__length_a", x_in)

    @property
    def length_b(self):
        """
Unit-cell lengths in angstroms corresponding to the structure
reported. 
The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html>`_
        """
        return getattr(self, "__length_b")
    @length_b.setter
    def length_b(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__length_b", x_in)

    @property
    def length_c(self):
        """
Unit-cell lengths in angstroms corresponding to the structure
reported. 

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html>`_
        """
        return getattr(self, "__length_c")
    @length_c.setter
    def length_c(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__length_c", x_in)

    @property
    def angle_alpha(self):
        """
Unit-cell angles of the reported structure in degrees.
The permitted range is 0.0 -> 180.0 

Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html>`_
        """
        return getattr(self, "__angle_alpha")
    @angle_alpha.setter
    def angle_alpha(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__angle_alpha", x_in)

    @property
    def angle_beta(self):
        """
Unit-cell angles of the reported structure in degrees.
The permitted range is 0.0 -> 180.0 
Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html>`_
        """
        return getattr(self, "__angle_beta")
    @angle_beta.setter
    def angle_beta(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__angle_beta", x_in)

    @property
    def angle_gamma(self):
        """
Unit-cell angles of the reported structure in degrees.
The permitted range is 0.0 -> 180.0 
Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html>`_
        """
        return getattr(self, "__angle_gamma")
    @angle_gamma.setter
    def angle_gamma(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__angle_gamma", x_in)

    @property
    def formula_units_z(self):
        """
The number of the formula units in the unit cell as specified
by _chemical_formula_structural, _chemical_formula_moiety or
_chemical_formula_sum.

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_formula_units_Z.html>`_
        """
        return getattr(self, "__formula_units_z")
    @formula_units_z.setter
    def formula_units_z(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__formula_units_z", x_in)



    @property
    def volume(self):
        """
Cell volume V in angstroms cubed.

.. math::

    V = a \\cdot b \\cdot c \\cdot \\sqrt{
        1 - \\cos^{2} \\alpha - \\cos^{2} \\beta  - \\cos^{2} \\gamma 
        + 2 \\cos \\alpha  \\cos \\beta  \\cos \\gamma  }

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_volume.html>`_
        """
        return getattr(self, "__volume")


    @property
    def reciprocal_length_a(self):
        """
The reciprocal-cell lengths in inverse angstroms.  These are
related to the real cell by:

.. math::

  a^{*} = b \\cdot c \\cdot \\sin(\\alpha) / V

where V is the cell volume.

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html>`_
        """
        return getattr(self, "__reciprocal_length_a")
    @property
    def reciprocal_length_b(self):
        """
The reciprocal-cell lengths in inverse angstroms.  These are
related to the real cell by:

.. math::

  b^{*} = c \\cdot a \\cdot \\sin(\\beta) / V

where V is the cell volume.

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html>`_
        """
        return getattr(self, "__reciprocal_length_b")
    @property
    def reciprocal_length_c(self):
        """
The reciprocal-cell lengths in inverse angstroms.  These are
related to the real cell by:

.. math::

  c^{*} = a \\cdot b \\cdot \\sin(\\gamma) / V

where V is the cell volume.

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> infinity

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html>`_
        """
        return getattr(self, "__reciprocal_length_c")


    @property
    def reciprocal_angle_alpha(self):
        """
The angles defining the reciprocal cell in degrees. These
are related to those in the real cell by:

.. math::

  \\cos(\\alpha^{*}) = \\frac{\\cos(\\beta)\\cos(\\gamma) - \\cos(\\alpha)}{\\sin(\\beta)\\sin(\\gamma)}

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> 180.0 
Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html>`_
        """
        return getattr(self, "__reciprocal_angle_alpha")


    @property
    def reciprocal_angle_beta(self):
        """
The angles defining the reciprocal cell in degrees. These
are related to those in the real cell by:

.. math::

  \\cos(\\beta^{*})  = \\frac{\\cos(\\gamma)\\cos(\\alpha) - \\cos(\\beta)}{\\sin(\\gamma)\\sin(\\alpha)}

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> 180.0 
Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html>`_
        """
        return getattr(self, "__reciprocal_angle_beta")

    @property
    def reciprocal_angle_gamma(self):
        """
The angles defining the reciprocal cell in degrees. These
are related to those in the real cell by:

.. math::

  \\cos(\\gamma^{*}) = \\frac{\\cos(\\alpha)\\cos(\\beta) - \\cos(\\gamma)}{\\sin(\\alpha)\\sin(\\beta)}

Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
     New York: John Wiley & Sons Inc.

The permitted range is 0.0 -> 180.0 
Enumeration default: 90.0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html>`_
        """
        return getattr(self, "__reciprocal_angle_gamma")

    @property
    def m_b(self):
        """
B matrix defined 
in Carthezian coordinate system where :math:`X` along axis :math:`a^{*}`,
:math:`Z` alond axis :math:`c` and :math:`Y = Z \\times X`

        """
        return getattr(self, "__m_b")
    @property
    def m_ib(self):
        """
Inversed B matrix
        """
        return getattr(self, "__m_ib")
    @property
    def m_ib_norm(self):
        """
Normalized inversed B matrix
        """
        return getattr(self, "__m_ib_norm")
    @property
    def m_m(self):
        """
M matrix defined 
in Carthezian coordinate system where :math:`X` along axis :math:`a`,
:math:`Z` alond axis :math:`c^{*}` and :math:`Y = Z \\times X`

        """
        return getattr(self, "__m_m")

    @property
    def type_cell(self):
        return getattr(self, "__type_cell")

    @property
    def it_coordinate_system_code(self):
        return getattr(self, "__it_coordinate_system_code")

    @property
    def cos_a(self):
        return getattr(self, "__cos_a")
    @property
    def cos_b(self):
        return getattr(self, "__cos_b")
    @property
    def cos_g(self):
        return getattr(self, "__cos_g")
    @property
    def sin_a(self):
        return getattr(self, "__sin_a")
    @property
    def sin_b(self):
        return getattr(self, "__sin_b")
    @property
    def sin_g(self):
        return getattr(self, "__sin_g")
    @property
    def cos_a_sq(self):
        return getattr(self, "__cos_a_sq")
    @property
    def cos_b_sq(self):
        return getattr(self, "__cos_b_sq")
    @property
    def cos_g_sq(self):
        return getattr(self, "__cos_g_sq")
    @property
    def sin_a_sq(self):
        return getattr(self, "__sin_a_sq")
    @property
    def sin_b_sq(self):
        return getattr(self, "__sin_b_sq")
    @property
    def sin_g_sq(self):
        return getattr(self, "__sin_g_sq")

    @property
    def cos_ia(self):
        return getattr(self, "__cos_ia")
    @property
    def cos_ib(self):
        return getattr(self, "__cos_ib")
    @property
    def cos_ig(self):
        return getattr(self, "__cos_ig")
    @property
    def sin_ia(self):
        return getattr(self, "__sin_ia")
    @property
    def sin_ib(self):
        return getattr(self, "__sin_ib")
    @property
    def sin_ig(self):
        return getattr(self, "__sin_ig")
    @property
    def cos_ia_sq(self):
        return getattr(self, "__cos_ia_sq")
    @property
    def cos_ib_sq(self):
        return getattr(self, "__cos_ib_sq")
    @property
    def cos_ig_sq(self):
        return getattr(self, "__cos_ig_sq")
    @property
    def sin_ia_sq(self):
        return getattr(self, "__sin_ia_sq")
    @property
    def sin_ib_sq(self):
        return getattr(self, "__sin_ib_sq")
    @property
    def sin_ig_sq(self):
        return getattr(self, "__sin_ig_sq")

    @property
    def is_variable(self):
        """
Output: True if there is any refined parameter
        """
        res = any([self.length_a.refinement,
                   self.length_b.refinement,
                   self.length_c.refinement,
                   self.angle_alpha.refinement,
                   self.angle_beta.refinement,
                   self.angle_gamma.refinement])
        return res
        
    def get_variables(self):
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.length_a.refinement: l_variable.append(self.length_a)
        if self.length_b.refinement: l_variable.append(self.length_b)
        if self.length_c.refinement: l_variable.append(self.length_c)
        if self.angle_alpha.refinement: l_variable.append(self.angle_alpha)
        if self.angle_beta.refinement: l_variable.append(self.angle_beta)
        if self.angle_gamma.refinement: l_variable.append(self.angle_gamma)
        return l_variable
    
    def estimate_type_cell(self):
        rel_tol, abs_tol = 0.0001, 0.0001
        flag_a_b = math.isclose(float(self.length_a), float(self.length_b), rel_tol=rel_tol, abs_tol=abs_tol)
        flag_a_c = math.isclose(float(self.length_a), float(self.length_c), rel_tol=rel_tol, abs_tol=abs_tol)
        flag_b_c = math.isclose(float(self.length_b), float(self.length_c), rel_tol=rel_tol, abs_tol=abs_tol)
        flag_alpha_beta = math.isclose(float(self.angle_alpha), float(self.angle_beta), rel_tol=rel_tol, abs_tol=abs_tol)
        flag_alpha_gamma = math.isclose(float(self.angle_alpha), float(self.angle_gamma), rel_tol=rel_tol, abs_tol=abs_tol)
        flag_alpha_90 = math.isclose(float(self.angle_alpha), 90., rel_tol=rel_tol, abs_tol=abs_tol)
        flag_beta_90 = math.isclose(float(self.angle_beta), 90., rel_tol=rel_tol, abs_tol=abs_tol)
        flag_gamma_90 = math.isclose(float(self.angle_gamma), 90., rel_tol=rel_tol, abs_tol=abs_tol)
        flag_gamma_120 = math.isclose(float(self.angle_gamma), 120., rel_tol=rel_tol, abs_tol=abs_tol)

        if all([flag_a_b, flag_a_c, flag_alpha_90, flag_beta_90, flag_gamma_90]):
            type_cell = "cP"
        elif all([flag_a_b, flag_a_c, flag_alpha_beta, flag_alpha_gamma, not(flag_alpha_90)]):
            type_cell = "hR"
        elif all([flag_a_b, flag_alpha_90, flag_beta_90, flag_gamma_120]):
            type_cell = "hR"
        elif all([flag_a_b, not(flag_a_c), flag_alpha_90, flag_beta_90, flag_gamma_90]):
            type_cell = "hP"
        elif all([not(flag_a_b), not(flag_a_c), not(flag_b_c), flag_alpha_90, flag_beta_90, flag_gamma_90]):
            type_cell = "oP"
        elif all([not(flag_a_b), not(flag_a_c), not(flag_b_c), flag_alpha_90, not(flag_beta_90), flag_gamma_90]):
            type_cell = "mP"
        else:
            type_cell = "aP"
        return type_cell

    @property
    def form_object(self):
        if self.is_defined_attribute("type_cell"):
            CONSTANTS_AND_FUNCTIONS.apply_constraint_on_cell_by_type_cell(self, self.type_cell, self.it_coordinate_system_code)
        flag = True
        rad=numpy.pi/180.

        c_a = numpy.cos(float(self.angle_alpha)*rad)
        c_b = numpy.cos(float(self.angle_beta)*rad)
        c_g = numpy.cos(float(self.angle_gamma)*rad)
        s_a = numpy.sin(float(self.angle_alpha)*rad)
        s_b = numpy.sin(float(self.angle_beta)*rad)
        s_g = numpy.sin(float(self.angle_gamma)*rad)
        setattr(self, "__cos_a", c_a)
        setattr(self, "__cos_b", c_b)
        setattr(self, "__cos_g", c_g)
        setattr(self, "__sin_a", s_a)
        setattr(self, "__sin_b", s_b)
        setattr(self, "__sin_g", s_g)

        c_a_sq, c_b_sq, c_g_sq = c_a**2, c_b**2, c_g**2
        s_a_sq, s_b_sq, s_g_sq = (1.-c_a_sq), (1.-c_b_sq), (1.-c_g_sq)
        setattr(self, "__cos_a_sq", c_a_sq)
        setattr(self, "__cos_b_sq", c_b_sq)
        setattr(self, "__cos_g_sq", c_g_sq)
        setattr(self, "__sin_a_sq", s_a_sq)
        setattr(self, "__sin_b_sq", s_b_sq)
        setattr(self, "__sin_g_sq", s_g_sq)

        a, b, c = float(self.length_a), float(self.length_b), float(self.length_c)

        vol = a*b*c*(1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)**0.5
        setattr(self, "__volume", vol)
    
        irad = 180./numpy.pi
        ialpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))*irad
        ibeta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))*irad
        igamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))*irad
        ia, ib, ic = b*c*s_a/vol, c*a*s_b/vol, a*b*s_g/vol

        setattr(self, "__reciprocal_length_a", ia)
        setattr(self, "__reciprocal_length_b", ib)
        setattr(self, "__reciprocal_length_c", ic)
        setattr(self, "__reciprocal_angle_alpha", ialpha)
        setattr(self, "__reciprocal_angle_beta", ibeta)
        setattr(self, "__reciprocal_angle_gamma", igamma)

        
        c_ia, c_ib, c_ig = numpy.cos(ialpha*rad), numpy.cos(ibeta*rad), numpy.cos(igamma*rad)
        s_ia, s_ib, s_ig = numpy.sin(ialpha*rad), numpy.sin(ibeta*rad), numpy.sin(igamma*rad)
        setattr(self, "__cos_ia", c_ia)
        setattr(self, "__cos_ib", c_ib)
        setattr(self, "__cos_ig", c_ig)
        setattr(self, "__sin_ia", s_ia)
        setattr(self, "__sin_ib", s_ib)
        setattr(self, "__sin_ig", s_ig)

        c_ia_sq, c_ib_sq, c_ig_sq = c_ia**2, c_ib**2, c_ig**2
        s_ia_sq, s_ib_sq, s_ig_sq = (1.-c_ia_sq), (1.-c_ib_sq), (1.-c_ig_sq)
        setattr(self, "__cos_ia_sq", c_ia_sq)
        setattr(self, "__cos_ib_sq", c_ib_sq)
        setattr(self, "__cos_ig_sq", c_ig_sq)
        setattr(self, "__sin_ia_sq", s_ia_sq)
        setattr(self, "__sin_ib_sq", s_ib_sq)
        setattr(self, "__sin_ig_sq", s_ig_sq)

        #B matrix
        m_b = numpy.array([[ia,  ib*c_ig,      ic*c_ib],
                           [0.,  ib*s_ig, -ic*s_ib*c_a],
                           [0.,       0.,         1./c]], dtype = float)
        setattr(self, "__m_b", m_b)

        #inversed B matrix
        m_ib = numpy.linalg.inv(m_b)
        setattr(self, "__m_ib", m_ib)

        m_ib_norm = numpy.copy(m_ib)
        m_ib_norm[0, :] /= a
        m_ib_norm[1, :] /= b
        m_ib_norm[2, :] /= c
        setattr(self, "__m_ib_norm", m_ib_norm)

        #M matrix (not sure)
        m_m = numpy.array([[ a,   b*c_g,       c*c_b],
                           [0.,   b*s_g, -c*s_b*c_ia],
                           [0.,      0.,       1./ic]], dtype = float)
        setattr(self, "__m_m", m_m)
        return flag


    def apply_constraint(self, type_cell:str, it_coordinate_system_code:str)->bool:
        setattr(self, "__type_cell", type_cell)
        setattr(self, "__it_coordinate_system_code", it_coordinate_system_code)
        CONSTANTS_AND_FUNCTIONS.apply_constraint_on_cell_by_type_cell(self, self.type_cell, self.it_coordinate_system_code)
        flag = False
        if self.is_defined:
            flag = self.form_object
        return flag



    def calc_sthovl(self, h=None, k=None, l=None, hkl=None, l_hkl=None):
        """
Calculate sin(theta)/lambda for list of hkl reflections.

Keyword arguments:

 - h, k, l --- Miller indices (example: h=1, k=0, l=2)

or

 - hkl --- a turple of Miller indices (example: hkl=(1, 0, 2))

or

 - l_hkl --- a list of turples of Miller indices (example: l_hkl=[(1, 0, 2), (1, 1, 3)])
 - f_print --- a flag to print output information in terminal
        """
        cond_1 = all([h is not None, k is not None, l is not None])
        cond_2 = hkl is not None
        cond_3 = l_hkl is not None
        if cond_1:
            np_h = numpy.array(h, dtype=int) # type: numpy.array or float
            np_k = numpy.array(k, dtype=int) # type: numpy.array or float
            np_l = numpy.array(l, dtype=int) # type: numpy.array or float
        elif cond_2:
            np_h = hkl[0] # type: float
            np_k = hkl[1] # type: float
            np_l = hkl[2] # type: float
        elif cond_3:
            np_h = numpy.array([hh[0] for hh in l_hkl], dtype=float) # type: numpy.array
            np_k = numpy.array([hh[1] for hh in l_hkl], dtype=float) # type: numpy.array
            np_l = numpy.array([hh[2] for hh in l_hkl], dtype=float) # type: numpy.array
        else: 
            self._show_message("Did not found correct input. Expected h, k, l or hkl or l_hkl")
            return
        a = float(self.length_a)
        b = float(self.length_b)
        c = float(self.length_c)
        c_a, c_b, c_g = self.cos_a, self.cos_b, self.cos_g
        c_a_sq, c_b_sq, c_g_sq = self.cos_a_sq, self.cos_b_sq, self.cos_g_sq
        s_a_sq, s_b_sq, s_g_sq = self.sin_a_sq, self.sin_b_sq, self.sin_g_sq

        A=( 1. - c_a_sq - c_b_sq - c_g_sq + 2.*c_a*c_b*c_g)
        B1 = (s_a_sq*(np_h*1./a)**2+s_b_sq*(np_k*1./b)**2+s_g_sq*(np_l*1./c)**2)
        B2 = 2.*(np_k*np_l*c_a)/(b*c)+2.*(np_h*np_l*c_b)/(a*c)+2.*(np_h*np_k*c_g)/(a*b)
        #it should be checked, I am not sure
        B = B1-B2
        inv_d = (B*1./A)**0.5
        res = 0.5*inv_d
        return res

    def calc_k_loc(self, h, k, l):
        """
Calculate unity scattering vector in Cartesian coordinate system (x||a*, z||c)

Keyword arguments:

    h, k, l: Miller indices

Output arguments:

    k_x, k_y, k_z: 1D numpy array of x, y, z components of unity scattering vector 
        """
        np_h = numpy.array(h, dtype=int)
        np_k = numpy.array(k, dtype=int)
        np_l = numpy.array(l, dtype=int)
        m_b = self.m_b
        k_x = m_b[0, 0]*np_h + m_b[0, 1]*np_k +m_b[0, 2]*np_l
        k_y = m_b[1, 0]*np_h + m_b[1, 1]*np_k +m_b[1, 2]*np_l
        k_z = m_b[2, 0]*np_h + m_b[2, 1]*np_k +m_b[2, 2]*np_l
        
        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        if not((type(h) is float)|(type(h) is int)|(type(h) is numpy.float64)):
            k_norm[k_norm == 0.] = 1.
        elif k_norm == 0.:
            k_norm = 1.
        
        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm
        
        return k_x, k_y, k_z
        
    def calc_m_t(self, h, k, l):
        """
Determine rotation matrix to have new z axis along kloc
Rotation matrix is defined by Euler angles
        """
        m_b = self.m_b
        k_x = m_b[0, 0]*h + m_b[0, 1]*k +m_b[0, 2]*l
        k_y = m_b[1, 0]*h + m_b[1, 1]*k +m_b[1, 2]*l
        k_z = m_b[2, 0]*h + m_b[2, 1]*k +m_b[2, 2]*l
        
        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        k_norm[k_norm == 0.] = 1.
        
        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm
        

        al = numpy.zeros(k_x.shape, dtype=float)
        
        be = numpy.arccos(k_z)
        sb = numpy.sin(be)
        flag = (sb != 0.)
        
        sa1 = k_x[flag]*1./sb[flag]
        ca2 = -1*k_y[flag]*1./sb[flag]
        sa1[sa1>1] = 1.
        sa1[sa1<-1] = -1.
            
        ca2[ca2>1] = 1.
        ca2[ca2<-1] = -1.

        al1 = numpy.arcsin(sa1)
        al2 = numpy.arccos(ca2)
        
        al_sh = numpy.copy(al1)
        al_sh[sa1 > 0.] = al2[sa1 > 0.]
        al_sh[sa1 <= 0.] = 2.*numpy.pi-al2[sa1 <= 0.]
        al_sh[numpy.abs(al2-al1)<0.00001] = al1[numpy.abs(al2-al1)<0.00001]

        al[flag] = al_sh
            
        ga=0.
        ca, cb, cg = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
        sa, sb, sg = numpy.sin(al), numpy.sin(be), numpy.sin(ga)
        t_11, t_12, t_13 = ca*cg-sa*cb*sg, -ca*sg-sa*cb*cg,  sa*sb
        t_21, t_22, t_23 = sa*cg+ca*cb*sg, -sa*sg+ca*cb*cg, -ca*sb
        t_31, t_32, t_33 =          sb*sg,           sb*cg,     cb
        
        flag = (((sa*sb-k_x)**2+(-ca*sb-k_y)**2+(cb-k_z)**2)>0.0001)
        if any(flag):
            warnings.warn("Mistake with k_loc\nProgram is stopped", UserWarning, stacklevel=2)
            quit()
        return t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 


    def calc_hkl(self, space_group, sthovl_min, sthovl_max):
        """
A list of reflections hkl for cell in the range sthovl_min, sthovl_max
taking into account the space group
        """
        if not(self.is_defined):
            warnings.warn("Object 'Cell' is not fully defined for calculations", UserWarning, stacklevel=2)
            return None
        lhkl,lmult=[],[]
        l_hklres=[]

        hmax = int(2.*self.length_a*sthovl_max)
        kmax = int(2.*self.length_b*sthovl_max)
        lmax = int(2.*self.length_c*sthovl_max)
        hmin, kmin, lmin = -1*hmax, -1*kmax, -1*lmax

        hmin=0
        shift = space_group.shift
        r_s_g_s = space_group.reduced_space_group_symop
        orig_x, orig_y, orig_z = zip(*shift)
        np_orig_x, np_orig_y, np_orig_z = numpy.array(orig_x, dtype=float), numpy.array(orig_y, dtype=float), numpy.array(orig_z, dtype=float), 
        
        r_11, r_12, r_13 = numpy.array(r_s_g_s.r_11, dtype=float), numpy.array(r_s_g_s.r_12, dtype=float), numpy.array(r_s_g_s.r_13, dtype=float)
        r_21, r_22, r_23 = numpy.array(r_s_g_s.r_21, dtype=float), numpy.array(r_s_g_s.r_22, dtype=float), numpy.array(r_s_g_s.r_23, dtype=float)
        r_31, r_32, r_33 = numpy.array(r_s_g_s.r_31, dtype=float), numpy.array(r_s_g_s.r_32, dtype=float), numpy.array(r_s_g_s.r_33, dtype=float)
        for h in range(hmin,hmax+1,1):
            for k in range(kmin,kmax+1,1):
                for l in range(lmin,lmax+1,1):
                    flag=(abs(sum(numpy.exp(2.*numpy.pi*1j*(np_orig_x*h+np_orig_y*k+np_orig_z*l))))>0.00001)
                    #flag=True
                    if (flag):
                        lhkls = [(_1, _2, _3) for _1, _2, _3 in zip(h*r_11+k*r_21+l*r_31, h*r_12+k*r_22+l*r_32, h*r_13+k*r_23+l*r_33)]
                        lhkls.extend([(-hkl[0],-hkl[1],-hkl[2]) for hkl in lhkls])
                        lhkls.sort(key=lambda x:10000*x[0]+100*x[1]+x[2])
                        if (not(lhkls[-1] in lhkl)):
                            lhkl.append(lhkls[-1])
                            lmult.append(len(set(lhkls)))
                            
        l_hklsthovl=[(hkl, self.calc_sthovl(hkl[0], hkl[1], hkl[2]), mult) for hkl, mult in zip(lhkl, lmult)]
        l_hklsthovl.sort(key=lambda x: x[1])
        l_hklres = [hklsthovl[0] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]
        l_multres = [hklsthovl[2] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]

        h = numpy.array([hh[0] for hh in l_hklres], dtype=int)
        k = numpy.array([hh[1] for hh in l_hklres], dtype=int)
        l = numpy.array([hh[2] for hh in l_hklres], dtype=int)
        mult = numpy.array(l_multres, dtype=int)
        return h, k, l, mult


    def calc_hkl_in_range(self, sthovl_min, sthovl_max):
        """
Give a list of reflections hkl for cell in the range sthovl_min, sthovl_max
        """
        if not(self.is_defined):
            print("Object 'Cell' is not fully defined for calculations")
            return None
        h_max = int(2.*self.length_a*sthovl_max)
        k_max = int(2.*self.length_b*sthovl_max)
        l_max = int(2.*self.length_c*sthovl_max)
        h_min, k_min, l_min = -1*h_max, -1*k_max, -1*l_max

        np_h = numpy.array(range(h_min, h_max+1, 1), dtype=int)
        np_k = numpy.array(range(k_min, k_max+1, 1), dtype=int)
        np_l = numpy.array(range(l_min, l_max+1, 1), dtype=int)
        h_3d, k_3d, l_3d = numpy.meshgrid(np_h, np_k, np_l, indexing="ij")
        
        sthovl_3d = self.calc_sthovl(h_3d, k_3d, l_3d)
        flag_1 = sthovl_3d >= sthovl_min
        flag_2 = sthovl_3d <= sthovl_max
        flag_12 = numpy.logical_and(flag_1, flag_2)
        
        h = h_3d[flag_12]
        k = k_3d[flag_12]
        l = l_3d[flag_12]
        mult = numpy.ones(h.size, dtype=int)
        sthovl = sthovl_3d[flag_12]
        arg_sort = numpy.argsort(sthovl)
        return h[arg_sort], k[arg_sort], l[arg_sort], mult[arg_sort] 
        
    def calc_position_by_coordinate(self, x ,y, z):
        """
Calculates position for coordinate :math:`(x,y,z)`
in Carthezian coordinate system where :math:`X` along axis :math:`a`,
:math:`Z` alond :math:`c^{*}` and :math:`Y = Z \\times X`


        """
        m_m = self.m_m
        p_x = m_m[0, 0]*x + m_m[0, 1]*y + m_m[0, 2]*z
        p_y = m_m[1, 0]*x + m_m[1, 1]*y + m_m[1, 2]*z
        p_z = m_m[2, 0]*x + m_m[2, 1]*y + m_m[2, 2]*z
        return p_x, p_y, p_z

    def calc_length_sq(self, x, y, z):
        """
According to IT_C Section 1.1.2 Lattice vectors, point rows and net planes

.. math::
    t^{2} = x^{2} a^{2} + y^{2} b^{2} + z^{2} c^{2} 
            + 2 x y a b \\cos \\gamma  + 2 x z a c \\cos \\beta  + 2 y z b c \\cos \\alpha 
        """
        a, b, c = float(self.length_a), float(self.length_b), float(self.length_c)
        c_a, c_b, c_g = self.cos_a, self.cos_b, self.cos_g
        t_sq = (x**2 * a**2 + y**2 * b**2 + z**2 * c**2 + 
                2.*x*y*a*b*c_g + 2.*x*z*a*c*c_b + 2.*y*z*b*c*c_a)
        return t_sq

    def calc_reciprocal_length_sq(self, h, k, l):
        """
According to IT_C Section 1.1.2 Lattice vectors, point rows and net planes

.. math::
    r^{*^{2}} = h^{2} a^{*^{2}} + k^{*^{2}} b^{*^{2}} + l^{2} c^{*^{2}} 
            + 2 h k a^{*} b^{*} \\cos \\gamma  + 2 h l a^{*} c^{*} \\cos \\beta  + 2 k l b^{*} c^{*} \\cos \\alpha 
        """
        a, b, c = float(self.reciprocal_length_a), float(self.reciprocal_length_b), float(self.reciprocal_length_c)
        c_a, c_b, c_g = self.cos_ia, self.cos_ib, self.cos_ig
        x, y, z = h, k, l
        t_sq = (x**2 * a**2 + y**2 * b**2 + z**2 * c**2 + 
                2.*x*y*a*b*c_g + 2.*x*z*a*c*c_b + 2.*y*z*b*c*c_a)
        return t_sq
