__author__ = 'ikibalin'
__version__ = "2019_12_06"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr


class Refln(ItemConstr):
    """
Refln
==============
Data items in the REFLN category record details about the
reflections used to determine the ATOM_SITE data items.
The REFLN data items refer to individual reflections and must
be included in looped lists.

The REFLNS data items specify the parameters that apply to all
reflections. The REFLNS data items are not looped.

Description in cif file:
-------------------------
::

  _refln_index_h         2
  _refln_index_k         0
  _refln_index_l         0
  _refln_d_spacing       0.13
  _refln_A_calc          2.1
  _refln_B_calc          0.
  _refln_chi_11_A_calc   0.
  _refln_chi_12_A_calc   0.
  _refln_chi_13_A_calc   0.
  _refln_chi_21_A_calc   0.
  _refln_chi_22_A_calc   0.
  _refln_chi_23_A_calc   0.
  _refln_chi_31_A_calc   0.
  _refln_chi_32_A_calc   0.
  _refln_chi_33_A_calc   0.
  _refln_chi_11_B_calc   0.
  _refln_chi_12_B_calc   0.
  _refln_chi_13_B_calc   0.
  _refln_chi_21_B_calc   0.
  _refln_chi_22_B_calc   0.
  _refln_chi_23_B_calc   0.
  _refln_chi_31_B_calc   0.
  _refln_chi_32_B_calc   0.
  _refln_chi_33_B_calc   0.
 
Reference:
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefln.html>`_

FIXME:
the attribute 'sint/lambda'is replaced by 'sintlambda'
    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("a_calc", "a_meas", "b_calc", "b_meas", 
                          "class_code", "crystal_id", "d_spacing",
                          "f_calc", "f_meas", "f_sigma", "f_squared_calc", "f_squared_meas", "f_squared_sigma",
                          "include_status", "intensity_calc", "intensity_meas", "intensity_sigma", "mean_path_length_tbar",
                          "phase_calc", "phase_meas", "refinement_status", "scale_group_code",
                          "sintlambda", "symmetry_epsilon", "symmetry_multiplicity", "wavelength", "wavelength_id")
    INTERNAL_ATTRIBUTE = ()
    ACCESIBLE_INCLUDE_STATUS = ("o", "<", "-", "x", "h", "r")
    ACCESIBLE_REFINEMENT_STATUS = ("incl", "excl", "extn")
    PREFIX = "refln"
    def __init__(self, index_h=None, index_k=None, index_l=None, 
    a_calc=None, a_meas=None, b_calc=None, b_meas=None, class_code=None, crystal_id=None, d_spacing=None,
    f_calc=None, f_meas=None, f_sigma=None, f_squared_calc=None, f_squared_meas=None, f_squared_sigma=None,
    include_status=None, intensity_calc=None, intensity_meas=None, intensity_sigma=None, mean_path_length_tbar=None,
    phase_calc=None, phase_meas=None, refinement_status=None, scale_group_code=None,
    sintlambda=None, symmetry_epsilon=None, symmetry_multiplicity=None, wavelength=None, wavelength_id=None):
        super(Refln, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE,
                                    optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                    internal_attribute=self.INTERNAL_ATTRIBUTE,
                                    prefix=self.PREFIX)
        self.index_h = index_h
        self.index_k = index_k
        self.index_l = index_l
        self.a_calc = a_calc
        self.a_meas = a_meas
        self.b_calc = b_calc
        self.b_meas = b_meas
        self.class_code = class_code
        self.crystal_id = crystal_id
        self.d_spacing = d_spacing
        self.f_calc = f_calc
        self.f_meas = f_meas
        self.f_sigma = f_sigma
        self.f_squared_calc = f_squared_calc
        self.f_squared_meas = f_squared_meas
        self.f_squared_sigma = f_squared_sigma
        self.include_status = include_status
        self.intensity_calc = intensity_calc
        self.intensity_meas = intensity_meas
        self.intensity_sigma = intensity_sigma
        self.mean_path_length_tbar = mean_path_length_tbar
        self.phase_calc = phase_calc
        self.phase_meas = phase_meas
        self.refinement_status = refinement_status
        self.scale_group_code = scale_group_code
        self.sintlambda = sintlambda
        self.symmetry_epsilon = symmetry_epsilon
        self.symmetry_multiplicity = symmetry_multiplicity
        self.wavelength = wavelength
        self.wavelength_id = wavelength_id

        if self.is_defined:
            self.form_object

    @property
    def index_h(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_h")
    @index_h.setter
    def index_h(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_h", x_in)
    @property
    def index_k(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_k")
    @index_k.setter
    def index_k(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_k", x_in)
    @property
    def index_l(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_l")
    @index_l.setter
    def index_l(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_l", x_in)
    @property
    def a_calc(self):
        return getattr(self, "__a_calc")
    @a_calc.setter
    def a_calc(self, x):
        """
The calculated structure-factor component A
(in electrons for X-ray diffraction).

   A =|F|cos(phase)

Appears in list containing _refln_index_

Type: numb

Category: refln        
        """
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__a_calc", x_in)
    @property
    def a_meas(self):
        """
The measured structure-factor component A
(in electrons for X-ray diffraction).

   A =|F|cos(phase)

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__a_meas")
    @a_meas.setter
    def a_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__a_meas", x_in)
    @property
    def b_calc(self):
        return getattr(self, "__b_calc")
    @b_calc.setter
    def b_calc(self, x):
        """
The calculated structure-factor component B
(in electrons for X-ray diffraction).

   B =|F|sin(phase)

Appears in list containing _refln_index_

Type: numb        
        """
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__b_calc", x_in)
    @property
    def b_meas(self):
        """
The measured structure-factor component B
(in electrons for X-ray diffraction).

   B =|F|sin(phase)

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__b_meas")
    @b_meas.setter
    def b_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__b_meas", x_in)
    @property
    def class_code(self):
        """
The code identifying the class to which this reflection has been
assigned. This code must match a value of _reflns_class_code.
Reflections may be grouped into classes for a variety of
purposes. For example, for modulated structures each reflection
class may be defined by the number m=sum|m~i~|, where the m~i~
are the integer coefficients that, in addition to h,k,l, index
the corresponding diffraction vector in the basis defined
for the reciprocal lattice.

Appears in list containing _refln_index_

Must match data name_reflns_class_code

Type: char        
        """
        return getattr(self, "__class_code")
    @class_code.setter
    def class_code(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__class_code", x_in)
    @property
    def crystal_id(self):
        """
Code identifying each crystal if multiple crystals are used. Is
used to link with _exptl_crystal_id in the _exptl_crystal_ list.

Appears in list containing _refln_index_

Must match data name_exptl_crystal_id

Type: char        
        """
        return getattr(self, "__crystal_id")
    @crystal_id.setter
    def crystal_id(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__crystal_id", x_in)
    @property
    def d_spacing(self):
        """
The d spacing in angstroms for this reflection. This is related
to the (sin theta)/lambda value by the expression
_refln_d_spacing = 2/(_refln_sint/lambda)

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__d_spacing")
    @d_spacing.setter
    def d_spacing(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__d_spacing", x_in)
    @property
    def f_calc(self):
        """
The calculated  of the structure factors (in electrons for
X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__f_calc")
    @f_calc.setter
    def f_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_calc", x_in)
    @property
    def f_meas(self):
        """
The measured of the structure factors (in electrons for
X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__f_meas")
    @f_meas.setter
    def f_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_meas", x_in)
    @property
    def f_sigma(self):
        """
The standard uncertainty (derived from
measurement) of the structure factors (in electrons for
X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """

        return getattr(self, "__f_sigma")
    @f_sigma.setter
    def f_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_sigma", x_in)
    @property
    def f_squared_calc(self):
        """
Calculated of the squared structure factors (in electrons
squared for X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__f_squared_calc")
    @f_squared_calc.setter
    def f_squared_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_squared_calc", x_in)
    @property
    def f_squared_meas(self):
        """
Measured of the squared structure factors (in electrons
squared for X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__f_squared_meas")
    @f_squared_meas.setter
    def f_squared_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_squared_meas", x_in)
    @property
    def f_squared_sigma(self):
        """
Estimated standard uncertainty (derived
from measurement) of the squared structure factors (in electrons
squared for X-ray diffraction).

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__f_squared_sigma")
    @f_squared_sigma.setter
    def f_squared_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_squared_sigma", x_in)
    @property
    def include_status(self):
        """
Classification of a reflection indicating its status with
respect to inclusion in the refinement and the calculation
of R factors.

Appears in list containing _refln_index_

Related item: _refln_observed_status (alternate)

The data value must be one of the following:


o	(lower-case letter o for 'observed') satisfies _refine_ls_d_res_high satisfies _refine_ls_d_res_low exceeds _reflns_threshold_expression
<	satisfies _refine_ls_d_res_high satisfies _refine_ls_d_res_low does not exceed _reflns_threshold_expression
-	systematically absent reflection
x	unreliable measurement -- not used
h	does not satisfy _refine_ls_d_res_high
l	does not satisfy _refine_ls_d_res_low

Enumeration default: o
Type: char        
        """
        return getattr(self, "__include_status")
    @include_status.setter
    def include_status(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_INCLUDE_STATUS):
                warnings.warn(f"include_status '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__include_status", x_in)
    @property
    def intensity_calc(self):
        """
The calculated of the intensity, all in the same arbitrary units
as _refln_intensity_meas.

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__intensity_calc")
    @intensity_calc.setter
    def intensity_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_calc", x_in)
    @property
    def intensity_meas(self):
        """
The measured of the intensity, all in the same arbitrary units
as _refln_intensity_meas.

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__intensity_meas")
    @intensity_meas.setter
    def intensity_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_meas", x_in)
    @property
    def intensity_sigma(self):
        """
The standard uncertainty (derived from
measurement) of the intensity, all in the same arbitrary units
as _refln_intensity_meas.

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__intensity_sigma")
    @intensity_sigma.setter
    def intensity_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_sigma", x_in)
    @property
    def mean_path_length_tbar(self):
        """
Mean path length in millimetres through the crystal for this
reflection.

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__mean_path_length_tbar")
    @mean_path_length_tbar.setter
    def mean_path_length_tbar(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__mean_path_length_tbar", x_in)
    @property
    def phase_calc(self):
        """
The calculated structure-factor phase in degrees.

Appears in list containing _refln_index_

Type: numb        
        """
        return getattr(self, "__phase_calc")
    @phase_calc.setter
    def phase_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phase_calc", x_in)
    @property
    def phase_meas(self):
        """
The measured structure-factor phase in degrees.

Appears in list containing _refln_index_

Type: numb

Category: refln

        """
        return getattr(self, "__phase_meas")
    @phase_meas.setter
    def phase_meas(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phase_meas", x_in)
    @property
    def refinement_status(self):
        """
Status of a reflection in the structure-refinement process.

Appears in list containing _refln_index_ 
The data value must be one of the following:

incl	included in ls process

excl	excluded from ls process

extn	excluded due to extinction

Enumeration default: incl        
        """
        return getattr(self, "__refinement_status")
    @refinement_status.setter
    def refinement_status(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)

            if not(x_in in self.ACCESIBLE_REFINEMENT_STATUS):
                warnings.warn(f"refinement_status '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None                    
            setattr(self, "__refinement_status", x_in)
    @property
    def scale_group_code(self):
        """
Code identifying the structure-factor scale. This code must
   correspond to one of the _reflns_scale_group_code values.

Examples:
1, 2, 3, s1
A, B, c1, c2, c3

Appears in list containing _refln_index_

Must match data name_reflns_scale_group_code

Type: char        
        """
        return getattr(self, "__scale_group_code")
    @scale_group_code.setter
    def scale_group_code(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__scale_group_code", x_in)
    @property
    def sintlambda(self):
        """
The (sin theta)/lambda value in reciprocal angstroms for this
   reflection.

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__sintlambda")
    @sintlambda.setter
    def sintlambda(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__sintlambda", x_in)
    @property
    def symmetry_epsilon(self):
        """
The symmetry reinforcement factor corresponding to the number of
times the reflection indices are generated identically from the
space-group symmetry operations.

Appears in list containing _refln_index_ 
The permitted range is 1 -> 48

Type: numb        
        """
        return getattr(self, "__symmetry_epsilon")
    @symmetry_epsilon.setter
    def symmetry_epsilon(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__symmetry_epsilon", x_in)
    @property
    def symmetry_multiplicity(self):
        """
The number of reflections symmetry-equivalent under the Laue
symmetry to the present reflection. In the Laue symmetry, Friedel
opposites (h k l and -h -k -l) are equivalent. Tables of
symmetry-equivalent reflections are available in International
Tables for Crystallography Volume A (2002), Chapter 10.1.

Appears in list containing _refln_index_ 
The permitted range is 1 -> 48

Type: numb        
        """
        return getattr(self, "__symmetry_multiplicity")
    @symmetry_multiplicity.setter
    def symmetry_multiplicity(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__symmetry_multiplicity", x_in)
    @property
    def wavelength(self):
        """
The mean wavelength in angstroms of the radiation used to measure
this reflection. This is an important parameter for data
collected using energy-dispersive detectors or the Laue method.

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__wavelength")
    @wavelength.setter
    def wavelength(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__wavelength", x_in)
    @property
    def wavelength_id(self):
        """
Code identifying the wavelength in the _diffrn_radiation_ list.
See _diffrn_radiation_wavelength_id.

Appears in list containing _refln_index_

Must match data name_diffrn_radiation_wavelength_id

Type: char        
        """
        return getattr(self, "__wavelength_id")
    @wavelength_id.setter
    def wavelength_id(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__wavelength_id", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append("Refln:")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []


class ReflnL(LoopConstr):
    """
ReflnL
==============
Data items in the REFLN category record details about the
reflections used to determine the ATOM_SITE data items.
The REFLN data items refer to individual reflections and must
be included in looped lists.

The REFLNS data items specify the parameters that apply to all
reflections. The REFLNS data items are not looped.

Description in cif file:
-------------------------
::

 loop_
  _refln_index_h
  _refln_index_k
  _refln_index_l
  _refln_d_spacing
  _refln_A_calc
  _refln_B_calc
  0 0 2 2.315 3.25  1.232
  2 2 0 4.213 5.00 -4.05
 
Reference:
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefln.html>`_
    """
    CATEGORY_KEY = ("index_h", "index_k", "index_l")
    ITEM_CLASS = Refln
    def __init__(self, item=[], loop_name=""):
        super(ReflnL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("ReflnL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
