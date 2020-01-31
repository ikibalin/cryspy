__author__ = 'ikibalin'
__version__ = "2020_01_21"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr


class RefineLs(ItemConstr):
    """
 Data items in the REFINE_LS category record details about the
 structure-refinement parameters.

Description in cif file::

  _refine_ls_goodness_of_fit_all        2.74
  _refine_ls_weighting_scheme           sigma
  _refine_ls_number_reflns           1408
  _refine_ls_number_parameters       272
  _refine_ls_number_restraints       0
  _refine_ls_number_constraints      0

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefine.html>`_

FIXME:
expression for goodness_of_fit_all is changed.
    """
    MANDATORY_ATTRIBUTE = ("number_reflns", "goodness_of_fit_all", "number_parameters")
    OPTIONAL_ATTRIBUTE = ("number_restraints", "number_constraints", "weighting_scheme")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("number_reflns", "goodness_of_fit_all", "number_parameters")
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("number_restraints", "number_constraints", "weighting_scheme")
    INTERNAL_ATTRIBUTE = ()
    ACCESIBLE_WEIGHTING_SCHEME = ("sigma", "unit", "calc")
    PREFIX = "refine_ls"
    def __init__(self, number_reflns=None, goodness_of_fit_all=None, number_parameters=None, 
    number_restraints=None, number_constraints=None, weighting_scheme=None):
        super(RefineLs, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE,
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)
        self.number_reflns = number_reflns
        self.goodness_of_fit_all = goodness_of_fit_all
        self.number_parameters = number_parameters
        self.number_restraints = number_restraints
        self.number_constraints = number_constraints
        self.weighting_scheme = weighting_scheme

        if self.is_defined:
            self.form_object

    @property
    def number_reflns(self):
        """
The number of unique reflections contributing to the
least-squares refinement calculation.

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Irefine_ls_number_reflns.html>`_
        """
        return getattr(self, "__number_reflns")
    @number_reflns.setter
    def number_reflns(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__number_reflns", x_in)

    @property
    def number_parameters(self):
        """
The number of parameters refined in the least-squares process.
If possible, this number should include some contribution from
the restrained parameters. The restrained parameters are
distinct from the constrained parameters (where one or more
parameters are linearly dependent on the refined value of
another). Least-squares restraints often depend on geometry or
energy considerations and this makes their direct contribution
to this number, and to the goodness-of-fit calculation,
difficult to assess.  

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Irefine_ls_number_parameters.html>`_
        """
        return getattr(self, "__number_parameters")
    @number_parameters.setter
    def number_parameters(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__number_parameters", x_in)

    @property
    def number_restraints(self):
        """
The number of restrained parameters. These are parameters which
are not directly dependent on another refined parameter.
Restrained parameters often involve geometry or energy
dependencies.


`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Irefine_ls_number_restraints.html>`_
        """
        return getattr(self, "__number_restraints")
    @number_restraints.setter
    def number_restraints(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__number_restraints", x_in)

    @property
    def number_constraints(self):
        """
The number of constrained (non-refined or dependent) parameters
in the least-squares process. These may be due to symmetry or any
other constraint process (e.g. rigid-body refinement). 

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Irefine_ls_number_constraints.html>`_
        """
        return getattr(self, "__number_constraints")
    @number_constraints.setter
    def number_constraints(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__number_constraints", x_in)


    @property
    def goodness_of_fit_all(self):
        """
The least-squares goodness-of-fit parameter S for all
reflections after the final cycle of refinement.
Ideally, account should be taken of parameters restrained
in the least-squares refinement. 

.. math::

    S = \\sqrt{\\frac{\\sum{w \\cdot (Y_{obs}-Y_{calc})^{2}}}{N_{ref}-N_{param}}}

Y(obs) 
    the observed coefficients (see _refine_ls_structure_factor_coef)

Y(calc) 
    the calculated coefficients (see _refine_ls_structure_factor_coef)

w  
    the least-squares reflection weight [1/(u^2^)]

u   
    the standard uncertainty

Nref
    the number of reflections used in the refinement

Nparam
    the number of refined parameters
    and the sum is taken over the specified reflections

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Irefine_ls_goodness_of_fit_all.html>`_
        """
        return getattr(self, "__goodness_of_fit_all")
    @goodness_of_fit_all.setter
    def goodness_of_fit_all(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__goodness_of_fit_all", x_in)
        
    @property
    def weighting_scheme(self):
        """
The weighting scheme applied in the least-squares process. The
standard code may be followed by a description of the weight
(but see _refine_ls_weighting_details for a preferred approach).

The data value must be one of the following:

sigma 
    based on measured s.u.'s 

unit 
    unit or no weights applied 

calc 
    calculated weights applied 

Enumeration default: sigma     
        """
        return getattr(self, "__weighting_scheme")
    @weighting_scheme.setter
    def weighting_scheme(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_WEIGHTING_SCHEME):
                warnings.warn(f"include_status '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__weighting_scheme", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append("RefineLs:")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

