__author__ = 'ikibalin'
__version__ = "2020_04_24"
import os
import math
import numpy
from pycifstar import Global, to_data

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable
import cryspy.corecif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS 
from cryspy.magneticcif.cl_atom_type_scat import AtomTypeScat

import scipy.optimize


def transs(lmax, nn, zeta, sthovl):
    """ 
Calculate integral( r**nn * exp(-zeta r) * j_l(Qr) r**2 dr) / Q**l 
for givel l

Q = 4 pi sint / lambda
    """
    Q = 4. * math.pi * sthovl
    a = [0.*Q for h in range(17)]
    ff = [0.*Q for h in range(lmax+1)]
    d = Q**2+zeta**2
    a[0] = 0.0*Q
    a[1] = 1./d
    n = nn-1
    tz = 2.*zeta
    ts = 4.*math.pi
    for l in range(lmax+1):
        ll = l+1
        if (ll != 1):
            a[ll] = a[ll-1] * ts * l * 1./d
            a[ll-1] = 0.0*Q
        for nx in range(ll, n+1):
            i1 = nx
            i2 = nx+1
            i3 = i2+1
            a[i3-1] = (tz*nx*a[i2-1] - (nx+l)*(nx-ll)*a[i1-1]) * 1./d
        ff[ll-1]=a[i3-1]
    return ff

def calc_GCF(ln1: list, ldzeta1: list, lcoeff1: list, kappa1: float, 
             ln2: list, ldzeta2: list, lcoeff2: list, kappa2: float, 
             sthovl: float, lmax: float):
    GCF = [0. for hh in range(lmax+1)]
    for n1, dzeta1, coeff1 in zip(ln1, ldzeta1, lcoeff1):
        for n2, dzeta2, coeff2 in zip(ln2, ldzeta2, lcoeff2):
            dzetakappa1 = dzeta1*kappa1
            dzetakappa2 = dzeta2*kappa2
            zetaA1 = dzetakappa1/0.529177 #transformation from atomic units to inverse angstrems
            zetaA2 = dzetakappa2/0.529177 #transformation from atomic units to inverse angstrems
            Norm1 = ((2*zetaA1)**(n1+0.5))/math.sqrt(math.factorial(2*n1))
            Norm2 = ((2*zetaA2)**(n2+0.5))/math.sqrt(math.factorial(2*n2))
            nn, zeta = (n1+n2), (zetaA1+zetaA2)
            GCFhelp = [Norm1*Norm2*hh for hh in transs(lmax, nn, zeta, sthovl)]
            GCF = [hh1+coeff1*coeff2*hh2 for hh1, hh2 in zip(GCF, GCFhelp)]
    return GCF

class AtomRhoOrbitalRadialSlater(ItemConstr):
    """
These items are used when the radial dependence of the
electron density, R(k(l), l, r), of the atom specified in
_atom_rho_orbital_atom_label is expressed as a Slater-type
function [Hansen & Coppens (1978), equation (3)]:

Description in cif file::

_atom_rho_orbital_radial_slater_n0         1 
_atom_rho_orbital_radial_slater_zeta0     22.77630  
_atom_rho_orbital_radial_slater_coeff_1s   0.95175 
_atom_rho_orbital_radial_slater_coeff_2s  -0.29620  
_atom_rho_orbital_radial_slater_coeff_3s   0.10614 
_atom_rho_orbital_radial_slater_coeff_4s  -0.02339
    """    
    MANDATORY_ATTRIBUTE = ("n0", "zeta0")
    OPTIONAL_ATTRIBUTE = ("coeff_1s", "coeff_2s", "coeff_3s", "coeff_4s", "coeff_5s",
                          "coeff_2p", "coeff_3p", "coeff_4p", "coeff_5p",
                          "coeff_3d", "coeff_4d", "coeff_5d",
                          "coeff_4f", "coeff_5f")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_rho_orbital_radial_slater"
    def __init__(self, n0=None, zeta0=None, 
                       coeff_1s=None, coeff_2s=None, coeff_3s=None, coeff_4s=None, coeff_5s=None,
                       coeff_2p=None, coeff_3p=None, coeff_4p=None, coeff_5p=None,
                       coeff_3d=None, coeff_4d=None, coeff_5d=None,
                       coeff_4f=None, coeff_5f=None):
        super(AtomRhoOrbitalRadialSlater, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                            optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                            internal_attribute=self.INTERNAL_ATTRIBUTE,
                                            prefix=self.PREFIX)

        self.n0 = n0 
        self.zeta0 = zeta0 
        self.coeff_1s = coeff_1s 
        self.coeff_2s = coeff_2s 
        self.coeff_3s = coeff_3s 
        self.coeff_4s = coeff_4s 
        self.coeff_5s = coeff_5s
        self.coeff_2p = coeff_2p 
        self.coeff_3p = coeff_3p 
        self.coeff_4p = coeff_4p 
        self.coeff_5p = coeff_5p
        self.coeff_3d = coeff_3d 
        self.coeff_4d = coeff_4d 
        self.coeff_5d = coeff_5d
        self.coeff_4f = coeff_4f 
        self.coeff_5f = coeff_5f

        if self.is_defined:
            self.form_object

    @property
    def n0(self)->str:
        return getattr(self, "__n0")
    @n0.setter
    def n0(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__n0", x_in)

    @property
    def zeta0(self)->str:
        """
Zeta are given in atomic units 

To transform it in inverse angstrems use:
    
    zeta0[Ang**-1] = zeta0[a.u.] / 0.52918
        """

        return getattr(self, "__zeta0")
    @zeta0.setter
    def zeta0(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__zeta0", x_in)

    @property
    def coeff_1s(self)->str:
        return getattr(self, "__coeff_1s")
    @coeff_1s.setter
    def coeff_1s(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_1s", x_in)

    @property
    def coeff_2s(self)->str:
        return getattr(self, "__coeff_2s")
    @coeff_2s.setter
    def coeff_2s(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_2s", x_in)

    @property
    def coeff_3s(self)->str:
        return getattr(self, "__coeff_3s")
    @coeff_3s.setter
    def coeff_3s(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_3s", x_in)

    @property
    def coeff_4s(self)->str:
        return getattr(self, "__coeff_4s")
    @coeff_4s.setter
    def coeff_4s(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_4s", x_in)

    @property
    def coeff_5s(self)->str:
        return getattr(self, "__coeff_5s")
    @coeff_5s.setter
    def coeff_5s(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_5s", x_in)

    @property
    def coeff_2p(self)->str:
        return getattr(self, "__coeff_2p")
    @coeff_2p.setter
    def coeff_2p(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_2p", x_in)

    @property
    def coeff_3p(self)->str:
        return getattr(self, "__coeff_3p")
    @coeff_3p.setter
    def coeff_3p(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_3p", x_in)

    @property
    def coeff_4p(self)->str:
        return getattr(self, "__coeff_4p")
    @coeff_4p.setter
    def coeff_4p(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_4p", x_in)

    @property
    def coeff_5p(self)->str:
        return getattr(self, "__coeff_5p")
    @coeff_5p.setter
    def coeff_5p(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_5p", x_in)

    @property
    def coeff_3d(self)->str:
        return getattr(self, "__coeff_3d")
    @coeff_3d.setter
    def coeff_3d(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_3d", x_in)

    @property
    def coeff_4d(self)->str:
        return getattr(self, "__coeff_4d")
    @coeff_4d.setter
    def coeff_4d(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_4d", x_in)

    @property
    def coeff_5d(self)->str:
        return getattr(self, "__coeff_5d")
    @coeff_5d.setter
    def coeff_5d(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_5d", x_in)

    @property
    def coeff_4f(self)->str:
        return getattr(self, "__coeff_4f")
    @coeff_4f.setter
    def coeff_4f(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_4f", x_in)

    @property
    def coeff_5f(self)->str:
        return getattr(self, "__coeff_5f")
    @coeff_5f.setter
    def coeff_5f(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__coeff_5f", x_in)

    def calc_normalized_rho(self, radius: numpy.array, kappa=1.) -> numpy.array:
        """
Gives normalized Slatter type function at given radius and defined n0 and zeta0:
    
    R = ((2*zeta0)**(n0+0.5)/(2*n0)!)*r**(n0-1)*exp(-zeta0*r) 
    (see Clementi_1974)

Input arguments:

    radius is numpy 1D array of float numbers in Angstrem
    kappa is describes contraction/extension of atomic orbital

Output arguments:
    normalized radial density given as numpy 1D array of float numbers

Example:

    >>> radius = numpy.linspace(0.,5,600)
    >>> rho = obj.calc_normalized_rho(radius)
        """
        #transformation from atomic units to inverse angstrems
        zeta0 = self.zeta0
        n0 = self.n0
        zeta_ang = kappa*zeta0/0.52918
        coeff_norm=((2.*zeta_ang)**(float(n0)+0.5))/math.sqrt(math.factorial(2*n0))
        norm_density = coeff_norm*(radius**(n0-1))*numpy.exp(-zeta_ang*radius)
        return norm_density

    @property
    def is_variable(self):
        _l = []
        if self.zeta0 is not None: _l.append(self.zeta0.refinement)
        if self.coeff_1s is not None: _l.append(self.coeff_1s.refinement)
        if self.coeff_2s is not None: _l.append(self.coeff_2s.refinement)
        if self.coeff_3s is not None: _l.append(self.coeff_3s.refinement)
        if self.coeff_4s is not None: _l.append(self.coeff_4s.refinement)
        if self.coeff_5s is not None: _l.append(self.coeff_5s.refinement)
        if self.coeff_2p is not None: _l.append(self.coeff_2p.refinement)
        if self.coeff_3p is not None: _l.append(self.coeff_3p.refinement)
        if self.coeff_4p is not None: _l.append(self.coeff_4p.refinement)
        if self.coeff_5p is not None: _l.append(self.coeff_5p.refinement)
        if self.coeff_3d is not None: _l.append(self.coeff_3d.refinement)
        if self.coeff_4d is not None: _l.append(self.coeff_4d.refinement)
        if self.coeff_5d is not None: _l.append(self.coeff_5d.refinement)
        if self.coeff_4f is not None: _l.append(self.coeff_4f.refinement)
        if self.coeff_5f is not None: _l.append(self.coeff_5f.refinement)
        res = any(_l)
        return res

    def get_variables(self):
        l_variable = []
        if self.zeta0 is not None:
            if self.zeta0.refinement: l_variable.append(self.zeta0)
        if self.coeff_1s is not None:
            if self.coeff_1s.refinement: l_variable.append(self.coeff_1s)
        if self.coeff_2s is not None:
            if self.coeff_2s.refinement: l_variable.append(self.coeff_2s)
        if self.coeff_3s is not None:
            if self.coeff_3s.refinement: l_variable.append(self.coeff_3s)
        if self.coeff_4s is not None:
            if self.coeff_4s.refinement: l_variable.append(self.coeff_4s)
        if self.coeff_5s is not None:
            if self.coeff_5s.refinement: l_variable.append(self.coeff_5s)
        if self.coeff_2p is not None:
            if self.coeff_2p.refinement: l_variable.append(self.coeff_2p)
        if self.coeff_3p is not None:
            if self.coeff_3p.refinement: l_variable.append(self.coeff_3p)
        if self.coeff_4p is not None:
            if self.coeff_4p.refinement: l_variable.append(self.coeff_4p)
        if self.coeff_5p is not None:
            if self.coeff_5p.refinement: l_variable.append(self.coeff_5p)
        if self.coeff_3d is not None:
            if self.coeff_3d.refinement: l_variable.append(self.coeff_3d)
        if self.coeff_4d is not None:
            if self.coeff_4d.refinement: l_variable.append(self.coeff_4d)
        if self.coeff_5d is not None:
            if self.coeff_5d.refinement: l_variable.append(self.coeff_5d)
        if self.coeff_4f is not None:
            if self.coeff_4f.refinement: l_variable.append(self.coeff_4f)
        if self.coeff_5f is not None:
            if self.coeff_5f.refinement: l_variable.append(self.coeff_5f)

        return l_variable


class AtomRhoOrbitalRadialSlaterL(LoopConstr):
    """
These items are used when the radial dependence of the
electron density, R(k(l), l, r), of the atom specified in
_atom_rho_orbital_atom_label is expressed as a Slater-type
function [Hansen & Coppens (1978), equation (3)]:

Description in cif file::

loop_V_s
_atom_rho_orbital_radial_slater_n0
_atom_rho_orbital_radial_slater_zeta0
_atom_rho_orbital_radial_slater_coeff_1s
_atom_rho_orbital_radial_slater_coeff_2s
_atom_rho_orbital_radial_slater_coeff_3s
_atom_rho_orbital_radial_slater_coeff_4s
    1 22.77630  0.95175 -0.29620  0.10614 -0.02339
    1 36.05340  0.02115  0.00019  0.00042 -0.00046
    2 19.54100  0.03441 -0.15795  0.06539 -0.01626
    2  9.37400  0.00415 10.04067 -0.42093  0.09958
    3  7.90503 -0.00366  0.12576 -0.25395  0.05408
    3  5.12985  0.00552 -0.02332  0.73273 -0.11001
    """
    CATEGORY_KEY = ()
    ITEM_CLASS = AtomRhoOrbitalRadialSlater
    def __init__(self, item=[], loop_name=""):
        super(AtomRhoOrbitalRadialSlaterL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def calc_normalized_rho(self, radius: numpy.array, shell: str, kappa=1., *argv) -> numpy.array:
        """
Gives normalized Slatter type function at given radius for given shell:
    
    R = sum_i coeff_i * ((2*zeta0)**(n0+0.5)/(2*n0)!)*r**(n0-1)*exp(-zeta0*r) 
    (see Clementi_1974)

Input arguments:

    radius is numpy 1D array of float numbers in Angstrem
    shell is it is orbital shell
    kappa is describes contraction/extension of atomic orbital

Output arguments:
    normalized radial density given as numpy 1D array of float numbers

Example:

    >>> radius = numpy.linspace(0.,5,600)
    >>> rho = obj.calc_rho(radius, "3d")
    >>> rho = obj.calc_rho(radius, "3d", "2s")
        """
        #transformation from atomic units to inverse angstrems
        l_shell = [shell]
        l_shell.extend(argv)
        n_shell = len(l_shell)
        n_item = len(self.item)
        coeff = numpy.zeros((n_shell, n_item), dtype=float)

        norm_density = numpy.zeros((radius.size, n_shell), dtype=float)

        l_attr = [f"coeff_{_.strip():}" for _ in l_shell]
        for _i, _item in enumerate(self.item):
            normalized_rho = _item.calc_normalized_rho(radius, kappa=kappa)
            for _j, _attr in enumerate(l_attr):
                coeff[_j, _i] = float(getattr(_item, _attr))
            norm_density += normalized_rho[:, numpy.newaxis] * (coeff[:, _i])[numpy.newaxis, :]
        return norm_density

    def take_objects_for_atom_type(self, atom_type: str) -> list:
        cls = type(self)
        l_arors = []
        f_dir = os.path.dirname(__file__)
        f_name = os.path.join(f_dir, "library.rcif")

        obj_rcif = to_data(f_name)
        s_atom_type = ("".join([_ for _ in atom_type if _.isalpha()])).lower()
        for loop in obj_rcif.loops:
            l_loop_rcif = cls.from_cif(str(loop))
            if l_loop_rcif is not None:
                loop_rcif = l_loop_rcif[0]
                s_loop_name = loop_rcif.loop_name
                s_atom_type_loop = (s_loop_name.split("_")[0]).lower()
                if s_atom_type_loop == s_atom_type:
                    l_arors.append(loop_rcif)
        return l_arors

    def calc_jl_by_radial_density(self, sthovl, lmax: int, shell: str, kappa=1.):
        n0 = numpy.array(self.n0, dtype=int)
        zeta0 = numpy.array(self.zeta0, dtype=float)
        coeff_h = getattr(self, f"coeff_{shell:}")

        coeff = numpy.array([_.value for _ in coeff_h], dtype=float)

        if any([_ is None for _ in coeff]): return None
        
        jl = calc_GCF(n0, zeta0, coeff, kappa, 
                      n0, zeta0, coeff, kappa, 
                      sthovl, lmax)
        return jl

    def refine_coefficients_by_jl(self, atom_type: str, shell: str, 
                                  sthovl_min=0.0, sthovl_max=2.0, lande=2., kappa=1., 
                                  sthvl_num=100):
        ats = AtomTypeScat.form_by_symbol(atom_type)

        sthovl = numpy.linspace(sthovl_min, sthovl_max, num=int(sthvl_num))
        ff = ats.calc_form_factor(sthovl, lande=lande, kappa=kappa)
        fitable = self.get_variables()

        def tempfunc(param):
            for _1, _2 in zip(fitable, param):
                _1.value = _2 
            j0 = self.calc_jl_by_radial_density(sthovl, 0, shell, kappa)[0]
            chi_sq = (numpy.square(ff - j0)).sum()
            return chi_sq

        param_0 = [_.value for _ in fitable]
        j0 = self.calc_jl_by_radial_density(sthovl, 0, shell, kappa)[0]

        res = scipy.optimize.basinhopping(tempfunc, param_0, niter=10, 
                                          T=0.1, stepsize=0.1, interval=20, 
                                          disp=True)
        return res
