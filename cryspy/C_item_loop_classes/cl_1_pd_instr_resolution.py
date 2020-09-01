import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class PdInstrResolution(ItemN):
    """
    PdInstrReflexAsymmetry describes asymmetry of Bragg reflections for
    1d powder diffractometer.

    Mandatory attributes:
        - ub_11, ub_12, ub_13, ub_21, ub_22, ub_23, ub_31, ub_32, ub_33

    Optional attributes:
        - occupancy
        - adp_type
        - u_iso_or_equiv
        - u_equiv_geom_mean
        - b_iso_or_equiv
        - multiplicity
        - wyckoff_symbol
        - cartn_x
        - cartn_y
        - cartn_z

    Internal attributes:
        - scat_length_neutron

    Internal protected attributes:
        - space_group_wyckoff
        - constr_number
    """
    ATTR_MANDATORY_NAMES = ("u", "v", "w", "x", "y")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float)
    ATTR_MANDATORY_CIF = ("U", "V", "W", "X", "Y")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("u", "v", "w", "x", "y")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"u": "{:.5f}", "v": "{:.5f}", "w": "{:.5f}", "x": "{:.2f}",
                 "y": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"u": 0., "v": 0., "w": 0., "x": 0., "y": 0.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "pd_instr_resolution"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdInstrResolution, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_tancos(self, th_hkl):
        """
        tth_hkl in radianas
        calculate tangenth (theta)
        """
        self.t_th = numpy.tan(th_hkl)
        self.t_th_sq = self.t_th**2
        res = numpy.cos(th_hkl)

        self.c_th = res
        self.ic_th = 1./res
        
    def calc_hg(self, i_g = 0.):
        """
        ttheta in radians, could be array
        gauss size
        """
        u, v, w = float(self.u), float(self.v), float(self.w)
        res_sq = (u*self.t_th_sq + v*self.t_th + w + 
                  i_g*self.ic_th**2)
        self.hg = numpy.sqrt(res_sq)
        
    def calc_hl(self):
        """
        ttheta in radians, could be array
        lorentz site
        """
        x, y = float(self.x), float(self.y)
        self.hl = x*self.t_th + y*self.ic_th

    def calc_hpveta(self):
        """
        ttheta in radians, could be array
        pseudo-Voight function
        """
        hg = self.hg
        hl = self.hl

        hg_2, hl_2 = hg**2, hl**2
        hg_3, hl_3 = hg_2*hg, hl_2*hl
        hg_4, hl_4 = hg_3*hg, hl_3*hl
        hg_5, hl_5 = hg_4*hg, hl_4*hl
        c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
        hpv = (hg_5 + c_2*hg_4*hl + c_3*hg_3*hl_2 + 
                       c_4*hg_2*hl_3 + c_5*hg*hl_4 + hl_5)**0.2
        hh = hl*1./hpv
        self.hpv = hpv 
        self.eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3

    def calc_agbg(self):
        hpv = self.hpv

        self.ag = (2./hpv)*(numpy.log(2.)/numpy.pi)**0.5
        self.bg = 4*numpy.log(2)/(hpv**2)
        
    def calc_albl(self):
        hpv = self.hpv
        self.al = 2./(numpy.pi*hpv )
        self.bl = 4./(hpv**2)
    
    def calc_resolution(self, tth_hkl, i_g = 0.):
        """
Calculate parameters for tth
tth_hkl in degrees

Output values:

h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l
        """
        self.calc_tancos(0.5*tth_hkl*numpy.pi/180.)
        self.calc_hg(i_g = i_g)
        self.calc_hl()
        self.calc_hpveta()
        self.calc_agbg()
        self.calc_albl()

        a_g = self.ag  
        b_g = self.bg  
        a_l = self.al  
        b_l = self.bl 
        h_g = self.hg  
        h_l = self.hl 
        h_pv = self.hpv 
        eta = self.eta
        
        return h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l
    


class PdInstrResolutionL(LoopN):
    """
    Description of AtomSite in loop.

    """
    ITEM_CLASS = PdInstrResolution
    ATTR_INDEX = None
    def __init__(self, loop_name = None) -> NoReturn:
        super(PdInstrResolutionL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
#  loop_

#  _pd_instr_resolution_u 
#  _pd_instr_resolution_v 
#  _pd_instr_resolution_w
#  _pd_instr_resolution_x 
#  _pd_instr_resolution_y 
# 16.9776(25) -2.8357 0.5763 0. 0.
# """

# obj = PdInstrResolutionL.from_cif(s_cont)
# print(obj, end="\n\n")
