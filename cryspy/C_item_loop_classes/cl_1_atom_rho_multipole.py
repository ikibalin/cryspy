

import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.charge_form_factor import get_n_zeta_coeff_for_atom_with_orbital, get_atom_name_ion_charge_shell

from cryspy.A_functions_base.orbital_functions import calc_jl

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

import matplotlib.pyplot as plt

PARAM_PLM = [
            "coeff_P00", "coeff_P11", "coeff_P1m1", "coeff_P10", 
            "coeff_P20", "coeff_P21", "coeff_P2m1", "coeff_P22", "coeff_P2m2",
            "coeff_P30", "coeff_P31", "coeff_P3m1", "coeff_P32", "coeff_P3m2", "coeff_P33", "coeff_P3m3",
            "coeff_P40", "coeff_P41", "coeff_P4m1", "coeff_P42", "coeff_P4m2", "coeff_P43", "coeff_P4m3", "coeff_P44", "coeff_P4m4"]

INDEX_PLM = numpy.array([[0,0],[1,1],[1,-1],[1,0],
                     [2,0],[2,1],[2,-1],[2,2],[2,-2],
                     [3,0],[3,1],[3,-1],[3,2],[3,-2],[3,3],[3,-3],
                     [4,0],[4,1],[4,-1],[4,2],[4,-2],[4,3],[4,-3],[4,4],[4,-4]],dtype=int)

class AtomRhoMultipole(ItemN):
    """   This category contains information about the multipole
   coefficients used to describe the electron density.

   This definition is not the standard one defined for Electron density CIF: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Catom_rho_multipole.html

    Mandatory attributes:
        - atom_label
        - atom_type_symbol
        - orbital

    Optional attributes:
        - kappa
        - P00
        - P11, P1-1, P10
        - P20, P21, P2-1, P22, P2-2
        - P30, P31, P3-1, P32, P3-2, P33, P3-3
        - P40, P41, P4-1, P42, P4-2, P43, P4-3, P44, P4-4
    """
    ATTR_MANDATORY_NAMES = ("atom_label", "atom_type_symbol", "orbital")
    ATTR_MANDATORY_TYPES = (str, str, str)
    ATTR_MANDATORY_CIF = ("atom_label", "atom_type_symbol", "orbital", )

    ATTR_OPTIONAL_NAMES = (
        "kappa", "coeff_P00",
        "coeff_P11", "coeff_P1m1", "coeff_P10",
        "coeff_P20", "coeff_P21", "coeff_P2m1", "coeff_P22", "coeff_P2m2",
        "coeff_P30", "coeff_P31", "coeff_P3m1", "coeff_P32", "coeff_P3m2", "coeff_P33", "coeff_P3m3",
        "coeff_P40", "coeff_P41", "coeff_P4m1", "coeff_P42", "coeff_P4m2", "coeff_P43", "coeff_P4m3", "coeff_P44", "coeff_P4m4",
        )
    ATTR_OPTIONAL_TYPES = (
        float, float,
        float, float, float, 
        float, float, float, float, float, 
        float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float,
        )
    ATTR_OPTIONAL_CIF = (
        "kappa", "coeff_P00",
        "coeff_P11", "coeff_P1-1", "coeff_P10",
        "coeff_P20", "coeff_P21", "coeff_P2-1", "coeff_P22", "coeff_P2-2",
        "coeff_P30", "coeff_P31", "coeff_P3-1", "coeff_P32", "coeff_P3-2", "coeff_P33", "coeff_P3-3",
        "coeff_P40", "coeff_P41", "coeff_P4-1", "coeff_P42", "coeff_P4-2", "coeff_P43", "coeff_P4-3", "coeff_P44", "coeff_P4-4")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = (
        "kappa", "coeff_P00",
        "coeff_P11", "coeff_P1m1", "coeff_P10",
        "coeff_P20", "coeff_P21", "coeff_P2m1", "coeff_P22", "coeff_P2m2",
        "coeff_P30", "coeff_P31", "coeff_P3m1", "coeff_P32", "coeff_P3m2", "coeff_P33", "coeff_P3m3",
        "coeff_P40", "coeff_P41", "coeff_P4m1", "coeff_P42", "coeff_P4m2", "coeff_P43", "coeff_P4m3", "coeff_P44", "coeff_P4m4",
        )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {
        "kappa": 1.,
        "coeff_P00": 0.,
        }
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_rho_multipole"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomRhoMultipole, self).__init__()

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

    def get_plm_lm(self):



        p_lm = numpy.zeros(shape=(25,), dtype=float)
        for i_param, param in enumerate(PARAM_PLM):
            if self.is_attribute(param):
                p_lm[i_param] = float(getattr(self, param))
        return p_lm, INDEX_PLM 

    def load_plm_from_dictionary(self, np_plm):

        for i_param, param in enumerate(PARAM_PLM):
            if self.is_attribute(param):
                setattr(self, param, float(np_plm[i_param]))

    def set_plm_sigma_by_index(self, ind_p, sigma):

        if ind_p < len(PARAM_PLM):
            param = PARAM_PLM[ind_p]
            if self.is_attribute(param):
                setattr(self, param+"_sigma", float(sigma))

    def get_flags_plm(self):
        flags_plm = numpy.zeros(shape=(25,), dtype=bool)
        for i_param, param in enumerate(PARAM_PLM):
            if self.is_attribute(param):
                flags_plm[i_param] = getattr(self, param+"_refinement") 
        return flags_plm

    def get_n_zeta_coeff(self):
        atom_name = get_atom_name_ion_charge_shell(self.atom_type_symbol)[0]
        n, zeta, coeff = get_n_zeta_coeff_for_atom_with_orbital(atom_name, self.orbital)
        return n, zeta, coeff

    def calc_jl(self, sthovl, l_max = 4):
        n, zeta, coeff = self.get_n_zeta_coeff()
        kappa = self.kappa
        jl = calc_jl(sthovl, coeff, n, zeta, kappa=kappa, l_max = l_max).transpose()
        return jl

    def calc_jl_p(self, sthovl, l_max = 4):
        jl = self.calc_jl(sthovl, l_max = l_max)
        plm, lm = self.get_plm_lm()
        l_res = []
        for p, (l, m) in zip(plm, lm):
            if numpy.isclose(p, 0):
                continue
            l_res.append( ((l,m),numpy.power(1j,l).real * p * jl[int(l)] ))
        return l_res

    def plots(self):
        return [self.plot_form_factor(),]

    def plot_form_factor(self):
        """Plot magnetic form factor.
        """
        x_min, x_max = 0, 1.5
        sthovl = numpy.linspace(x_min, x_max, 100)
        try:
            l_res = self.calc_jl_p(sthovl, l_max = 4)
        except AttributeError as e:
            print(f"Error occurred while calculating form factor for {self.label}: {e}")
            return
        label = self.atom_label
        fig, ax = plt.subplots()
        ax.set_xlabel('sin theta / lambda')
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel('Magnetic form factor')
        for (l, m), res in l_res:
            ax.plot(sthovl, res, label=label+f" {l},{m}") 
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)
    
class AtomRhoMultipoleL(LoopN):
    """Details about properties of the atoms.
    """
    ITEM_CLASS = AtomRhoMultipole
    ATTR_INDEX = "atom_label"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomRhoMultipoleL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
    def get_dictionary(self):
        d_out = super(AtomRhoMultipoleL, self).get_dictionary()
        np_plm, lm  = self.get_plm_lm()
        flags_np_plm = self.get_flags_plm()
        flags_kappa = self.get_flags_kappa()
        d_out["atom_rho_multipole_plm"] = np_plm
        d_out["flags_atom_rho_multipole_plm"] = flags_np_plm
        d_out["atom_rho_multipole_lm"] = lm
        d_out["atom_rho_multipole_kappa"] = numpy.array(self.kappa, dtype=float)
        d_out["flags_atom_rho_multipole_kappa"] = flags_kappa
        d_out["atom_rho_multipole_orbital"] = numpy.array(self.orbital, dtype=str)
        d_out["atom_rho_multipole_label"] = numpy.array(self.atom_label, dtype=str)
        d_out["atom_rho_multipole_type_symbol"] = numpy.array(self.atom_type_symbol, dtype=str)
        l_n, l_zeta, l_coeff = [],[],[]
        for item in self.items:
            n, zeta, coeff = item.get_n_zeta_coeff()
            l_n.append(n)
            l_zeta.append(zeta)
            l_coeff.append(coeff)
        d_out["atom_rho_multipole_n_zeta_coeff"] = (l_n, l_zeta, l_coeff)
        return d_out

    def get_plm_lm(self):

        l_p_lm = []
        for item_arm in self.items:
            plm, lm = item_arm.get_plm_lm()
            l_p_lm.append(plm)
        np_plm = numpy.stack(l_p_lm, axis=1)
        return np_plm, lm 

    def get_flags_plm(self):

        l_flags_plm = []
        for item_arm in self.items:
            flags_plm = item_arm.get_flags_plm()
            l_flags_plm.append(flags_plm)
        np_flags_plm = numpy.stack(l_flags_plm, axis=1)
        return np_flags_plm
    def get_flags_kappa(self):

        l_flags_kappa = []
        for item_arm in self.items:
            flags_kappa = item_arm.kappa_refinement
            l_flags_kappa.append(flags_kappa)
        np_flags_kappa = numpy.array(l_flags_kappa, dtype=bool)
        return np_flags_kappa
    def plots(self):
        return [self.plot_form_factors(),]

    def plot_form_factors(self):
        x_min, x_max = 0, 1.5
        sthovl = numpy.linspace(x_min, x_max, 100)
        fig, ax = plt.subplots()
        ax.set_xlabel('sin theta / lambda')
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel('Magnetic form factor )')

        for item in self.items:
            try:
                label = item.atom_label
                l_res = item.calc_jl_p(sthovl, l_max = 4)
            except AttributeError as e:
                print(f"Error occurred while calculating form factor for {self.atom_label}: {e}")
                continue
            for (l, m), res in l_res:
                ax.plot(sthovl, res, label=label+f" j_{l},{m} / (4*pi)") 
        
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)
        

# s_cont = """
#  loop_
#  _atom_type_symbol
#  _atom_type_oxidation_number
#  _atom_type_number_in_cell
#  _atom_type_scat_dispersion_real
#  _atom_type_scat_dispersion_imag
#  _atom_type_scat_source
#    C 0 72  .017  .009  International_Tables_Vol_IV_Table_2.2B
#    H 0 100  0     0    International_Tables_Vol_IV_Table_2.2B
#    O 0 12  .047  .032  International_Tables_Vol_IV_Table_2.2B
#    N 0 4   .029  .018  International_Tables_Vol_IV_Table_2.2B
#   """

# obj = AtomTypeL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["O"].scat_length_neutron, end="\n\n")
