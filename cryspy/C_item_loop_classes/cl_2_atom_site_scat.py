"""Classes AtomSiteScat, AtomSiteScatL."""
from typing import NoReturn
import numpy
import matplotlib.pyplot as plt

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_atom_type_scat import AtomTypeScat, AtomTypeScatL


class AtomSiteScat(ItemN):
    """Description of magnetic structure factor.

    Attributes
    ----------
        - label (mandatory)
        - lande, kappa, type_symbol (optional)
        - atom_type_scat (internal protected attribute)

    Method
    ------
        - calc_form_factor(sthovl, flag_only_orbital=False)
        - load_atom_type_scat_by_symbol(symbol:str)
    """

    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = ("type_symbol", "lande", "kappa")
    ATTR_OPTIONAL_TYPES = (str, float, float)
    ATTR_OPTIONAL_CIF = ("type_symbol", "Lande", "kappa")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ("atom_type_scat", )

    # parameters considered are refined parameters
    ATTR_REF = ("lande", "kappa")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {"modulation_flag": ["yes", "y", "no", "n"],
                     "refinement_flags_magnetic": ["S", "M", "A", "SM", "SA",
                                                   "MA", "SMA"]}

    # default values for the parameters
    D_DEFAULT = {"lande": 2., "kappa": 1.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_scat"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteScat, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"lande": 0, "kappa": 0}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_form_factor(self, sthovl, flag_only_orbital=False):
        """Calculate form factor."""
        if self.is_attribute("atom_type_scat"):
            atom_type_scat = self.atom_type_scat
        else:
            self.load_atom_type_scat_by_symbol(self.type_symbol)
            atom_type_scat = self.atom_type_scat

        form_factor = atom_type_scat.calc_form_factor(
            sthovl, lande=self.lande, kappa=self.kappa,
            flag_only_orbital=flag_only_orbital)
        return form_factor

    def load_atom_type_scat_by_symbol(self, symbol: str):
        """Load details about atom type scattering."""
        self.type_symbol = symbol
        _a_t_s = AtomTypeScat.form_by_symbol(symbol)
        self.__dict__["atom_type_scat"] = _a_t_s

    def plots(self):
        return [self.plot_form_factor()]

    def plot_form_factor(self):
        """Plot magnetic form factor.
        """
        x_min, x_max = 0, 1.5
        sthovl = numpy.linspace(x_min, x_max, 100)
        try:
            res = self.calc_form_factor(sthovl, flag_only_orbital=False)
        except AttributeError:
            return
        label = self.label
        fig, ax = plt.subplots()
        ax.set_xlabel('sin theta / lambda')
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel('Magnetic form factor')
        ax.plot(sthovl, res, label=label) 
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def get_flags_lande(self):
        res = self.lande_refinement
        return res

    def get_flags_kappa(self):
        res = self.kappa_refinement
        return res

class AtomSiteScatL(LoopN):
    """
    Description of magnetic structure factor by Lande factor and kappa.

    Methods
    -------
        - calc_form_factor(sthovl, flag_only_orbital=False)
        - load_atom_type_scat_by_atom_site(atom_site)

    """

    ITEM_CLASS = AtomSiteScat
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteScatL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def calc_form_factor(self, sthovl, flag_only_orbital=False):
        """Calculate form factor."""
        form_factor = [item.calc_form_factor(
            sthovl, flag_only_orbital=flag_only_orbital)
            for item in self.items]
        return numpy.array(list(zip(*form_factor)), dtype=float)

    def load_atom_type_scat_by_atom_site(self, atom_site) -> NoReturn:
        """Load details about atom type scattering."""
        all([item.load_atom_type_scat_by_symbol(
            atom_site[item.label].type_symbol) for item in self.items])
   
    def plots(self):
        return [self.plot_form_factor()]

    def plot_form_factor(self):
        """Plot magnetic form factor.
        """
        x_min, x_max = 0, 1.5
        sthovl = numpy.linspace(x_min, x_max, 100)
        try:
            res = self.calc_form_factor(sthovl, flag_only_orbital=False)
        except AttributeError:
            return
        labels = self.label
        fig, ax = plt.subplots()
        ax.set_title("Magnetic form factor")
        ax.set_xlabel('sin theta / lambda')
        ax.set_xlim(x_min, x_max)
        ax.set_ylabel('Magnetic form factor')
        for ind in range(res.shape[1]):
            ax.plot(sthovl, res[:, ind], label=labels[ind], alpha=0.7) 
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def report(self):
        s_out = ""
        l_ats = []
        for item in self.items:
            try:
                ats = item.atom_type_scat
                l_ats.append(ats)
            except AttributeError:
                pass
        if len(l_ats) != 0:
            obj = AtomTypeScatL()
            obj.items=l_ats
            s_out = str(obj)
        return s_out
    
    def get_flags_lande(self):
        flags_lande = numpy.stack([item.get_flags_lande() for item in self.items], axis=0)
        return flags_lande

    def get_flags_kappa(self):
        flags_kappa = numpy.stack([item.get_flags_kappa() for item in self.items], axis=0)
        return flags_kappa

# s_cont = """
#  loop_
#  _atom_site_scat_label
#  _atom_site_scat_Lande
#  _atom_site_scat_kappa
#   Fe3A    2.00 1.00
#   Fe3B    2.00 1.00(12)
#   """

# obj = AtomSiteScatL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["Fe3B"], end="\n\n")
