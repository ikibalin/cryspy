"""Described classes AtomTypeScat, AtomTypeScatL."""
import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_magnetic import \
    get_j0_j2_by_symbol
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

class AtomTypeScat(ItemN):
    """
    AtomTypeScat class.

    Data items in the ATOM_TYPE_SCAT category describe atomic
    scattering information used in crystallographic structure studies.
    This category is fully defined in the core CIF dictionary.

    Attributes
    ----------
        - symbol (mandatory)
        - neutron_magnetic_j0_a1, neutron_magnetic_j0_a2
        - neutron_magnetic_j0_b1, neutron_magnetic_j0_b2
        - neutron_magnetic_j0_c1, neutron_magnetic_j0_c2
        - neutron_magnetic_j0_d
        - neutron_magnetic_j2_a1, neutron_magnetic_j2_a2
        - neutron_magnetic_j2_b1, neutron_magnetic_j2_b2
        - neutron_magnetic_j2_c1, neutron_magnetic_j2_c2
        - neutron_magnetic_j2_d
        - neutron_magnetic_j4_a1, neutron_magnetic_j4_a2
        - neutron_magnetic_j4_b1, neutron_magnetic_j4_b2
        - neutron_magnetic_j4_c1, neutron_magnetic_j4_c2
        - neutron_magnetic_j4_d
        - neutron_magnetic_j6_a1, neutron_magnetic_j6_a2
        - neutron_magnetic_j6_b1, neutron_magnetic_j6_b2
        - neutron_magnetic_j6_c1, neutron_magnetic_j6_c2
        - neutron_magnetic_j6_d

    Methods
    -------
        - calc_form_factor
        - calc_j0
        - calc_j2
        - load_by_symbol
        - form_by_symbol (classmethod)
    """

    ATTR_MANDATORY_NAMES = ("symbol", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("symbol", )

    ATTR_OPTIONAL_NAMES = ("neutron_magnetic_j0_a1", "neutron_magnetic_j0_a2",
                           "neutron_magnetic_j0_b1", "neutron_magnetic_j0_b2",
                           "neutron_magnetic_j0_c1", "neutron_magnetic_j0_c2",
                           "neutron_magnetic_j0_d",  "neutron_magnetic_j0_e",
                           "neutron_magnetic_j2_a1", "neutron_magnetic_j2_a2",
                           "neutron_magnetic_j2_b1", "neutron_magnetic_j2_b2",
                           "neutron_magnetic_j2_c1", "neutron_magnetic_j2_c2",
                           "neutron_magnetic_j2_d",  "neutron_magnetic_j2_e",
                           "neutron_magnetic_j4_a1", "neutron_magnetic_j4_a2",
                           "neutron_magnetic_j4_b1", "neutron_magnetic_j4_b2",
                           "neutron_magnetic_j4_c1", "neutron_magnetic_j4_c2",
                           "neutron_magnetic_j4_d",  "neutron_magnetic_j4_e",
                           "neutron_magnetic_j6_a1", "neutron_magnetic_j6_a2",
                           "neutron_magnetic_j6_b1", "neutron_magnetic_j6_b2",
                           "neutron_magnetic_j6_c1", "neutron_magnetic_j6_c2",
                           "neutron_magnetic_j6_d",  "neutron_magnetic_j6_e",
                           "neutron_magnetic_source")
    ATTR_OPTIONAL_TYPES = (
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, float, float, float,
        float, float, str)
    ATTR_OPTIONAL_CIF = ("neutron_magnetic_j0_a1", "neutron_magnetic_j0_a2",
                         "neutron_magnetic_j0_b1", "neutron_magnetic_j0_b2",
                         "neutron_magnetic_j0_c1", "neutron_magnetic_j0_c2",
                         "neutron_magnetic_j0_d",  "neutron_magnetic_j0_e",
                         "neutron_magnetic_j2_a1", "neutron_magnetic_j2_a2",
                         "neutron_magnetic_j2_b1", "neutron_magnetic_j2_b2",
                         "neutron_magnetic_j2_c1", "neutron_magnetic_j2_c2",
                         "neutron_magnetic_j2_d",  "neutron_magnetic_j2_e",
                         "neutron_magnetic_j4_a1", "neutron_magnetic_j4_a2",
                         "neutron_magnetic_j4_b1", "neutron_magnetic_j4_b2",
                         "neutron_magnetic_j4_c1", "neutron_magnetic_j4_c2",
                         "neutron_magnetic_j4_d",  "neutron_magnetic_j4_e",
                         "neutron_magnetic_j6_a1", "neutron_magnetic_j6_a2",
                         "neutron_magnetic_j6_b1", "neutron_magnetic_j6_b2",
                         "neutron_magnetic_j6_c1", "neutron_magnetic_j6_c2",
                         "neutron_magnetic_j6_d",  "neutron_magnetic_j6_e",
                         "neutron_magnetic_source")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("j0_parameters", "j2_parameters")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_type_scat"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomTypeScat, self).__init__()

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

    def calc_form_factor(self, sthovl, lande=2., kappa=1.,
                         flag_only_orbital=False):
        r"""Calculate magnetic form factor in frame of Spherical model.

        (Int.Tabl.C.p.592)

        :LFactor:  is Lande factor
        :lsthovl:  is list :math:`sin(\\theta)/\\lambda` in :math:`\\A^{-1}`

        Calculation of magnetic form factor <j0>, <j2>, <j4>, <j6>

        :coeff: is a list [A, a, B, b, C, c, D] at n = 0, 2, 4, 6

        `<https://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.140405>`
        mismatch with international tables where (1.0-2.0/np_factor_lande)
        """
        # not sure about kappa, it is here just for test, by default it is 1.0
        j2_av = self.calc_j2(sthovl, kappa=kappa)

        if flag_only_orbital:
            form_factor = (2.0/float(lande)-1.0)*j2_av
        else:
            j0_av = self.calc_j0(sthovl, kappa=kappa)
            form_factor = j0_av+(2.0/float(lande)-1.0)*j2_av
        return form_factor

    def calc_j0(self, sthovl, kappa: float = 1):
        """Calculate j0."""
        j0_A = self.neutron_magnetic_j0_a1
        j0_a = self.neutron_magnetic_j0_a2
        j0_B = self.neutron_magnetic_j0_b1
        j0_b = self.neutron_magnetic_j0_b2
        j0_C = self.neutron_magnetic_j0_c1
        j0_c = self.neutron_magnetic_j0_c2
        j0_D = self.neutron_magnetic_j0_d
        _h = (sthovl/float(kappa))**2
        j0_av = (j0_A*numpy.exp(-j0_a*_h) +
                 j0_B*numpy.exp(-j0_b*_h) +
                 j0_C*numpy.exp(-j0_c*_h)+j0_D)
        return j0_av

    def calc_j2(self, sthovl, kappa: float = 1):
        """Calculate j2."""
        j2_A = self.neutron_magnetic_j2_a1
        j2_a = self.neutron_magnetic_j2_a2
        j2_B = self.neutron_magnetic_j2_b1
        j2_b = self.neutron_magnetic_j2_b2
        j2_C = self.neutron_magnetic_j2_c1
        j2_c = self.neutron_magnetic_j2_c2
        j2_D = self.neutron_magnetic_j2_d
        _h = (sthovl/float(kappa))**2
        j2_av = (j2_A*numpy.exp(-j2_a*_h) +
                 j2_B*numpy.exp(-j2_b*_h) +
                 j2_C*numpy.exp(-j2_c*_h) + j2_D)*_h
        return j2_av

    def load_by_symbol(self, symbol: str = None) -> NoReturn:
        """Load by symbol."""
        if symbol is None:
            symbol = self.symbol
        else:
            self.symbol = symbol
        j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, j2_B1, \
            j2_b2, j2_C1, j2_c2, j2_D = get_j0_j2_by_symbol(symbol)
        self.neutron_magnetic_j0_a1 = j0_A1
        self.neutron_magnetic_j0_a2 = j0_a2
        self.neutron_magnetic_j0_b1 = j0_B1
        self.neutron_magnetic_j0_b2 = j0_b2
        self.neutron_magnetic_j0_c1 = j0_C1
        self.neutron_magnetic_j0_c2 = j0_c2
        self.neutron_magnetic_j0_d = j0_D
        self.neutron_magnetic_j2_a1 = j2_A1
        self.neutron_magnetic_j2_a2 = j2_a2
        self.neutron_magnetic_j2_b1 = j2_B1
        self.neutron_magnetic_j2_b2 = j2_b2
        self.neutron_magnetic_j2_c1 = j2_C1
        self.neutron_magnetic_j2_c2 = j2_c2
        self.neutron_magnetic_j2_d = j2_D

        self.j0_parameters = numpy.array([j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D], dtype=float)
        self.j2_parameters = numpy.array([j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D], dtype=float)

    @classmethod
    def form_by_symbol(cls, symbol: str):
        """Form by symbol (classmethod)."""
        j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, j2_B1,\
            j2_b2, j2_C1, j2_c2, j2_D = get_j0_j2_by_symbol(symbol)
        item = cls(
            symbol=symbol,
            neutron_magnetic_j0_a1=j0_A1, neutron_magnetic_j0_a2=j0_a2,
            neutron_magnetic_j0_b1=j0_B1, neutron_magnetic_j0_b2=j0_b2,
            neutron_magnetic_j0_c1=j0_C1, neutron_magnetic_j0_c2=j0_c2,
            neutron_magnetic_j0_d=j0_D,
            neutron_magnetic_j2_a1=j2_A1, neutron_magnetic_j2_a2=j2_a2,
            neutron_magnetic_j2_b1=j2_B1, neutron_magnetic_j2_b2=j2_b2,
            neutron_magnetic_j2_c1=j2_C1, neutron_magnetic_j2_c2=j2_c2,
            neutron_magnetic_j2_d=j2_D)
        item.j0_parameters = numpy.array([j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D], dtype=float)
        item.j2_parameters = numpy.array([j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D], dtype=float)
        return item


class AtomTypeScatL(LoopN):
    """Description of AtomTypeScat in loop."""

    ITEM_CLASS = AtomTypeScat
    ATTR_INDEX = "symbol"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomTypeScatL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def calc_form_factor(self, sthovl, l_lande: list, l_kappa: list,
                         flag_only_orbital: bool = False):
        """Calculate form-factor."""
        form_factor = [item.calc_form_factor(
            sthovl, float(lande), float(kappa),
            flag_only_orbital=flag_only_orbital) for item, lande, kappa in
            zip(self.items, l_lande, l_kappa)]
        return numpy.array(form_factor, dtype=float)

    @classmethod
    def form_by_symbols(cls, symbols: list):
        """Form by symbols (classmethod)."""
        items = []
        item_class = cls.ITEM_CLASS
        for symbol in symbols:
            items.append(item_class.form_by_symbol(symbol))
        obj = cls()
        obj.items = items
        return obj

# obj = AtomTypeScatL.form_by_symbols(["O2-", "Fe3+"])
# print(obj, end="\n\n")
