from typing import NoReturn
import re

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomElectronConfiguration(ItemN):
    """
    Define the electron configuration for atoms.
    The elcetron configuration is separated on the core part which is spherical
    and non-spherical valence part. 

    The orientation and population of valence part is described in the object
    AtomRhoOrbitalValence of cryspy library.

    """
    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = ("core", "valence")
    ATTR_OPTIONAL_TYPES = (str, str)
    ATTR_OPTIONAL_CIF = ("core", "valence")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
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

    PREFIX = "atom_electron_configuration"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomElectronConfiguration, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"mosaicity": 0., "radius": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_core_shells_populations(self):
        """
        Give list of shels and their population for core part
        """
        l_res = []
        s_core = self.core
        if s_core is not None:
            for _s_shell in s_core.split():
                m = re.match(r"(\d*)([s,p,d,f])(\d*)", _s_shell)
                if m is not None:
                    res = (f"{m.group(1):}{m.group(2):}", int(m.group(3)))
                    l_res.append(res)
        return l_res

    def get_valence_shells(self):
        """
        Give list of shels for valence part
        """
        l_res = []
        s_valence = self.valence
        if s_valence is not None:
            for _s_shell in s_valence.split():
                m = re.match(r"(\d*)([s,p,d,f])", _s_shell)
                if m is not None:
                    res = f"{m.group(1):}{m.group(2):}"
                    l_res.append(res)
        return l_res

class AtomElectronConfigurationL(LoopN):
    """
    Define the electron configuration for atoms.
    The elcetron configuration is separated on the core part which is spherical
    and non-spherical valence part. 

    The orientation and population of valence part is described in the object
    AtomRhoOrbitalValence of cryspy library.
    """
    ITEM_CLASS = AtomElectronConfiguration
    ATTR_INDEX = "label"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomElectronConfigurationL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
        
    def get_core_shells_populations(self):
        """
        Give list of shels and their population for core part of all atoms 
        """
        ll_res = [item.get_core_shells_populations() for item in self.items]
        return ll_res

# s_cont = """
# loop_
# _atom_electron_configuration_label
# _atom_electron_configuration_core
# _atom_electron_configuration_valence
#  O1 "1s2" "2s 2p"
#  O2 "1s2" "2s 2p"
#  O3 "1s2" "2s 2p"
# Co1 "1s2 2s2 2p6 3s2 3p6" .
# Co2 "1s2 2s2 2p6 3s2 3p6" .
#  V1 "1s2 2s2 2p6 3s2 3p6" .
# """

# obj = AtomElectronConfigurationL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["V1"].get_core_shells_populations(), end="\n\n")
