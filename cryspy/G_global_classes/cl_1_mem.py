"""Description of MEM class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from typing import NoReturn

from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn

from cryspy.F_functions_data.script_1_mem import maximize_entropy


class MEM(GlobalN):
    """
    Class to describe RhoChi class.

    Attributes
    ----------
        - crystal_#name (mandatory)
        - diffrn_#name
        - density_point_#name

    Methods
    -------
        - crystals()
        - experiments()
        - maximize_entropy()
        - apply_constraints()
    """

    CLASSES_MANDATORY = (Crystal, Diffrn)
    CLASSES_OPTIONAL = (DensityPointL, )
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "mem"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, global_name=None, **kwargs) -> NoReturn:
        super(MEM, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["global_name"] = global_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def form_object(self):
        """Form object."""
        self.apply_constraint()

    def apply_constraint(self):
        """Apply constraints."""
        for item in self.items:
            if isinstance(item, Crystal):
                item.apply_constraint()

    def experiments(self):
        """List of expreiments."""
        return [item for item in self.items if isinstance(item, Diffrn)]

    def crystals(self):
        """List of crystals."""
        return [item for item in self.items if isinstance(item, Crystal)]

    def maximize_entropy(self, chi_iso_ferro: float = 0., 
                         chi_iso_antiferro: float = 0.,
                         n_x: int = 48, n_y: int = 48, n_z: int = 48,
                         prior_density: str = "uniform", disp: bool = True):
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        c_lambda = 1e-7
        n_cycle = 30000

        density_point = maximize_entropy(
            crystal, l_diffrn, c_lambda=c_lambda, n_cycle=n_cycle,
            chi_iso_ferro=chi_iso_ferro, chi_iso_antiferro=chi_iso_antiferro,
            n_x=n_x, n_y=n_y, n_z=n_z, prior_density=prior_density, disp=disp)

        self.add_items([density_point])

    def refine_susceptibility(self, chi_iso_ferro: float = 0., 
                         chi_iso_antiferro: float = 0.,
                         n_x: int = 48, n_y: int = 48, n_z: int = 48,
                         prior_density: str = "uniform", disp: bool = True):
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        density_point = self.density_point
