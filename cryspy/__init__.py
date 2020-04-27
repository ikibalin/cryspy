"""
To run rhochi refinement from a command-line::

    python -m cryspy file_name
    
where file_name is a rcif file name

To create a template in the folder type in a command-line::

    python -m cryspy
"""
name = "cryspy"
__version__ = '0.2.0'

from .common.cl_global_constr import GlobalConstr
from .common.cl_data_constr import DataConstr
from .common.cl_loop_constr import LoopConstr
from .common.cl_item_constr import ItemConstr
from .common.cl_fitable import Fitable

from .cif_like.cl_crystal import Crystal
from .cif_like.cl_pd import Pd
from .cif_like.cl_pd2d import Pd2d
from .cif_like.cl_diffrn import Diffrn

from .pd2dcif_like.cl_pd2d_meas import Pd2dMeas
from .pd2dcif_like.cl_pd2d_proc import Pd2dProc

from .cif_like.cl_chi2 import Chi2 
from .cif_like.cl_diffrn_radiation import DiffrnRadiation 
from .cif_like.cl_extinction import Extinction 
from .cif_like.cl_range import Range 
from .cif_like.cl_setup import Setup 
from .cif_like.cl_texture import Texture 
from .cif_like.cl_diffrn_refln import DiffrnRefln, DiffrnReflnL 
from .cif_like.cl_exclude import Exclude, ExcludeL 
from .cif_like.cl_phase import Phase, PhaseL 

from .corecif.cl_atom_site import AtomSite, AtomSiteL
from .corecif.cl_atom_site_aniso import AtomSiteAniso, AtomSiteAnisoL
from .corecif.cl_atom_type import AtomType, AtomTypeL
from .corecif.cl_cell import Cell
from .corecif.cl_diffrn_orient_matrix import DiffrnOrientMatrix
from .corecif.cl_refine_ls import RefineLs
from .corecif.cl_refln import Refln, ReflnL

from .magneticcif.cl_atom_site_moment import AtomSiteMoment, AtomSiteMomentL
from .magneticcif.cl_atom_site_scat import AtomSiteScat, AtomSiteScatL
from .magneticcif.cl_atom_site_susceptibility import AtomSiteSusceptibility, AtomSiteSusceptibilityL
from .magneticcif.cl_atom_type_scat import AtomTypeScat, AtomTypeScatL
from .magneticcif.cl_refln_susceptibility import ReflnSusceptibility, ReflnSusceptibilityL

from .symcif.cl_space_group import SpaceGroup
from .symcif.cl_space_group_symop import SpaceGroupSymop, SpaceGroupSymopL
from .symcif.cl_space_group_wyckoff import SpaceGroupWyckoff, SpaceGroupWyckoffL

from .pd1dcif_like.cl_pd_background import PdBackground, PdBackgroundL
from .pd1dcif_like.cl_pd_instr_reflex_asymmetry import PdInstrReflexAsymmetry
from .pd1dcif_like.cl_pd_instr_resolution import PdInstrResolution
from .pd1dcif_like.cl_pd_meas import PdMeas, PdMeasL
from .pd1dcif_like.cl_pd_peak import PdPeak, PdPeakL
from .pd1dcif_like.cl_pd_proc import PdProc, PdProcL

from .rhocif.cl_atom_local_axes import AtomLocalAxes, AtomLocalAxesL
from .rhocif.cl_atom_rho_orbital_radial_slater import AtomRhoOrbitalRadialSlater, AtomRhoOrbitalRadialSlaterL
from .rhocif.cl_atom_electron_configuration import AtomElectronConfiguration, AtomElectronConfigurationL

from .scripts.cl_rhochi import RhoChi
