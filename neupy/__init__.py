__version__ = '1.0.0'
__author__ = 'Iurii Kibalin'
__contact__ = 'iurii.kibalin@cea.fr'

from neupy.f_common.cl_variable import Variable

from neupy.f_rcif import RCif

from neupy.f_crystal import (
    SpaceGroup,
    Cell,
    Extinction,
    AtomType,
    Fract,
    ADP,
    Magnetism,
    AtomSite,
    Crystal,
    calc_mRmCmRT
    )

from neupy.f_experiment import (
    BeamPolarization,
    AsymmetryPowder1D,
    BackgroundPowder1D,
    CalculatedDataPowder1D,
    FactorLorentzPowder1D,
    ObservedDataPowder1D,
    ResolutionPowder1D,
    SetupPowder1D,
    ExperimentPowder1D,
    AsymmetryPowder2D,
    BackgroundPowder2D,
    CalculatedDataPowder2D,
    FactorLorentzPowder2D,
    ObservedDataPowder2D,
    ResolutionPowder2D,
    SetupPowder2D,
    ExperimentPowder2D,
    ExperimentPowderTexture2D,
    CalculatedDataSingle,
    ExperimentSingle,
    ObservedDataSingle,
    SetupSingle,
    ExperimentSingleDomain,
    ObservedDataSingleDomain,
    Model,
    MemReconstruction,
    ObservedDataMem,
    CellDensity
    )



from neupy.f_common.error_simplex import error_estimation_simplex

from neupy.f_api_rcif.api_rcif_model import (
    conv_rcif_to_model,
    conv_model_to_rcif
    )

from neupy.f_api_rcif.api_rcif_crystal import (
    conv_data_to_crystal
    )

from neupy.f_api_rcif.api_rcif_mem import (
    conv_rcif_to_mem_reconstruction
    )