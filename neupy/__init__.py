__version__ = '1.0.0'
__author__ = 'Iurii Kibalin'
__contact__ = 'iurii.kibalin@cea.fr'

from neupy.f_common.cl_fitable import Fitable


from neupy.f_crystal import (
    SpaceGroup,
    Cell,
    AtomType,
    AtomSite,
    Crystal,
    calc_mRmCmRT,
    AtomSiteMagnetism,
    AtomSiteMagnetismAniso,
    AtomSiteAniso
    )

from neupy.f_experiment import (
    BeamPolarization,
    DiffrnRefln, 
    Diffrn
    )



#from neupy.f_common.error_simplex import error_estimation_simplex

#from neupy.f_api_rcif.api_rcif_model import (
#    conv_rcif_to_model,
#    conv_model_to_rcif
#    )

#from neupy.f_api_rcif.api_rcif_crystal import (
#    conv_data_to_crystal
#    )
#
#from neupy.f_api_rcif.api_rcif_mem import (
#    conv_rcif_to_mem_reconstruction
#    )