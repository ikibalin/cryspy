import os
import pycifstar
from cryspy.scripts.cl_rhochi import RhoChi

f_dir = os.path.dirname(__file__)
f_rcif = os.path.join(f_dir, "main.rcif")

STR_FROM_CIF = str(pycifstar.to_global(f_rcif))

_obj = RhoChi.from_cif(STR_FROM_CIF)
_obj.apply_constraint()
d_res = _obj.refine()
