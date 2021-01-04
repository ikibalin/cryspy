import os
import pycifstar
from cryspy_future import RhoChi

f_dir = os.path.dirname(__file__)
f_rcif = os.path.join(f_dir, "main.rcif")

STR_FROM_CIF = str(pycifstar.to_global(f_rcif))

obj = RhoChi.from_cif(STR_FROM_CIF)
d_res = obj.refine()
print(d_res)
with open("test.rcif", "w") as fid:
    fid.write(str(obj))
