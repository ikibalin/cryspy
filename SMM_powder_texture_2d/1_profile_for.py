import os
from cryspy_future import RhoChi

f_dir = os.path.dirname(__file__)
f_rcif = os.path.join(f_dir, "main.rcif")

obj = RhoChi.from_cif_file(f_rcif)
d_res = obj.refine()
print(d_res)