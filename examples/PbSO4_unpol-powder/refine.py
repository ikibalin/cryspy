import os, sys

from pycifstar import to_global
from cryspy.scripts.cl_rhochi import RhoChi

current_dir = os.path.dirname(sys.argv[0])
path_to_project_dir = os.path.join(current_dir)
main_rcif_path = os.path.join(path_to_project_dir, "main.rcif")

rhochi = RhoChi.from_cif(str(to_global(main_rcif_path)))
rhochi.apply_constraint()
#print(rhochi.experiments)

res = rhochi.refine()
print(res)
