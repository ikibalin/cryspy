import os
import sys
from cryspy_future import RhoChi

current_dir = os.path.dirname(sys.argv[0])
path_to_project_dir = os.path.join(current_dir)
main_rcif_path = os.path.join(path_to_project_dir, "main.rcif")

rhochi = RhoChi.from_cif_file(main_rcif_path)
res = rhochi.refine()
print(res)
