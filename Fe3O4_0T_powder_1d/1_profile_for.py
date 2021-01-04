import os
import pycifstar
from cryspy_future import RhoChi

f_dir = os.path.dirname(__file__)
f_rcif = os.path.join(f_dir, "main.rcif")

STR_FROM_CIF = str(pycifstar.to_global(f_rcif))

obj = RhoChi.from_cif(STR_FROM_CIF)
d_res = obj.refine()
# print(d_res)
print(obj.crystal_fe3O4, end=2*"\n")
# for name in obj.get_variable_names():
#     print(f"{name[-1][0].ljust(12):}: {obj.get_variable_by_name(name): .5f}")

# for obj_e in obj.experiments():
#     for item in obj_e.items:
#         print(item.get_name())
