# %%
import os
import cryspy

# dir = os.path.dirname(__file__)
f_name = os.path.join(".", "main.rcif")
rhochi = cryspy.file_to_globaln(f_name)

# %%
