import os
import cryspy
f_name = "/home/ikibalin/links/working_comp/cryspy/cryspy-examples/rhochi_single_susceptibility_Ho2Ti2O7_simulation_und/rhochi_single_susceptibility_Ho2Ti2O7.rcif"
rcif = cryspy.file_to_globaln(f_name)
res = cryspy.rhochi_rietveld_refinement(rcif)
