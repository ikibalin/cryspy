
import os
import math
import numpy


from cryspy.rhocif.cl_atom_rho_orbital_radial_slater import AtomRhoOrbitalRadialSlater, AtomRhoOrbitalRadialSlaterL


f_data = os.path.join(os.path.dirname(__file__), "data.rcif")

with open(f_data, "r") as fid:
    STR_FROM_CIF_1 = fid.read()

objects = AtomRhoOrbitalRadialSlaterL.from_cif(STR_FROM_CIF_1)

print(len(objects))
