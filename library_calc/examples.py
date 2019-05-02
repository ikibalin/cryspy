# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""

import matplotlib.pyplot
import numpy

from crystal import *
from calculated_data import *
from experiment import *
from observed_data import *
from variable import *

from fitting import *

from read_rcif import *

#crystal definition
cell = Cell(a=8.5502, singony="Cubic")

fe_a = AtomType(type_n="Fe", type_m="Fe3", x=0.125, y=0.125, z=0.125,
                chi_11=-3.616, chi_22=-3.616, chi_33=-3.616)
fe_b = AtomType(type_n="Fe", type_m="Fe3", x=0.500, y=0.500, z=0.500,
                chi_11=3.1846, chi_22=3.1846, chi_33=3.1846)
o = AtomType(type_n="O", x=0.25223, y=0.25223, z=0.25223)


atom_site = AtomSite()
atom_site.add_atom(fe_a)
atom_site.add_atom(fe_b)
atom_site.add_atom(o)

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m", spgr_choice="2")

crystal = Crystal(space_groupe, cell, atom_site)



# single crystal polarized neutron diffraction
calculated_data_single = CalculatedDataSingle(crystal=crystal)
experiment_single= ExperimentSingle()
experiment_single.add_calculated_data(calculated_data_single)

observed_data_single = ObservedDataSingle()
f_inp = os.path.join("HoTi_single", "full.dat")
observed_data_single.read_data(f_inp)
experiment_single.set_val(observed_data=observed_data_single)

h = observed_data_single.get_val("h")
k = observed_data_single.get_val("k")
l = observed_data_single.get_val("l")

f_r_exp = observed_data_single.get_val("flip_ratio")
sf_r_exp = observed_data_single.get_val("sflip_ratio")

chi_sq, n = experiment_single.calc_chi_sq()


f_r_mod = experiment_single.calc_iint_u_d_flip_ratio(h, k, l)[2]
matplotlib.pyplot.scatter(f_r_mod, f_r_exp, f_r_exp)






# 2 dimensional polarized powder diffraction
calculated_data_powder_2d = CalculatedDataPowder2D(scale=0.0232, crystal=crystal)
experiment_powder_2d= ExperimentPowder2D()
experiment_powder_2d.add_calculated_data(calculated_data_powder_2d )

observed_data_powder_2d = ObservedDataPowder2D()
f_inp = os.path.join("Fe3O4_150K_6T_2d", "full.dat")
observed_data_powder_2d.read_data(f_inp)
experiment_powder_2d.set_val(observed_data=observed_data_powder_2d)

chi_sq, n = experiment_powder_2d.calc_chi_sq()

np_tth = observed_data_powder_2d.get_val('tth')
np_phi = observed_data_powder_2d.get_val('phi')
int_u_exp = observed_data_powder_2d.get_val("int_u")
sint_u_exp = observed_data_powder_2d.get_val("sint_u")
int_d_exp = observed_data_powder_2d.get_val("int_d")
sint_d_exp = observed_data_powder_2d.get_val("sint_d")


int_u, int_d = experiment_powder_2d.calc_profile(np_tth, np_phi)







# 1 dimensional polarized powder diffraction
calculated_data_powder_1d = CalculatedDataPowder1D(scale=0.0232, crystal=crystal)
experiment_powder_1d = ExperimentPowder1D()
experiment_powder_1d.add_calculated_data(calculated_data_powder_1d)


setup_powder_1d = ExperimentPowder1D().get_val("setup")
setup_powder_1d.set_val(wave_length=0.84, zero_shift=-0.412)
resolution = setup_powder_1d.get_val("resolution")
resolution.set_val(u=17.487, v=-2.8357, w=0.576309)



observed_data_powder_1d = ObservedDataPowder1D()
f_inp = os.path.join("Fe3O4_0T", "full.dat")
observed_data_powder_1d.read_data(f_inp)
experiment_powder_1d.set_val(observed_data=observed_data_powder_1d)


chi_sq, n = experiment_powder_1d.calc_chi_sq()

np_tth = observed_data_powder_1d.get_val('tth')
int_u_exp = observed_data_powder_1d.get_val("int_u")
int_d_exp = observed_data_powder_1d.get_val("int_d")

int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp+int_d_exp, "k-", np_tth, int_u+int_d, "b-")







#data refinement
dd_1 = Variable(val=0, ref=True)
setup_powder_1d.set_val(zero_shift=dd_1)

dd_2 = Variable(val=0.023, ref=True)
calculated_data_powder_1d.set_val(scale=dd_2)
 

fitting = Fitting()
fitting.add_experiment(experiment_powder_1d)
fitting.add_variable(dd_1)
fitting.add_variable(dd_2)
res = fitting.refinement()



int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp+int_d_exp, "k-", np_tth, int_u+int_d, "b-")


int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp-int_d_exp, "k-", np_tth, int_u-int_d, "b-")




rcif = RCif()
f_inp = os.path.join("HoTi_single", "full.rcif")
rcif.load_from_file(f_inp)
data = rcif.data
len(data)
data[0]["loops"]
data[1]["loops"]
data[2]["loops"]
