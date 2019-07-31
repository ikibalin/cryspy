# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:08 2019

@author: Iurii Kibalin
"""
import os
import matplotlib.pyplot
import numpy


import f_crystal.cl_crystal 
import f_experiment.f_single.cl_calculated_data_single 
import f_experiment.f_powder_1d.cl_calculated_data_powder_1d 
import f_experiment.f_powder_2d.cl_calculated_data_powder_2d 
import f_experiment.f_single.cl_experiment_single 
import f_experiment.f_powder_1d.cl_experiment_powder_1d 
import f_experiment.f_powder_2d.cl_experiment_powder_2d 
import f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d 
import f_experiment.f_single.cl_observed_data_single 
import f_experiment.f_powder_1d.cl_observed_data_powder_1d 
import f_experiment.f_powder_2d.cl_observed_data_powder_2d 
import f_common.cl_variable 
import f_rhochi_model.cl_model 

import f_rcif.cl_rcif 
import f_api_rcif.api_rcif_crystal
import f_api_rcif.api_rcif_mem


spgr = f_crystal.cl_crystal.SpaceGroupe("Fd-3m", "2")

n_a, n_b, n_c = 30, 30, 30
l_assymm, l_symm = spgr.calc_assymmetric_cell(n_a, n_b, n_c)

n_a, n_b, n_c = 24, 24, 24
np_x = [float(el[0])/]

x_at, y_at, z_at = 2.3, 0.2, 0.32

np_x = [float(hh[0])/float(n_a) for hh in l_assym]
np_y = [float(hh[1])/float(n_b) for hh in l_assym]
np_z = [float(hh[2])/float(n_c) for hh in l_assym]

np_diff_x_at = np_x - x_at
np_diff_y_at = np_y - y_at
np_diff_z_at = np_z - z_at


l_n_x = [hh for hh in range(10)]
for i_n_x, n_x in enumerate(l_n_x):
      print("1: ", i_n_x, "  ", n_x)
      n_x_2 = n_x+1

      if n_x_2 in l_n_x:
            l_n_x.remove(n_x_2)


len(l_assymm)

len(l_symm)
len(l_symm[0])


f_name = r"C:\docs\working_comp\RhoChi\version_next\examples\mem\full.rcif"

rcif = f_rcif.cl_rcif.RCif()

rcif.load_from_file(f_name)

mem_reconstruction = f_api_rcif.api_rcif_mem.conv_rcif_to_mem_reconstruction(rcif)

rcif_2 = f_api_rcif.api_rcif_mem.conv_mem_reconstruction_to_rcif(mem_reconstruction)

cell_density = mem_reconstruction.get_val("cell_density")



dir_name = r"C:\Users\ikibalin\Documents\documents\working_comp\RhoChi\version_next\examples"
os.chdir(dir_name)
# single crystal polarized neutron diffraction
rcif_single = f_rcif.cl_rcif.RCif()
f_inp = os.path.join("Fe3O4_150K_6T_2d", "full.rcif")
f_inp = os.path.join("Fe3O4_0T", "full.rcif")
f_inp = os.path.join("HoTi_single", "full.rcif")
f_inp = os.path.join("EuTiO3_single_domain", "full.rcif")
f_inp = os.path.join("SMM_powder_texture_2d", "full.rcif")
rcif_single.load_from_file(f_inp)

model = api_rcif_model.conv_rcif_to_model(rcif_single)
model._list_experiment

d_info_in = {}
res = model.refine_model(d_info_in)


rcif_2 = api_rcif_model.conv_model_to_rcif(model)

experiment = model._list_experiment[0]
l_crystal = model._list_crystal
experiment.save_exp_mod_data(l_crystal)
s_out = experiment.print_report(l_crystal)
print(s_out)

experiment.get_val("scale_domain")
observed_data = experiment.get_val("observed_data")







experiment_single = model_single._list_experiment[0]



observed_data_single = experiment_single.get_val("observed_data")

h = observed_data_single.get_val("h")
k = observed_data_single.get_val("k")
l = observed_data_single.get_val("l")

f_r_exp = observed_data_single.get_val("flip_ratio")
sf_r_exp = observed_data_single.get_val("sflip_ratio")

chi_sq, n_point = experiment_single.calc_chi_sq()


f_r_mod = experiment_single.calc_iint_u_d_flip_ratio(h, k, l)[2]
matplotlib.pyplot.scatter(f_r_mod, f_r_exp, f_r_exp)





# 2 dimensional polarized powder diffractionrcif = RCif()
rcif_powder_2d = RCif()
f_inp = os.path.join("Fe3O4_150K_6T_2d", "full.rcif")
rcif_powder_2d.load_from_file(f_inp)

model_powder_2d = rcif_powder_2d.trans_to_model()

d_map_powder_2d = {}
model_powder_2d.refine_model(d_map_powder_2d)

model_powder_2d.get_variables()



observed_data_powder_2d = experiment_powder_2d.get_val("observed_data")

np_tth = observed_data_powder_2d.get_val('tth')
np_phi = observed_data_powder_2d.get_val('phi')
int_u_exp = observed_data_powder_2d.get_val("int_u")
sint_u_exp = observed_data_powder_2d.get_val("sint_u")
int_d_exp = observed_data_powder_2d.get_val("int_d")
sint_d_exp = observed_data_powder_2d.get_val("sint_d")


matplotlib.pyplot.imshow(int_u_exp-int_d_exp)


int_u_mod, int_d_mod = experiment_powder_2d.calc_profile(np_tth, np_phi, d_map)

matplotlib.pyplot.imshow(int_u_mod-int_d_mod)
print(" chi_sq/n:   {:.3f}\n root of it: {:.3f}".format(chi_sq/n, 
      (chi_sq/n)**0.5))

matplotlib.pyplot.imshow(int_u_mod+int_d_mod)
matplotlib.pyplot.imshow(int_u_exp+int_d_exp)

matplotlib.pyplot.imshow(int_u_mod+int_d_mod-int_u_exp-int_d_exp)
matplotlib.pyplot.imshow(int_u_mod-int_d_mod-int_u_exp+int_d_exp)






# 1 dimensional polarized powder diffraction
rcif_powder_1d = RCif()
f_inp = os.path.join("Fe3O4_0T", "full.rcif")
rcif_powder_1d.load_from_file(f_inp)

model_powder_1d = rcif_powder_1d.trans_to_model()
d_map_powder_1d = {}
model_powder_1d.refine_model(d_map_powder_1d)

model_powder_1d.get_variables()

del model_powder_1d
del rcif_powder_1d

experiment_powder_1d = model_powder_1d._list_experiment[0]
calculated_data_powder_1d = experiment_powder_1d._list_calculated_data[0]


observed_data_powder_1d = experiment_powder_1d.get_val("observed_data")


chi_sq, n = experiment_powder_1d.calc_chi_sq()

np_tth = observed_data_powder_1d.get_val('tth')
int_u_exp = observed_data_powder_1d.get_val("int_u")
int_d_exp = observed_data_powder_1d.get_val("int_d")

int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp+int_d_exp, "k.", 
                       np_tth, int_u+int_d, "b-",
                       np_tth, int_u+int_d-int_u_exp-int_d_exp, "k-")
print(" chi_sq/n:   {:.3f}\n root of it: {:.3f}".format(chi_sq/n, 
      (chi_sq/n)**0.5))






#data refinement
dd_1 = Variable(val=0, ref=True)
setup_powder_1d = experiment_powder_1d.get_val("setup")
setup_powder_1d.set_val(zero_shift=dd_1)

dd_2 = Variable(val=0.023, ref=True)
calculated_data_powder_1d.set_val(scale=dd_2)

model.add_variable(dd_1)
model.add_variable(dd_2)
res = model.refine_model()





int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp+int_d_exp, "k-", np_tth, int_u+int_d, "b-")


int_u, int_d = experiment_powder_1d.calc_profile(np_tth)
matplotlib.pyplot.plot(np_tth, int_u_exp-int_d_exp, "k-", np_tth, int_u-int_d, "b-")







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


