import numpy
from cryspy import file_to_globaln

from cryspy.A_functions_base.function_1_matrices import calc_chi_sq

f_name = "main.rcif"

mem_obj = file_to_globaln(f_name)

mem_parameters = mem_obj.mem_parameters
points_a = mem_parameters.points_a
points_b = mem_parameters.points_b
points_c = mem_parameters.points_c

density_point = mem_obj.density_point
diffrn = mem_obj.diffrn_exp_mnbite
diffrn_refln = diffrn.diffrn_refln

crystal = mem_obj.crystals()[0]
space_group = crystal.space_group
full_space_group_symop = space_group.full_space_group_symop

np_r_11 = full_space_group_symop.numpy_r_11
np_r_12 = full_space_group_symop.numpy_r_12
np_r_13 = full_space_group_symop.numpy_r_13
np_r_21 = full_space_group_symop.numpy_r_21
np_r_22 = full_space_group_symop.numpy_r_22
np_r_23 = full_space_group_symop.numpy_r_23
np_r_31 = full_space_group_symop.numpy_r_31
np_r_32 = full_space_group_symop.numpy_r_32
np_r_33 = full_space_group_symop.numpy_r_33
np_b_1 = full_space_group_symop.numpy_b_1
np_b_2 = full_space_group_symop.numpy_b_2
np_b_3 = full_space_group_symop.numpy_b_3

number_symmetry = len(full_space_group_symop.items)

for item in density_point.items:
    ind_x = item.index_x
    ind_y = item.index_y
    ind_z = item.index_z
    mult = item.multiplicity
    
    ind_x_n = np_r_11*ind_x+np_r_12*ind_y+np_r_13*ind_z+np_b_1*points_a
    ind_y_n = np_r_21*ind_x+np_r_22*ind_y+np_r_23*ind_z+np_b_2*points_b
    ind_z_n = np_r_31*ind_x+np_r_32*ind_y+np_r_33*ind_z+np_b_3*points_c

    f_o_x, f_o_y, f_o_z = ind_x_n[0], ind_y_n[1], ind_z_n[2]
    i_mult = 0
    for i_x, i_y, i_z in zip(ind_x_n, ind_y_n, ind_z_n):
        flag_1 = (f_o_x%points_a) == (i_x%points_a)
        flag_2 = (f_o_y%points_b) == (i_y%points_b)
        flag_3 = (f_o_z%points_c) == (i_z%points_c)
        flag = ((flag_1 & flag_2) & flag_3)
        if flag:
            i_mult += 1
    
    hh = number_symmetry%i_mult
    print(f"{ind_x: 3}{ind_y: 3}{ind_z: 3}   {mult: 3}{number_symmetry//i_mult: 3}")
    if hh != 0:
        print("ALARM   "+50*"*")
    elif mult != number_symmetry//i_mult:
        print("ALARM   "+50*"*")