import numpy
from cryspy import file_to_globaln

from cryspy.A_functions_base.function_1_matrices import calc_chi_sq

f_name = "main.rcif"

mem_obj = file_to_globaln(f_name)

mem_obj.calc_fr()

density_point = mem_obj.density_point
diffrn = mem_obj.diffrn_exp_mnbite
diffrn_refln = diffrn.diffrn_refln

# Analytical Derivatives
fr_e = numpy.array(diffrn_refln.fr, dtype=float)
fr_s = numpy.array(diffrn_refln.fr_sigma, dtype=float)

crystal = mem_obj.crystals()[0]
l_diffrn = mem_obj.experiments()

mem_parameters = mem_obj.mem_parameters

density_point.volume_unit_cell = crystal.cell.volume
density_point.number_unit_cell = \
    mem_parameters.points_a * mem_parameters.points_b * \
    mem_parameters.points_c


chi_iso_ferro = mem_parameters.chi_ferro
chi_iso_antiferro = mem_parameters.chi_antiferro
points_a = mem_parameters.points_a
points_b = mem_parameters.points_b
points_c = mem_parameters.points_c
prior_density = mem_parameters.prior_density
flag_two_channel = mem_parameters.method == "2channel"
gof_desired = mem_parameters.gof_desired

cell = crystal.cell
space_group = crystal.space_group
atom_site = crystal.atom_site
space_group_symop = space_group.full_space_group_symop
atom_site_susceptibility = crystal.atom_site_susceptibility
l_magnetic_labes = atom_site_susceptibility.label


l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s = [], [], [], []
total_peaks = 0
for diffrn in l_diffrn:
    diffrn_orient_matrix = diffrn.diffrn_orient_matrix
    e_up = diffrn_orient_matrix.calc_e_up()
    setup = diffrn.setup
    field = float(setup.field)
    h_loc = (field*e_up[0], field*e_up[1], field*e_up[2])
    diffrn_refln = diffrn.diffrn_refln
    ind_h = numpy.array(diffrn_refln.index_h, dtype=int)
    ind_k = numpy.array(diffrn_refln.index_k, dtype=int)
    ind_l = numpy.array(diffrn_refln.index_l, dtype=int)
    total_peaks += ind_h.size
    hkl = (ind_h, ind_k, ind_l)
    fr_e = numpy.array(diffrn_refln.fr, dtype=float)
    fr_s = numpy.array(diffrn_refln.fr_sigma, dtype=float)
    v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro = \
        density_point.calc_factor_in_front_of_density_for_fm_perp(
            hkl, space_group_symop, cell, atom_site_susceptibility, h_loc,
            chi_iso_ferro=chi_iso_ferro,
            chi_iso_antiferro=chi_iso_antiferro,
            flag_two_channel=flag_two_channel)
    f_nucl = crystal.calc_f_nucl(*hkl)
    l_f_nucl.append(f_nucl)
    l_v_2d_i.append((v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro))
    l_fr_e.append(fr_e)
    l_fr_s.append(fr_s)

def temp_func(numpy_den=None):
    l_chi_sq, l_der_chi_sq = [], []
    l_der_chi_sq_f, l_der_chi_sq_a = [], []
    for diffrn, f_nucl, v_2d_i, fr_e, fr_s in \
            zip(l_diffrn, l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s):

        f_m_perp, delta_f_m_perp, delta_f_m_perp_f, delta_f_m_perp_a = \
            density_point.calc_fm(*v_2d_i)
        fr_m, delta_fr_m = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                          delta_f_nucl=None,
                                          delta_f_m_perp=delta_f_m_perp)
        delta_fr_m_f = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                      delta_f_nucl=None,
                                      delta_f_m_perp=delta_f_m_perp_f)[1]
        delta_fr_m_a = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                      delta_f_nucl=None,
                                      delta_f_m_perp=delta_f_m_perp_a)[1]

        diffrn.diffrn_refln.numpy_fr_calc = fr_m
        chi_sq, der_chi_sq = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m)
        der_chi_sq_f = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_f)[1]
        der_chi_sq_a = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_a)[1]
        l_chi_sq.append(chi_sq)
        l_der_chi_sq.append(der_chi_sq)
        l_der_chi_sq_f.append(der_chi_sq_f)
        l_der_chi_sq_a.append(der_chi_sq_a)
    # print(" ".join([f" {val:10.2f}" for val in l_chi_sq]))
    return sum(l_chi_sq), sum(l_der_chi_sq), sum(l_der_chi_sq_f), \
        sum(l_der_chi_sq_a)

chi_sq, delta_chi_sq, delta_chi_sq_f, delta_chi_sq_a = temp_func()


print("lengths: ", len(delta_chi_sq), len(density_point.items))
# fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
# 
# delta_fr_m = 
# 
# chi_sq_0, der_chi_sq = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m)


# Numerical Derivatives
delta_den = 1e-5
atom_label = "Mn1"
print("indexes_xyz numerical_deriv analytical_deriv")
for item, der_anal, der_anal_f, der_anal_a in zip(density_point.items, delta_chi_sq, delta_chi_sq_f, delta_chi_sq_a):
    if item.basin_atom_label == atom_label:
        # Minus delta
        den_orig = float(item.density)
        item.density = den_orig - delta_den
        den_minus = item.density
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_minus = (numpy.square((fr_e-fr_m)/fr_s)).sum()

        # Plus delta
        item.density = den_orig + delta_den
        den_plus = item.density
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_plus = (numpy.square((fr_e-fr_m)/fr_s)).sum()
        delta_chi_sq = (chi_sq_plus-chi_sq_minus)/(den_plus-den_minus)
        print(f"{item.index_x:3}{item.index_y:3}{item.index_z:3} {delta_chi_sq:12.3f} {der_anal:12.3f}")
        if abs(der_anal-delta_chi_sq) > 0.002:
            print("ALARM  "+50*"*")
        item.density = den_orig
       
        # Minus delta
        den_orig = float(item.density_ferro)
        item.density_ferro = den_orig - delta_den
        den_minus = item.density_ferro
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_minus = (numpy.square((fr_e-fr_m)/fr_s)).sum()

        # Plus delta
        item.density_ferro = den_orig + delta_den
        den_plus = item.density_ferro
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_plus = (numpy.square((fr_e-fr_m)/fr_s)).sum()
        delta_chi_sq = (chi_sq_plus-chi_sq_minus)/(den_plus-den_minus)
        print(f"          {delta_chi_sq:12.3f} {der_anal_f:12.3f} - ferro")
        if abs(der_anal_f-delta_chi_sq) > 0.002:
            print("ALARM  "+50*"*")
        
        item.density_ferro = den_orig
       
        # Minus delta
        den_orig = float(item.density_antiferro)
        item.density_antiferro = den_orig - delta_den
        den_minus = item.density_antiferro
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_minus = (numpy.square((fr_e-fr_m)/fr_s)).sum()

        # Plus delta
        item.density_antiferro = den_orig + delta_den
        den_plus = item.density_antiferro
        mem_obj.calc_fr()
        fr_m = numpy.array(diffrn_refln.fr_calc, dtype=float)
        chi_sq_plus = (numpy.square((fr_e-fr_m)/fr_s)).sum()
        delta_chi_sq = (chi_sq_plus-chi_sq_minus)/(den_plus-den_minus)
        print(f"          {delta_chi_sq:12.3f} {der_anal_a:12.3f} - antiferro")
        if abs(der_anal_a-delta_chi_sq) > 0.002:
            print("ALARM  "+50*"*")
        
        item.density_antiferro = den_orig
       
        