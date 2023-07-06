import os
import numpy
import cryspy
import cryspy.A_functions_base.structure_factor as structure_factor
import cryspy.A_functions_base.extinction as extinction
import cryspy.A_functions_base.flip_ratio as flip_ratio
import cryspy.A_functions_base.unit_cell as unit_cell

from cryspy.H_functions_global.form_dictionary import form_dictionary_by_crystal, \
  form_dictionary_by_diffrn

f_nucl = numpy.array([
        25.32575171+1.71321270e-14j,   6.46451557+4.27787562e-15j,
         6.46451557+7.98266487e-15j, -14.39094569-1.81227764e-16j,
       -14.39094569-4.03490903e-15j,  29.95022452+0.00000000e+00j,
        14.91614642+3.22131211e-17j,  -7.17462107-1.41980682e-15j,
        -7.17462107-1.41980682e-15j,   7.17462107+4.43122378e-15j,
        -7.17462107+1.41980682e-15j,  -4.62228025+1.21406786e-15j,
        -4.62228025+6.98657921e-16j,  -4.62228025+1.21406786e-15j,
        -4.62228025+6.98657921e-16j, -15.79750462-6.99988380e-15j,
       -15.79750462-6.14086724e-15j,  15.79750462+4.13650382e-15j,
       -15.15393358-5.76301755e-15j, -15.15393358-2.93217829e-15j,
       -15.15393358-4.55663063e-15j,  -4.67315933-4.66603941e-15j,
        -4.67315933-6.94249991e-15j,   4.67315933+4.66906284e-15j,
        -8.31404207-8.05624966e-15j,  -8.31404207-2.44528612e-15j,
        20.21441843+4.21584414e-15j,  20.21441843+3.78633586e-15j,
       -20.21441843-2.75676700e-15j, -20.21441843+4.32563285e-16j,
         9.03919406+2.78969243e-15j,  -9.03919406-1.16098813e-14j,
        -9.03919406-1.26704033e-14j,   9.03919406+2.77494127e-15j,
        -3.35703226-1.11994029e-14j,  -3.35703226-6.49347834e-15j,
        -7.12374199-4.60105427e-15j,  -7.12374199-7.42070770e-16j,
         2.13603032+3.85538087e-15j,   2.13603032-2.03208742e-15j,
        -7.81819211-1.45821676e-15j,  -7.81819211-1.86624962e-15j,
       -23.1332461 -9.56680990e-16j, -23.1332461 -5.48695344e-15j,
         6.91800303+3.70447267e-15j,   6.91800303+1.46812268e-15j,
         6.91800303+3.70447267e-15j,   6.91800303+1.90299981e-15j,
         7.3152    +4.11042271e-15j,   7.3152    +3.59501277e-15j,
         7.3152    +1.00363785e-14j,   7.3152    +1.05517884e-14j,
        -1.53362782+1.30733740e-15j,   1.53362782-8.45068121e-16j,
         1.53362782+1.95881123e-15j,  -1.53362782-2.25797731e-15j,
        -9.46005795-8.07420293e-15j,  -9.46005795-3.30378553e-15j,
         7.3152    +5.80321703e-15j,   7.3152    -1.17675512e-15j,
         8.02173273+9.12070402e-15j,   8.02173273+3.46554990e-15j,
        -2.75770727-3.36779658e-15j,  -2.75770727-6.45440904e-16j,
        -2.75770727-2.76648498e-15j,  -2.75770727-4.73637592e-16j,
       -28.19373094-1.37667560e-14j, -28.19373094-9.44588638e-15j,
       -17.71295669-7.55307817e-15j, -17.71295669-3.10884060e-15j,
        13.93293164+3.04145384e-15j,  13.93293164+2.35424059e-15j,
       -13.93293164-4.92013640e-15j, -13.93293164+5.21013060e-15j,
        17.01850657+4.00265911e-15j, -17.01850657-5.25879488e-15j,
       -17.01850657-8.22932241e-15j,  17.01850657+3.83868657e-15j,
         6.53773232+1.67863034e-14j,   6.53773232+1.99646647e-14j,
        -6.53773232-9.65601459e-15j,  15.46397622+7.81458859e-15j,
        15.46397622+6.24106950e-15j,  15.46397622+8.29778541e-15j,
         9.98525868-1.64834779e-15j,   9.98525868-1.61613466e-15j,
       -25.10815601-8.42659440e-15j, -25.10815601-9.62921759e-15j,
        -5.8432822 -5.37104221e-15j,  -5.8432822 -5.69200029e-15j,
        -5.8432822 -7.94809190e-15j,  -5.8432822 -3.71626220e-15j,
         7.3152    +5.01239010e-15j,   7.3152    +6.10763621e-15j], dtype=complex)

fr = numpy.array([
       0.67769245, 1.7400608 , 1.61654164, 0.92303798, 0.93245852,
       0.73714893, 0.40084516, 0.50707637, 0.50795978, 0.48982221,
       0.50707637, 1.32039548, 0.59300708, 0.5916218 , 1.20424258,
       0.26216475, 0.26492688, 0.82986437, 0.33546764, 0.33383482,
       1.12074955, 0.17855907, 0.17960172, 0.57088846, 0.17511859,
       0.17401576, 0.91998266, 0.92101488, 0.463532  , 0.46219096,
       0.93118655, 1.31636915, 1.2325219 , 0.92126681, 2.02369069,
       2.10804872, 1.2464412 , 1.1946929 , 1.5758018 , 1.57551678,
       1.0475533 , 1.02447433, 0.69226483, 0.69191805, 1.2232208 ,
       0.28342698, 0.28216954, 1.18997359, 0.14681338, 0.14755406,
       0.08987179, 0.09044616, 0.89966622, 1.5618064 , 1.43274175,
       0.91577276, 0.88814321, 0.90020169, 0.07214861, 0.07263155,
       1.06119519, 1.0503256 , 1.27122221, 0.58716649, 0.58641253,
       1.16313051, 0.7322256 , 0.73278371, 0.56721618, 0.56665512,
       0.91421339, 0.92236371, 0.30185899, 0.30081766, 0.96430564,
       1.05840772, 1.04005396, 0.96652207, 0.81260982, 0.83355983,
       0.1146756 , 0.47731122, 0.47819234, 0.51865214, 1.11041753,
       1.07842329, 0.6923462 , 0.69290315, 1.14664124, 0.07892507,
       0.07801952, 1.08325336, 0.0901045 , 0.09113749], dtype=float)

def test_calc_flip_ratio_intensity_f_nucl():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "single_crystal.rcif")
    obj = cryspy.load_file(f_name)
    crystal = obj.crystal_ho2ti2o7
    diffrn = obj.diffrn_ho2ti2o7

    dict_crystal = form_dictionary_by_crystal(crystal)
    dict_diffrn = form_dictionary_by_diffrn(diffrn)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    full_symm_elems = dict_crystal["full_symm_elems"]
    atom_fract_xyz = dict_crystal["atom_fract_xyz"]
    atom_occupancy = dict_crystal["atom_occupancy"]
    atom_scat_length_neutron = dict_crystal["atom_scat_length_neutron"] 
    atom_b_iso = dict_crystal["atom_b_iso"]
    atom_beta = dict_crystal["atom_beta"]
    mag_atom_b_iso = dict_crystal["mag_atom_b_iso"]
    mag_atom_beta = dict_crystal["mag_atom_beta"]
    mag_atom_occupancy = dict_crystal["mag_atom_occupancy"]
    mag_atom_fract_xyz = dict_crystal["mag_atom_fract_xyz"]
    lande_factor = dict_crystal["lande_factor"]
    kappa = dict_crystal["kappa"]
    j0_parameters = dict_crystal["j0_parameters"]
    j2_parameters = dict_crystal["j2_parameters"]
    mag_atom_susceptibility = dict_crystal["mag_atom_susceptibility"]

    index_hkl = dict_diffrn["index_hkl"]
    magnetic_field = dict_diffrn["magnetic_field"]

    f_nucl_new = structure_factor.calc_f_nucl(
        index_hkl, full_symm_elems, unit_cell_parameters,
        atom_fract_xyz, atom_occupancy, atom_b_iso, atom_beta, atom_scat_length_neutron)[0]
    
    print("Nuclear structure factor")

    assert numpy.all(numpy.isclose(f_nucl_new, f_nucl))

    mag_atom_multiplicity = None
    flag_unit_cell_parameters = True
    flag_mag_atom_fract_xyz = True
    flag_mag_atom_occupancy = True
    flag_mag_atom_b_iso = True
    flag_mag_atom_beta = True
    flag_lande_factor = True
    flag_kappa = True
    flag_susceptibility = True

    sthovl = unit_cell.calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]

    sft_ccs, dder_sft_ccs = structure_factor.calc_structure_factor_tensor_ccs(index_hkl, full_symm_elems, unit_cell_parameters,
      mag_atom_fract_xyz, mag_atom_occupancy, mag_atom_b_iso, mag_atom_beta,
      lande_factor, kappa, mag_atom_susceptibility, j0_parameters, j2_parameters,
      mag_atom_multiplicity = mag_atom_multiplicity,
      flag_unit_cell_parameters = flag_unit_cell_parameters,
      flag_mag_atom_fract_xyz = flag_mag_atom_fract_xyz,
      flag_mag_atom_occupancy = flag_mag_atom_occupancy,
      flag_mag_atom_b_iso = flag_mag_atom_b_iso,
      flag_mag_atom_beta = flag_mag_atom_beta,
      flag_lande_factor = flag_lande_factor,
      flag_kappa = flag_kappa,
      flag_susceptibility = flag_susceptibility)

    eq_ccs, dder_eq_ccs = unit_cell.calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)
    flag_sft_ccs = True
    flag_magnetic_field = True
    flag_eq_ccs = True

    f_m_perp, dder_f_m_perp = structure_factor.calc_f_m_perp_by_sft(
      sft_ccs, magnetic_field, eq_ccs,
      flag_sft_ccs=flag_sft_ccs,
      flag_magnetic_field=flag_magnetic_field,
      flag_eq_ccs=flag_eq_ccs)

    def func_extinction(f_sq, flag_f_sq: bool = False):
      model = dict_diffrn["extinction_model"]
      radius = dict_diffrn["extinction_radius"]
      mosaicity = dict_diffrn["extinction_mosaicity"]
      volume_unit_cell = unit_cell.calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters)[0]
      wavelength = dict_diffrn["wavelength"]
      cos_2theta = 1.-2*numpy.square(sthovl * wavelength)
      flag_radius = True
      flag_mosaicity = True
      flag_volume_unit_cell = True
      flag_cos_2theta = True
      flag_wavelength = True
      return extinction.calc_extinction_sphere(
          f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
          model, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
          flag_mosaicity=flag_mosaicity,
          flag_volume_unit_cell=flag_volume_unit_cell,
          flag_cos_2theta=flag_cos_2theta,
          flag_wavelength=flag_wavelength)
    
    polarization = diffrn.diffrn_radiation.polarization
    flipper = diffrn.diffrn_radiation.efficiency
    m_u = diffrn.diffrn_orient_matrix.u

    flag_polarization = True
    flag_flipper = True
    flag_f_n = True
    flag_f_m_perp = True

    c_lambda2 = None
    f_n_2h = None
    f_m_perp_2h = None
    flag_c_lambda2 = False
    flag_f_n_2h = False
    flag_f_m_perp_2h = False



    fr_new, dder_new = flip_ratio.calc_flip_ratio_by_structure_factors(
      polarization, flipper, f_nucl_new, f_m_perp, m_u, func_extinction=func_extinction,
      c_lambda2=c_lambda2, f_n_2h=f_n_2h, f_m_perp_2h=f_m_perp_2h,
      flag_polarization=flag_polarization, flag_flipper=flag_flipper,
      flag_f_n=flag_f_n, flag_f_m_perp=flag_f_m_perp,
      flag_c_lambda2=flag_c_lambda2,
      flag_f_n_2h=flag_f_n_2h, flag_f_m_perp_2h=flag_f_m_perp_2h)

    print("flip ratios")
    assert numpy.all(numpy.isclose(fr_new, fr))

    def numerical_derivative(function, param, delta):
      f_p, dder = function(param)
      f_p, dder_p = function(param+delta)
      f_m, dder_m = function(param-delta)
      dder_num = (f_p-f_m)/(2.*delta)
      return dder_num, dder

    def func_temp_fr(param):
      np_f_m_perp_new = numpy.copy(f_m_perp)
      ind_vec = 1 # x, y, z
      ind_hkl = 0 # list of hkl
      np_f_m_perp_new[ind_vec, ind_hkl] = numpy.real(param)
      fr_hkl, dder_fr = flip_ratio.calc_flip_ratio_by_structure_factors(
          polarization, flipper, f_nucl_new, np_f_m_perp_new, m_u, func_extinction=func_extinction,
          c_lambda2=c_lambda2, f_n_2h=f_n_2h, f_m_perp_2h=f_m_perp_2h,
          flag_polarization=flag_polarization, flag_flipper=flag_flipper,
          flag_f_n=flag_f_n, flag_f_m_perp=flag_f_m_perp,
          flag_c_lambda2=flag_c_lambda2,
          flag_f_n_2h=flag_f_n_2h, flag_f_m_perp_2h=flag_f_m_perp_2h)
      res = fr_hkl[ind_hkl]
      dder = dder_fr["f_m_perp_real"][ind_vec, ind_hkl]
      return res, dder

    a,b = numerical_derivative(func_temp_fr, 1, 0.0001)
    print("derivatives of flip ratio")
    assert numpy.isclose(a, b, atol=1e-5)

    def func_temp_iint(param):
      np_f_m_perp_new = numpy.copy(f_m_perp)
      ind_vec = 1 # x, y, z
      ind_hkl = 0 # list of hkl
      np_f_m_perp_new[ind_vec, ind_hkl] = numpy.real(param)
      axis_z = m_u[6:9]
      i_p_hkl, i_m_hkl, dder_p, dder_m = flip_ratio.calc_iint(
          polarization, flipper, f_nucl_new, np_f_m_perp_new, axis_z, func_extinction=func_extinction,
          flag_polarization=flag_polarization, flag_flipper=flag_flipper,
          flag_f_n=flag_f_n, flag_f_m_perp=flag_f_m_perp)
      res = i_p_hkl[ind_hkl]
      dder = dder_p["f_m_perp_real"][ind_vec, ind_hkl]
      return res, dder
    a,b = numerical_derivative(func_temp_iint, 1, 0.0001)
    print("derivatives of integrated intensity up")
    assert numpy.isclose(a, b, atol=1e-5)
