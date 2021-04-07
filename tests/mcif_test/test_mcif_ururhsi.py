import os
import numpy
import cryspy

from cryspy.D_functions_item_loop.function_1_calc_for_magcrystal import calc_f_mag

def test_f_mag():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "1.0.2_URu0.96Rh0.04Si2.mcif")
    globaln = cryspy.file_to_globaln(f_name)
    mag_cif = globaln.mag_crystal_5yOhtAoR
    
    index_hkl = numpy.array(
    [[1, 2, 0, 1, 0, 3, 2, 2, 4],
     [0, 0, 0, 1, 1, 0, 0, 1, 0],
     [1, 0, 2, 0, 1, 1, 2, 1, 0]], dtype=int)
    space_group_symop_magn_operation = mag_cif.space_group_symop_magn_operation
    space_group_symop_magn_centering = mag_cif.space_group_symop_magn_centering
    cell = mag_cif.cell
    atom_site = mag_cif.atom_site
    try:
        atom_site_aniso = mag_cif.atom_site_aniso
    except AttributeError as e:
        atom_site_aniso = None
    atom_site_scat = mag_cif.atom_site_scat
    atom_site_scat.load_atom_type_scat_by_atom_site(atom_site)
    atom_site_moment = mag_cif.atom_site_moment
    flag_derivatives = False
    flag_only_orbital= False
    f_mag, dder = calc_f_mag(index_hkl, space_group_symop_magn_operation,
    space_group_symop_magn_centering, cell, atom_site,
    atom_site_aniso, atom_site_scat, atom_site_moment,
    flag_derivatives=flag_derivatives, flag_only_orbital=flag_only_orbital)
    assert abs(f_mag[2][0].real - (-2.22872*0.2695)) < 0.001

def test_refine_mcif():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "pd_test.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    chi_sq, n_points = rhochi.calc_chi_sq()
    assert chi_sq < 3036.0
    assert int(n_points) == 4901
