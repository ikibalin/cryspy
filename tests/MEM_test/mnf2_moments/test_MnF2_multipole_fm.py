import os 
import cryspy


from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

def test_fm_by_multipole_for_mnf2():
    f_name = os.path.join(os.path.dirname(__file__), "MnF2_multipole_prior.rcif")
    rcif_obj = cryspy.load_file(f_name)
    d_rcif = rcif_obj.get_dictionary()
    dict_crystal = d_rcif["crystal_mnf2"]
    dict_diffrn = d_rcif["diffrn_exp1"]
    dict_in_out = {}
    chi_sq, n_hkl, der_chi_sq, dder_chi_sq, l_parameter_name = calc_chi_sq_for_diffrn_by_dictionary(dict_diffrn, dict_crystal, dict_in_out)
    # res =calculate_moments_by_multipole_for_mnf2(rcif_obj)
    assert chi_sq < 800.

