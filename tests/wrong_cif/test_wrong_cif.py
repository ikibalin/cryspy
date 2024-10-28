import os
import cryspy

def test_wrong_cif():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "wrong_cif.rcif")
    globaln = cryspy.file_to_globaln(f_name)
    assert globaln.data_Tb2Ti2O7.geom_angle[0].label is None
    assert globaln.data_tb2ti2o7.atom_type["TI"].scat_length_neutron.real == -0.3438
    assert globaln.data_global.journal.name_full == "Journal Name"
