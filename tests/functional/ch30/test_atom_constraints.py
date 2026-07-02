import os
import cryspy

def test_atom_constraints_ch30():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "ch30.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    flag = all([hh1 == hh2 for hh1, hh2 in zip(rhochi.crystal_ch30.atom_site.wyckoff_symbol, ['e', 'e', 'c', 'c', 'c', 'e', 'f', 'f', 'f', 'f'])])
    assert flag


test_atom_constraints_ch30()