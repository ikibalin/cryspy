import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.magneticcif.cl_atom_type_scat import AtomTypeScat, AtomTypeScatL

STR_FROM_CIF_1 = """
loop_
_atom_type_scat_symbol
_atom_type_scat_neutron_magnetic_j0_A1
_atom_type_scat_neutron_magnetic_j0_a2
_atom_type_scat_neutron_magnetic_j0_B1
_atom_type_scat_neutron_magnetic_j0_b2
_atom_type_scat_neutron_magnetic_j0_C1
_atom_type_scat_neutron_magnetic_j0_c2
_atom_type_scat_neutron_magnetic_j0_D
 O2   0.99895  12.09652  0.28854  0.12914  0.11425 -0.22968 -0.40685 
 N2   1.00581  13.37218 -0.05868  0.07792 -0.00444  0.01678  0.05146

    """

def test_init():
    try:
        _object = AtomTypeScatL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomTypeScatL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomTypeScat()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomTypeScat()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.0001
    l_obj = AtomTypeScatL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.symbol == ["O2", "N2"]
    assert math.isclose(_obj["O2"].neutron_magnetic_j0_b1, 0.28854, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

