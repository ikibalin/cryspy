import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_atom_type import AtomType, AtomTypeL

STR_FROM_CIF_1 = """
 loop_
 _atom_type_symbol
 _atom_type_oxidation_number
 _atom_type_number_in_cell
 _atom_type_scat_dispersion_real
 _atom_type_scat_dispersion_imag
 _atom_type_scat_source
   C 0 72  .017  .009  International_Tables_Vol_IV_Table_2.2B
   H 0 100  0     0    International_Tables_Vol_IV_Table_2.2B
   O 0 12  .047  .032  International_Tables_Vol_IV_Table_2.2B
   N 0 4   .029  .018  International_Tables_Vol_IV_Table_2.2B
    """

def test_init():
    try:
        _object = AtomTypeL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomTypeL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomType()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomType()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = AtomTypeL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.symbol == ["C", "H", "O", "N"]
    assert _obj["C"].number_in_cell == 72
    val = _obj.to_cif(separator=".")

