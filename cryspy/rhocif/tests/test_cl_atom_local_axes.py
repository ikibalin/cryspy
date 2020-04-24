import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.rhocif.cl_atom_local_axes import AtomLocalAxes, AtomLocalAxesL

STR_FROM_CIF_1 = """
loop_
_atom_local_axes_atom_label
_atom_local_axes_atom0
_atom_local_axes_ax1
_atom_local_axes_atom1
_atom_local_axes_atom2
_atom_local_axes_ax2
    Ni2+(1)  DUM0      +Z    Ni2+(1)  N(1)      -x
    """

def test_init():
    try:
        _object = AtomLocalAxesL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomLocalAxesL()
        _str = _object.to_cif()
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomLocalAxes()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomLocalAxes()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = AtomLocalAxesL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.atom0 == ["DUM0"]
    assert _obj.ax2 == ["-X"]
    assert  _obj["Ni2+(1)"].ax1 == "Z"
    val = _obj.to_cif(separator=".")

