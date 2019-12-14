import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.magneticcif.cl_atom_site_scat import AtomSiteScat, AtomSiteScatL

STR_FROM_CIF_1 = """
 loop_
 _atom_site_scat_label
 _atom_site_scat_lande
 _atom_site_scat_kappa
  Fe3A    2.00 1.00
  Fe3B    2.00 1.00
    """

def test_init():
    try:
        _object = AtomSiteScatL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomSiteScatL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomSiteScat()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomSiteScat()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.0001
    l_obj = AtomSiteScatL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.label == ["Fe3A", "Fe3B"]
    assert math.isclose(_obj["Fe3A"].lande.value, 2.000, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(_obj["Fe3B"].kappa.value, 1.000, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

