import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_proc import PdProc, PdProcL

STR_FROM_CIF_1 = """
 loop_
 _pd_proc_ttheta
 _pd_proc_ttheta_corrected
 _pd_proc_d_spacing
 _pd_proc_intensity_up_net
 _pd_proc_intensity_down_net
 _pd_proc_intensity_up_total
 _pd_proc_intensity_down_total
 _pd_proc_intensity_bkg_calc
 _pd_proc_intensity_up
 _pd_proc_intensity_up_sigma
 _pd_proc_intensity_down
 _pd_proc_intensity_down_sigma
  4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000
"""

def test_init():
    try:
        _object = PdProcL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdProcL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = PdProc()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = PdProc()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = PdProcL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.ttheta == [4.00]
    assert _obj[4.00].intensity_up == 465.80000
    val = _obj.to_cif(separator=".")

