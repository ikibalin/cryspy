import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd2dcif_like.cl_pd2d_proc import Pd2dProc

STR_FROM_CIF_1 = """
 _pd2d_proc_ttheta_phi_intensity_up_net
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_down_net
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_up_total
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_down_total
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_bkg_calc
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_up
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_up_sigma
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_down
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_ttheta_phi_intensity_down_sigma
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
"""

def test_init():
    try:
        _object = Pd2dProc()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Pd2dProc()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = Pd2dProc()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = Pd2dProc()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    _obj = Pd2dProc.from_cif(STR_FROM_CIF_1)
    assert _obj.phi[0] == -3.00
    assert _obj.phi[1] == 41.00
    assert _obj.ttheta[1] == 40.00
    assert _obj.intensity_up_net[0][1] == -357.00
