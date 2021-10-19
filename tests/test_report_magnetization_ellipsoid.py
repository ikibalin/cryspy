import os
import sys
if not(os.getcwd()  in sys.path):
    sys.path.append(os.getcwd())

from cryspy.D_functions_item_loop.function_1_report_magnetization_ellipsoid import \
        magnetization_ellipsoid_by_u_ij, AtomSiteSusceptibility, Cell


def test_uij():


    cell = Cell(
        length_a = 9.90710, length_b = 12.70190, length_c = 14.15560,
        angle_alpha = 90., angle_beta = 92.82, angle_gamma = 90.)
    
    a_s_s = AtomSiteSusceptibility(
        chi_11=0.350, chi_22=0.454, chi_33=2.620,
        chi_12=-0.299, chi_13=0.233, chi_23=1.066)
    
    u_ij = magnetization_ellipsoid_by_u_ij(cell, a_s_s)
    print(u_ij)
    assert True
