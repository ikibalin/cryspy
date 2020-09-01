"""
Module realized some mathematical functions 
"""

from typing import Tuple, Union
import numpy

def calc_m_sigma(np_M, np_sigma):
    """
Estimate errors for equation:

    np_M * np_val * np_M^T
    """
    np_sig_1 = numpy.zeros((3, 3), dtype=float)
    for _i in range(3):
        for _j in range(3):
            res = []
            for _k in range(3):
                for _l in range(3):
                    res.append((np_sigma[_k, _l]*np_M[_i, _k]*np_M[_j, _l])**2)
            np_sig_1[_i, _j] = sum(res)
    return numpy.sqrt(np_sig_1)



def calc_scalar_product_by_vectors(vector_1, vector_2, flag_vector_1=False,
                                   flag_vector_2=False):
    """
    The vector is given as tuple of its coordinates defined in Chartezian coordinate system
    """
    product = tuple([v_1*v_2 for v_1, v_2 in zip(vector_1, vector_2)])
    flag_derivatives = any([flag_vector_1, flag_vector_2])
    if flag_derivatives:
        der_vector_1, der_vector_2 = [], []
        for v_1, v_2 in zip(vector_1, vector_2):
            der_vector_1.append((v_1*0.+1.)*v_2)
            der_vector_2.append(v_1*(v_2*0.+1.))
        dder = {"der__vector_1": tuple(der_vector_1), "der__vector_2": tuple(der_vector_2)}
    else:
        dder = {}
    return product, dder


def calc_scalar_product_by_complex_vectors(vector_1, vector_2, 
                                           flag_vector_1=False, flag_vector_2=False):
    """
    The vector is given as tuple of its coordinates defined in Chartezian coordinate system
    """
    product_comp = tuple([v_1*v_2 for v_1, v_2 in zip(vector_1, vector_2)])
    product = product_comp[0] + product_comp[1] + product_comp[2]
    flag_derivatives = any([flag_vector_1, flag_vector_2])
    dder = {}
    if flag_derivatives:
        der_vector_1_re, der_vector_1_im = [], []
        der_vector_2_re, der_vector_2_im = [], []
        for v_1, v_2 in zip(vector_1, vector_2):
            der_vector_1_re.append((v_1*0.+1.)*v_2)
            der_vector_1_im.append((v_1*0.+1.)*v_2)
            der_vector_2_re.append(v_1*(v_2*0.+1.))
            der_vector_2_im.append(v_1*(v_2*0.+1.))
        if flag_vector_1:
            dder["der__vector_1_re"] = tuple(der_vector_1_re)
            dder["der__vector_1_im"] = tuple(der_vector_1_im) 
        if flag_vector_2:
            dder["der__vector_2_re"] = tuple(der_vector_2_re)
            dder["der__vector_2_im"] = tuple(der_vector_2_im) 
    return product, dder


def calc_modulus_sq_by_complex_vector(vector_1, flag_vector_1=False):
    """
    The vector is given as tuple of its coordinates defined in Chartezian coordinate system
    """
    modulus_comp = tuple([(v_1*v_1.conjugate()).real for v_1 in vector_1])
    modulus = modulus_comp[0]+modulus_comp[1]+modulus_comp[2]
    flag_derivatives = any([flag_vector_1])
    if flag_derivatives:
        der_vector_1_re, der_vector_1_im = [], []
        for v_1 in vector_1:
            der_vector_1_re.append(2.*v_1.real)
            der_vector_1_im.append(2.*v_1.imag)
        dder = {"der__vector_1_re": tuple(der_vector_1_re), 
                "der__vector_1_im": tuple(der_vector_1_im)}
    else:
        dder = {}
    return modulus, dder


def tri_linear_interpolation(ind_xyz, den_3d):
    ind_x, ind_y, ind_z = ind_xyz[:, 0], ind_xyz[:, 1], ind_xyz[:, 2]
    (n_x, n_y, n_z) = den_3d.shape 
    ind_x0 = (numpy.floor(ind_x)).astype(int)
    ind_x1 = ind_x0+1
    ind_y0 = (numpy.floor(ind_y)).astype(int)
    ind_y1 = ind_y0+1
    ind_z0 = (numpy.floor(ind_z)).astype(int)
    ind_z1 = ind_z0+1

    xd = (ind_x-ind_x0)/(ind_x1-ind_x0)
    yd = (ind_y-ind_y0)/(ind_y1-ind_y0)
    zd = (ind_z-ind_z0)/(ind_z1-ind_z0)

    Vx0y0z0, Vx1y0z0 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y0, n_y), numpy.mod(ind_z0, n_z)], den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y0, n_y), numpy.mod(ind_z0, n_z)]
    Vx0y0z1, Vx1y0z1 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y0, n_y), numpy.mod(ind_z1, n_z)], den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y0, n_y), numpy.mod(ind_z1, n_z)]
    Vx0y1z0, Vx1y1z0 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y1, n_y), numpy.mod(ind_z0, n_z)], den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y1, n_y), numpy.mod(ind_z0, n_z)]
    Vx0y1z1, Vx1y1z1 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y1, n_y), numpy.mod(ind_z1, n_z)], den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y1, n_y), numpy.mod(ind_z1, n_z)]

    c00 = Vx0y0z0*(1.-xd) + Vx1y0z0*xd
    c01 = Vx0y0z1*(1.-xd) + Vx1y0z1*xd
    c10 = Vx0y1z0*(1.-xd) + Vx1y1z0*xd
    c11 = Vx0y1z1*(1.-xd) + Vx1y1z1*xd

    c0 = c00*(1.-yd) + c10*yd
    c1 = c01*(1.-yd) + c11*yd
    
    c = c0*(1.-zd) + c1*zd

    return c
    
    

FUNCTIONS = [
calc_scalar_product_by_vectors,
calc_scalar_product_by_complex_vectors,
calc_modulus_sq_by_complex_vector,
tri_linear_interpolation
]    