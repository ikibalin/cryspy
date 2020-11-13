"""Module realized some mathematical functions.

Functions
---------
    - calc_m_sigma
    - calc_scalar_product_by_vectors
    - calc_scalar_product_by_complex_vectors
    - calc_modulus_sq_by_complex_vector

"""
import numpy


def calc_m_sigma(np_M, np_sigma):
    """Estimate errors for equation.

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
    """Scalar product by vectors.

    The vector is given as tuple of its coordinates defined in Chartezian
    coordinate system
    """
    product = tuple([v_1*v_2 for v_1, v_2 in zip(vector_1, vector_2)])
    flag_derivatives = any([flag_vector_1, flag_vector_2])
    if flag_derivatives:
        der_vector_1, der_vector_2 = [], []
        for v_1, v_2 in zip(vector_1, vector_2):
            der_vector_1.append((v_1*0.+1.)*v_2)
            der_vector_2.append(v_1*(v_2*0.+1.))
        dder = {"der__vector_1": tuple(der_vector_1),
                "der__vector_2": tuple(der_vector_2)}
    else:
        dder = {}
    return product, dder


def calc_scalar_product_by_complex_vectors(
        vector_1, vector_2, flag_vector_1=False, flag_vector_2=False):
    """Scalar product by complex vectors.

    The vector is given as tuple of its coordinates defined in Chartezian
    coordinate system.
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
    """Square of modulus of complex vector.

    The vector is given as tuple of its coordinates defined in Chartezian
    coordinate system.
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
