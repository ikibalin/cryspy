"""Module realized error estimation after minimization procedure by simplex.

Functions
---------
    - error_estimation_simplex

"""
import numpy


def error_estimation_simplex(vertex_vector_h, vertex_chi_sq_h, func):
    """
    Error estimation.

    Calculations according to
    ANALYTICAL CHEMISTRY, VOL. 60, NO. 8, APRIL 15, 1988

    notations
    ---------
        theta_i = vertex_vector[i, :]
        chi_sq_i = vertex_chi_sq[i]
    """
    # print("\nvertex_vector")
    # print(vertex_vector_h)
    # print("\nvertex_chi_sq")
    # print(vertex_chi_sq_h)

    # temporary solution
    k, hh = vertex_vector_h.shape  # hh = k-1
    theta_0 = vertex_vector_h[0, :]
    m_q = numpy.zeros((k-1, k-1))
    vertex_vector = numpy.zeros(vertex_vector_h.shape, dtype=float)
    vertex_vector[0, :] = theta_0
    max_radius = numpy.zeros(k-1, dtype=float)
    for i in range(1, k):
        theta_i = vertex_vector_h[i, :]
        rand_radius = numpy.abs(theta_i-theta_0)
        max_radius = numpy.max(numpy.vstack([max_radius, rand_radius]), axis=0)
    # print("max_radius ", max_radius)
    for i in range(1, k):
        radius_h = numpy.zeros(k-1, dtype=float)
        radius_h[i-1] = max_radius[i-1]
        vertex_vector[i, :] = theta_0+radius_h

    l_chi_sq = []
    for i in range(0, k):
        theta_i = vertex_vector_h[i, :]
        chi_sq = func(theta_i)
        l_chi_sq.append(chi_sq)
    vertex_chi_sq = numpy.array(l_chi_sq, dtype=float)

    # print("hh, k: ", hh, k)
    # print("theta_0: ", theta_0)
    chi_sq_0 = vertex_chi_sq[0]
    # print("chi_sq_0: ", chi_sq_0)
    v_a = numpy.zeros(k-1)
    m_b = numpy.zeros((k-1, k-1))
    m_q = numpy.zeros((k-1, k-1))
    m_chi_sq_0i = numpy.zeros(k-1)
    # print("step 1")
    for i in range(1, k):
        theta_i = vertex_vector[i, :]
        theta_0i = 0.5*(theta_0+theta_i)
        chi_sq_0i = func(theta_0i)
        # print("ii: {:}     {:}".format(i, chi_sq_0i))
        m_chi_sq_0i[i-1] = chi_sq_0i
        m_q[i-1, :] = theta_i-theta_0

    # print("step 2")
    for i in range(1, k):
        chi_sq_i = vertex_chi_sq[i]
        theta_i = vertex_vector[i, :]
        chi_sq_0i = m_chi_sq_0i[i-1]

        a_i = 4.*chi_sq_0i - chi_sq_i - 3.*chi_sq_0
        v_a[i-1] = a_i

        b_ii = 2.*(chi_sq_i + chi_sq_0 - 2.*chi_sq_0i)
        m_b[i-1, i-1] = b_ii

        for j in range(i+1, k):
            chi_sq_0j = m_chi_sq_0i[j-1]
            theta_j = vertex_vector[j, :]
            theta_ij = 0.5*(theta_i+theta_j)
            chi_sq_ij = func(theta_ij)
            # print("ij: {:} {:}    {:}".format(i, j, chi_sq_ij))
            b_ij = 2.*(chi_sq_ij + chi_sq_0 - chi_sq_0i - chi_sq_0j)
            m_b[i-1, j-1] = b_ij
            m_b[j-1, i-1] = b_ij
    # print("step 3")
    m_ib = numpy.linalg.inv(m_b)
    m_qib = numpy.matmul(m_q, m_ib)
    v_qiba = numpy.matmul(m_qib, v_a)
    # theta_min = theta_0 - v_qiba
    m_qibqt = numpy.matmul(m_qib, m_q.transpose())
    m_error = 2.*chi_sq_0*m_qibqt

    # print("\nm_q")
    # print(m_q)
    # print("\nm_b")
    # print(m_b)
    # print("\nm_ib")
    # print(m_ib)
    # print("\nv_a")
    # print(v_a)
    # print("\ntheta_min: ", theta_min)
    # print("\ntheta_0: ", theta_0)

    # print("\nm_error: ", m_error)
    # print(50*"*")
    return m_error, numpy.abs(v_qiba)
