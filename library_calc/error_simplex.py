# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:14:11 2019

@author: ikibalin
"""
import numpy

def error_estimation_simplex(vertex_vector, vertex_chi_sq, func):
    """
    Calculations according to 
    ANALYTICAL CHEMISTRY, VOL. 60, NO. 8, APRIL 15, 1988
    
    notations:
    theta_i = vertex_vector[i, :]
    chi_sq_i = vertex_chi_sq[i]
    """
    k, hh = vertex_vector.shape#hh = k-1
    theta_0 = vertex_vector[0, :]
    #print("hh, k: ", hh, k)
    #print("theta_0: ", theta_0)
    chi_sq_0 = vertex_chi_sq[0]
    #print("chi_sq_0: ", chi_sq_0)
    v_a = numpy.zeros(k-1)
    m_b = numpy.zeros((k-1, k-1))
    m_q = numpy.zeros((k-1, k-1))
    m_chi_sq_0i = numpy.zeros(k-1)
    #print("step 1")
    for i in range(1, k):
        theta_i = vertex_vector[i, :]
        theta_0i = 0.5*(theta_0+theta_i)
        chi_sq_0i = func(theta_0i)
        m_chi_sq_0i[i-1] = chi_sq_0i
        
    #print("step 2")
    for i in range(1, k):
        chi_sq_i = vertex_chi_sq[i]
        theta_i = vertex_vector[i, :]
        chi_sq_0i = m_chi_sq_0i[i-1]
        
        a_i = 4.*chi_sq_0i - chi_sq_i - 3.*chi_sq_0
        v_a[i-1] = a_i
        
        b_ii = 2.*(chi_sq_i + chi_sq_0 - 2.*chi_sq_0i)
        m_b[i-1, i-1] = b_ii
        for j in range(1, k):
            theta_j = vertex_vector[j, :]
            chi_sq_0j = m_chi_sq_0i[j-1]

            theta_ij = 0.5*(theta_i+theta_j)
            chi_sq_ij = func(theta_ij)

            b_ij = 2.*(chi_sq_ij + chi_sq_0 - chi_sq_0i - chi_sq_0j)
            m_b[i-1, j-1] = b_ij
            m_b[j-1, i-1] = b_ij
            q_ij = chi_sq_ij - chi_sq_0j
            m_q[i-1, j-1] = q_ij
            m_q[j-1, i-1] = q_ij
    #print("step 3")
    m_ib = numpy.linalg.inv(m_b)
    m_qib = numpy.matmul(m_q, m_ib)
    v_qiba = numpy.matmul(m_qib, v_a)
    theta_min = theta_0 - v_qiba
    #print("\ntheta_min: ", theta_min)
    #print("\ntheta_0: ", theta_0)
    m_qibqt = numpy.matmul(m_qib, m_q.transpose())
    m_error = 2.*chi_sq_0*m_qibqt
    #print("\nm_error: ", m_error)
    #print(50*"*")
    return m_error