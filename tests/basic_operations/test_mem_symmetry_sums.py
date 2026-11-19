"""Sums over the symmetry elements in MEM (asymmetric unit, chi and moment
channels) compared with direct loops on small arrays."""
import numpy

from cryspy.A_functions_base.local_susceptibility import calc_m_r_inv_m
from cryspy.A_functions_base.mempy import calc_mem_chi
from cryspy.A_functions_base.symmetry_elements import calc_asymmetric_unit_cell_indexes
from cryspy.A_functions_base.unit_cell import calc_eq_ccs_by_unit_cell_parameters

# P 4/m: 4-fold along c, inversion, translations 0 (rows b1, b2, b3,
# denominator, r_11 ... r_33)
R_4Z = numpy.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])


def form_p4m_elements():
    l_r = [numpy.linalg.matrix_power(R_4Z, i) for i in range(4)]
    l_r += [-r for r in l_r]
    return numpy.array([[0, 0, 0, 1] + list(r.flatten()) for r in l_r], dtype=int).transpose()


def test_asymmetric_unit_cell_indexes():
    elems = form_p4m_elements()
    n_abc = numpy.array([8, 8, 6])
    point_index, multiplicity = calc_asymmetric_unit_cell_indexes(n_abc, elems)
    assert multiplicity.sum() == numpy.prod(n_abc)
    # orbits by brute force: the representative has the smallest flat index
    for (i_a, i_b, i_c), mult in zip(point_index.transpose(), multiplicity):
        orbit = set()
        for i_elem in range(elems.shape[1]):
            r = elems[4:, i_elem].reshape(3, 3)
            orbit.add(tuple(numpy.mod(r @ numpy.array([i_a, i_b, i_c]), n_abc)))
        flat = [n_abc[2] * n_abc[1] * p[0] + n_abc[2] * p[1] + p[2] for p in orbit]
        assert min(flat) == n_abc[2] * n_abc[1] * i_a + n_abc[2] * i_b + i_c
        assert len(orbit) == mult


def test_mem_chi_equals_direct_sum():
    rng = numpy.random.default_rng(2)
    elems = form_p4m_elements()
    ucp = numpy.array([4.1, 4.1, 5.3, 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])
    index_hkl = numpy.array([[1, 0, 0], [1, 1, 0], [0, 1, 2], [2, 1, 1]], dtype=int).transpose()
    symm_elem_points = numpy.array([[1, 3, 0, 8], [2, 1, 3, 8], [5, 0, 1, 8]]).transpose()
    n_point = symm_elem_points.shape[1]
    point_susceptibility = rng.normal(size=(9, n_point, 1))
    h_ccs = numpy.array([0.3, -0.2, 0.9])
    multiplicity = numpy.array([8, 8, 4])
    mem = calc_mem_chi(index_hkl, ucp, h_ccs, elems, symm_elem_points, point_susceptibility,
                       70., 2, point_multiplicity=multiplicity, dict_in_out={})

    eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, ucp)[0]
    rs_xyz = calc_m_r_inv_m(ucp, elems)[0]
    coeff = 0.2695 * 70. / (elems.shape[1] * 2)
    expected = numpy.zeros((3, index_hkl.shape[1], n_point, 1), dtype=complex)
    for i_h in range(index_hkl.shape[1]):
        e = eq_ccs[:, i_h]
        for i_p in range(n_point):
            xyz = symm_elem_points[:3, i_p] / symm_elem_points[3, i_p]
            chi = point_susceptibility[:, i_p, 0].reshape(3, 3)
            for i_s in range(elems.shape[1]):
                r_ccs = rs_xyz[:, i_s].reshape(3, 3)
                m = r_ccs @ chi @ numpy.linalg.inv(r_ccs) @ h_ccs
                m_perp = numpy.cross(numpy.cross(e, m), e)
                xyz_s = elems[4:, i_s].reshape(3, 3) @ xyz
                phase = numpy.exp(-2j * numpy.pi * index_hkl[:, i_h] @ xyz_s)
                expected[:, i_h, i_p, 0] += coeff * multiplicity[i_p] * phase * m_perp
    assert numpy.allclose(mem, expected, rtol=1e-10, atol=1e-12)
