"""Group of k, one-dimensional representations, allowed Fourier components,
orbit splitting and magnetic domains (A_functions_base.magnetic_k_symmetry)."""
import numpy
import pytest

from cryspy.A_functions_base.magnetic_k_symmetry import (
    calc_axial_moment_matrix, calc_domain_matrices, calc_fourier_basis,
    calc_fourier_projector, calc_k_group_flags, calc_one_dimensional_irreps,
    form_moment_uvw, form_operation_xyz, split_atom_orbit)
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    parse_operation_xyz)

NO_CENTRING = numpy.zeros((3, 1))


def take_operations(l_xyz):
    l_op = [parse_operation_xyz(xyz) for xyz in l_xyz]
    return (numpy.stack([op[0] for op in l_op], axis=1),
            numpy.stack([op[1] for op in l_op], axis=1))


def test_form_operation_roundtrip():
    for xyz in ("x,y,z", "-y,x-y,z+1/3", "x,y,-z-1/2", "-x+1/2,y+1/2,-z"):
        rotation, translation = parse_operation_xyz(xyz)
        assert form_operation_xyz(rotation, translation) == xyz


def test_axial_moment_matrix():
    rotation, _ = take_operations(["-x,-y,-z", "x,y,-z", "-y,x,z"])
    moment = calc_axial_moment_matrix(rotation)
    assert [form_moment_uvw(m) for m in moment.transpose()] == [
        "u,v,w", "-u,-v,w", "-v,u,w"]


def test_k_group():
    rotation, _ = take_operations(["x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z"])
    flags = calc_k_group_flags(rotation, [0.3, 0., 0.])
    numpy.testing.assert_array_equal(flags, [True, False, False, True])
    # k = (1/2, 0, 0) is equivalent to -k: all operations keep it
    assert calc_k_group_flags(rotation, [0.5, 0., 0.]).all()


def check_irreps(rotation, translation, k_vector, l_chi):
    irreps = calc_one_dimensional_irreps(rotation, translation, NO_CENTRING, k_vector)
    assert len(irreps) == len(l_chi)
    for chi in l_chi:
        assert any(numpy.allclose(chi, irrep) for irrep in irreps)


def test_irreps_screw_axis_at_zone_boundary():
    """P 21: (2_1)^2 = translation (0, 1, 0), exp(-2 pi i k.t) = -1 for
    k = (0, 1/2, 0): the characters of 2_1 are ±i (Izyumov, Naish, Ozerov)."""
    rotation, translation = take_operations(["x,y,z", "-x,y+1/2,-z"])
    check_irreps(rotation, translation, [0., 0.5, 0.], [[1., 1j], [1., -1j]])
    # incommensurate k along the screw axis: chi^2 = exp(-2 pi i k_y)
    k_y = 0.3
    chi = numpy.exp(-1j * numpy.pi * k_y)
    check_irreps(rotation, translation, [0., k_y, 0.], [[1., chi], [1., -chi]])


def test_irreps_of_mmm_group_of_k():
    rotation, translation = take_operations(
        ["x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z",
         "-x,-y,-z", "x,y,-z", "x,-y,z", "-x,y,z"])
    irreps = calc_one_dimensional_irreps(
        rotation, translation, NO_CENTRING, [0.5, 0., 0.])
    assert len(irreps) == 8
    for chi in irreps:
        numpy.testing.assert_allclose(numpy.abs(chi.imag), 0., atol=1e-10)


def test_fourier_components_on_twofold_axis():
    """P 2 (2 || b), atom on the axis: Γ1 allows S || b, Γ2 S in the ac plane."""
    rotation, translation = take_operations(["x,y,z", "-x,y,-z"])
    k_vector = numpy.array([0.5, 0., 0.])
    moment = calc_axial_moment_matrix(rotation)
    xyz = numpy.array([0., 0.37, 0.])
    for phase, n_free, direction in (([0., 0.], 1, [0., 1., 0.]),
                                     ([0., 0.5], 2, None)):
        basis = calc_fourier_basis(calc_fourier_projector(
            xyz, rotation, translation, moment, numpy.array(phase),
            NO_CENTRING, k_vector))
        assert basis.shape[1] == n_free
        if direction is not None:
            numpy.testing.assert_allclose(numpy.abs(basis[:, 0]), direction, atol=1e-10)
        else:
            numpy.testing.assert_allclose(basis[1, :], 0., atol=1e-10)
    # the atom at x = 1/2 is moved by the axis to -1/2 = 1/2 - 1: the lattice
    # phase exp(2 pi i k.L) = -1 exchanges the two cases
    xyz = numpy.array([0.5, 0.37, 0.])
    basis = calc_fourier_basis(calc_fourier_projector(
        xyz, rotation, translation, moment, numpy.array([0., 0.]),
        NO_CENTRING, k_vector))
    assert basis.shape[1] == 2


def test_split_atom_orbit():
    """Immm, k = (1/2, 0, 1/2) keeps all rotations: no splitting. P 4/m m m
    with k = (1/2, 0, 0): the group of k is mmm, the 4-fold orbit of
    (x, 0, 0) (4 positions) splits into (x, 0, 0) and (0, x, 0)."""
    parent_xyz = ["x,y,z", "-x,-y,z", "-y,x,z", "y,-x,z", "-x,y,-z", "x,-y,-z",
                  "y,x,-z", "-y,-x,-z", "-x,-y,-z", "x,y,-z", "y,-x,-z", "-y,x,-z",
                  "x,-y,z", "-x,y,z", "-y,-x,z", "y,x,z"]
    rotation, translation = take_operations(parent_xyz)
    flags = calc_k_group_flags(rotation, [0.5, 0., 0.])
    assert flags.sum() == 8
    l_out = split_atom_orbit(
        numpy.array([0.2, 0., 0.]), rotation, translation,
        rotation[:, flags], translation[:, flags], NO_CENTRING)
    assert len(l_out) == 2
    numpy.testing.assert_allclose(l_out[1][0], [0., 0.2, 0.], atol=1e-10)

    domains = calc_domain_matrices(rotation, rotation[:, flags])
    # two k-domains: k = (1/2, 0, 0) and (0, 1/2, 0)
    assert domains.shape[1] == 2
    numpy.testing.assert_array_equal(domains[:, 0], numpy.eye(3).flatten())


def test_projector_is_idempotent():
    rotation, translation = take_operations(["x,y,z", "-x,y+1/2,-z"])
    moment = calc_axial_moment_matrix(rotation)
    projector = calc_fourier_projector(
        numpy.array([0.1, 0.2, 0.3]), rotation, translation, moment,
        numpy.array([0., 0.25]), NO_CENTRING, numpy.array([0., 0.5, 0.]))
    numpy.testing.assert_allclose(projector @ projector, projector, atol=1e-12)
    assert numpy.trace(projector).real == pytest.approx(3.)
