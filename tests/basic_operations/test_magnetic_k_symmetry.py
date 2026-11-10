"""Group of k, one-dimensional representations, superspace groups, allowed
Fourier components (symmform), orbit splitting and magnetic domains
(A_functions_base.magnetic_k_symmetry)."""
import numpy
import pytest

from cryspy.A_functions_base.magnetic_k_symmetry import (
    calc_axial_moment_matrix, calc_domain_matrices, calc_fourier_site_projector,
    calc_k_group_flags, calc_minus_k_flags, calc_one_dimensional_irreps,
    form_fourier_symmform, form_operation_xyz, form_superspace_group,
    form_symmform_constraint, parse_symmform, split_atom_orbit)
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_superspace_operation_arrays, parse_operation_xyz)

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
    moment = calc_axial_moment_matrix(rotation).reshape(3, 3, -1)
    numpy.testing.assert_array_equal(moment[:, :, 0], numpy.eye(3))
    numpy.testing.assert_array_equal(moment[:, :, 1], numpy.diag([-1, -1, 1]))
    numpy.testing.assert_array_equal(
        moment[:, :, 2], [[0, -1, 0], [1, 0, 0], [0, 0, 1]])


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


def take_arrays(l_operation, k_vector):
    return calc_superspace_operation_arrays(l_operation, [], k_vector)


def calc_symmform(position, l_operation, k_vector):
    arrays = take_arrays(l_operation, k_vector)
    projector = calc_fourier_site_projector(
        position, k_vector, arrays["rotation"], arrays["translation"],
        arrays["moment"], arrays["phase"], arrays["conjugate"])
    return form_fourier_symmform(projector), projector


def test_symmform_on_twofold_axis():
    """P 2 (2 || b), k = (0, 0.3, 0), atom on the axis: the representation
    with tau = 0 allows M || b, the one with time reversal M in the ac plane."""
    k_vector = numpy.array([0., 0.3, 0.])
    position = numpy.array([0., 0.37, 0.])
    symmform, _ = calc_symmform(
        position, ["x1,x2,x3,x4,+1", "-x1,x2,-x3,x4,+1"], k_vector)
    assert symmform == ["0", "myc1", "0", "0", "mys1", "0"]
    symmform, _ = calc_symmform(
        position, ["x1,x2,x3,x4,+1", "-x1,x2,-x3,x4,-1"], k_vector)
    assert symmform == ["mxc1", "0", "mzc1", "mxs1", "0", "mzs1"]
    # the same with tau = 1/2 instead of time reversal
    symmform_tau, _ = calc_symmform(
        position, ["x1,x2,x3,x4,+1", "-x1,x2,-x3,x4+1/2,+1"], k_vector)
    assert symmform_tau == symmform


def test_symmform_with_minus_k_operation():
    """The inversion maps k to -k: at the inversion centre M(-x4) = M(x4),
    so the sin components vanish; shifted by 1/4 in x4 the cos ones do."""
    k_vector = numpy.array([0., 0., 0.31])
    symmform, _ = calc_symmform(
        numpy.zeros(3), ["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4,+1"], k_vector)
    assert symmform == ["mxc1", "myc1", "mzc1", "0", "0", "0"]
    symmform, _ = calc_symmform(
        numpy.zeros(3), ["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4+1/2,+1"], k_vector)
    assert symmform == ["0", "0", "0", "mxs1", "mys1", "mzs1"]


def test_projector_is_idempotent():
    k_vector = numpy.array([0., 0.5, 0.])
    _, projector = calc_symmform(
        numpy.array([0.1, 0.2, 0.3]), ["x1,x2,x3,x4,+1", "-x1,x2+1/2,-x3,x4+1/4,+1"],
        k_vector)
    numpy.testing.assert_allclose(projector @ projector, projector, atol=1e-12)
    assert numpy.trace(projector) == pytest.approx(6.)


def test_parse_symmform_and_constraint():
    assert parse_symmform("-mxs1") == {"mxs1": -1.}
    assert parse_symmform("2mys1+1/2mzs1") == {"mys1": 2., "mzs1": 0.5}
    assert parse_symmform("0") == {}
    constraint, flag_free = form_symmform_constraint(
        ["mxc1", "-mxc1", "0", "mxc1"], numpy.array([0, 0, 0, 1]))
    values = constraint @ numpy.array([1.5, 7., 3., 2.])
    # the symbols are local to the atom (scope): the fourth one is free
    numpy.testing.assert_allclose(values, [1.5, -1.5, 0., 2.])
    numpy.testing.assert_array_equal(flag_free, [True, False, False, True])


def test_superspace_group_of_pnma_as_magndata():
    """Pnma, k = (0, b, 0), Γ1 with the operations mapping k to -k: the
    group of MAGNDATA 1.1.1 (Cs2CuCl4, Pnma1'(0b0)000s)."""
    parent = ["x,y,z", "-x+1/2,-y,z+1/2", "-x,y+1/2,-z", "x+1/2,-y+1/2,-z+1/2",
              "-x,-y,-z", "x+1/2,y,-z+1/2", "x,-y+1/2,z", "-x+1/2,y+1/2,z+1/2"]
    rotation, translation = take_operations(parent)
    k_vector = numpy.array([0., 0.472, 0.])
    flags = calc_k_group_flags(rotation, k_vector)
    i_minus = int(numpy.flatnonzero(calc_minus_k_flags(rotation, k_vector))[0])
    # Γ1: the characters are exp(-2 pi i k.t)
    phase = k_vector @ translation[:, flags]
    ssg = form_superspace_group(
        rotation[:, flags], translation[:, flags], phase,
        NO_CENTRING, k_vector, rotation[:, i_minus], translation[:, i_minus])
    magndata = {"x1,x2,x3,x4,+1", "-x1+1/2,-x2,x3+1/2,-x4,+1",
                "-x1,x2+1/2,-x3,x4,+1", "x1+1/2,-x2+1/2,-x3+1/2,-x4,+1",
                "-x1,-x2,-x3,-x4,+1", "x1+1/2,x2,-x3+1/2,x4,+1",
                "x1,-x2+1/2,x3,-x4,+1", "-x1+1/2,x2+1/2,x3+1/2,x4,+1"}
    assert set(ssg["operations"]) == magndata
    assert ssg["centrings"] == ["x1,x2,x3,x4,+1", "x1,x2,x3,x4+1/2,-1"]
    # the Cu site (MAGNDATA symmform: cos (0, myc1, 0), sin (mxs1, 0, mzs1))
    symmform, _ = calc_symmform(
        numpy.array([0.2302, 0.25, 0.4182]), ssg["operations"], k_vector)
    assert symmform == ["0", "myc1", "0", "mxs1", "0", "mzs1"]


def test_superspace_group_is_not_formed_for_complex_representation():
    """P 4/m, k = (0, 0, g): the group of k is 4 with chi(4) = i; the
    inversion commutes with 4 but maps k to -k, so it needs chi = chi* and
    cannot be added; chi(4) = -1 allows it."""
    rotation, translation = take_operations(["x,y,z", "-y,x,z", "-x,-y,z", "y,-x,z"])
    k_vector = numpy.array([0., 0., 0.31])
    inversion = -numpy.eye(3).flatten()
    ssg = form_superspace_group(
        rotation, translation, numpy.array([0., 0.75, 0.5, 0.25]),  # chi(4) = i
        NO_CENTRING, k_vector, inversion, numpy.zeros(3))
    assert ssg is None
    ssg = form_superspace_group(
        rotation, translation, numpy.array([0., 0.5, 0., 0.5]),  # chi(4) = -1
        NO_CENTRING, k_vector, inversion, numpy.zeros(3))
    assert len(ssg["operations"]) == 8


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


def test_snap_k_vector():
    from cryspy.A_functions_base.magnetic_k_symmetry import snap_k_vector
    k_vector, l_note = snap_k_vector([0.3333, 0.28, 0.5])
    numpy.testing.assert_allclose(k_vector, [1. / 3., 0.28, 0.5])
    assert l_note == ["0.3333 → 1/3"]
    # an incommensurate value far from simple fractions is kept
    assert snap_k_vector([0.0197, 0.0, 0.0])[1] == []


def test_irrep_numbers_and_superposition():
    from cryspy.D_functions_item_loop.function_1_magnetic_k import (
        form_irrep_name, take_irrep_numbers)
    from cryspy.exceptions import CrysPyError
    assert take_irrep_numbers(2, 4) == [2]
    assert take_irrep_numbers("Γ3, 1", 4) == [1, 3]
    assert form_irrep_name([1, 3]) == "Γ1⊕Γ3"
    with pytest.raises(CrysPyError):
        take_irrep_numbers("5", 4)
    with pytest.raises(CrysPyError):
        take_irrep_numbers("1", 0)


def test_classify_moment_modulation():
    """M = Re(A exp(-i t)): A real is collinear, A = (1, -i, 0) rotates in
    the ab plane: a helix for k || c, a cycloid for k || a."""
    from cryspy.D_functions_item_loop.function_1_magnetic_k import (
        classify_moment_modulation)
    cell = numpy.array([4., 4., 4., 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])
    m_min, m_max, kind = classify_moment_modulation(
        numpy.array([2., 0., 0.], dtype=complex), [0., 0., 0.2137], cell)
    assert (m_min, m_max) == pytest.approx((0., 2.), abs=1e-4)
    assert kind.startswith("collinear, amplitude-modulated")
    amplitude = numpy.array([2., -2j, 0.])
    m_min, m_max, kind = classify_moment_modulation(amplitude, [0., 0., 0.2137], cell)
    assert (m_min, m_max) == pytest.approx((2., 2.))
    assert kind.startswith("circular helix")
    assert classify_moment_modulation(amplitude, [0.2137, 0., 0.], cell)[2].startswith(
        "circular cycloid")
    # k = 1/2: only +-A.real
    assert classify_moment_modulation(
        numpy.array([1.5, 0., 0.]), [0.5, 0., 0.], cell)[2] == "collinear, sign alternates"
