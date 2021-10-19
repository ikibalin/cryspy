import os
import math
import numpy

from cryspy.A_functions_base.symmetry_elements import \
    calc_symm_flags, \
    calc_multiplicity_by_atom_symm_elems

na = numpy.newaxis

symm_elems = numpy.array([
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1, 1, 1, 1],
    [1,-1,-1, 1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1,-1, 1,-1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1,-1, 1,-1]], dtype=int)

atom_symm_elems_1 = numpy.array([
    [1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0],
    [0, 2, 0, 0, 0],
    [4, 3, 2, 1, 1],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0]], dtype=int)

atom_symm_elems_2 = numpy.array([
    [1, 0, 0, 0, 0],
    [0, 1, 2, 0, 0],
    [0, 2, 0, 0, 0],
    [4, 3, 4, 2, 2]], dtype=int)

symm_flags = numpy.array(
    [[ True,  True,  True,  True,  True],
     [False, False,  True,  True,  True],
     [False,  True,  True,  True,  True],
     [ True, False,  True,  True,  True]], dtype=bool)

multiplicity = numpy.array([2, 2, 1, 1, 1], dtype=int)


def test_calc_symm_flags():
    res = calc_symm_flags(symm_elems[:, :, na], atom_symm_elems_1[:, na, :])
    print("res: ", res)
    print(symm_flags)

    assert numpy.all(numpy.isclose(res, symm_flags))

    res = calc_symm_flags(symm_elems[:, :, na], atom_symm_elems_2[:, na, :])
    print("res: ", res)
    print(symm_flags)

    assert numpy.all(numpy.isclose(res, symm_flags))


def test_calc_multiplicity_by_atom_symm_elems():
    res = calc_multiplicity_by_atom_symm_elems(symm_elems, atom_symm_elems_1)
    print("res: ", res)
    print(multiplicity)

    assert numpy.all(numpy.isclose(res, multiplicity))

    res = calc_multiplicity_by_atom_symm_elems(symm_elems, atom_symm_elems_2)
    print("res: ", res)
    print(multiplicity)

    assert numpy.all(numpy.isclose(res, multiplicity))

