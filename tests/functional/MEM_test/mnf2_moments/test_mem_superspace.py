"""MEM for a crystal described by a magnetic superspace group with k = 0.

The superspace description (Crystal.convert_to_superspace) is the same
structure as the magnetic space group: MEM must give the same result.
"""
import os

import numpy
import pytest

import cryspy
from cryspy.procedure_mempy.mempy_by_dictionary import form_mem_crystal_dictionary

F_NAME = os.path.join(os.path.dirname(__file__), "main.rcif")


def take_crystal(obj):
    return [item for item in obj.items if isinstance(item, cryspy.Crystal)][0]


def check_same_magnetic_space_group(dict_msg, dict_mem):
    elems_msg = set(map(tuple, dict_msg["full_mcif_elems"].transpose()))
    elems_mem = set(map(tuple, dict_mem["full_mcif_elems"].transpose()))
    assert elems_msg == elems_mem
    assert numpy.array_equal(dict_msg["atom_multiplicity"], dict_mem["atom_multiplicity"])
    assert list(dict_msg["atom_ordered_label"]) == list(dict_mem["atom_ordered_label"])
    assert numpy.allclose(dict_msg["atom_ordered_moment_crystalaxis_dn"],
                          dict_mem["atom_ordered_moment_crystalaxis_dn"])
    assert numpy.allclose(dict_msg["atom_ordered_sc_m"], dict_mem["atom_ordered_sc_m"])
    assert not any(key.startswith("magnetic_k_") for key in dict_mem.keys())


def test_mem_dictionary_of_superspace_k0():
    dict_msg = take_crystal(cryspy.load_file(F_NAME)).get_dictionary()
    obj = cryspy.load_file(F_NAME)
    take_crystal(obj).convert_to_superspace()
    dict_mem = form_mem_crystal_dictionary(take_crystal(obj).get_dictionary())
    check_same_magnetic_space_group(dict_msg, dict_mem)


def test_mem_dictionary_time_reversal_by_x4_shift():
    """x4 + 1/2 reverses the moments of a single harmonic: the operation
    'x2,x1,x3,x4+1/2,+1' is 'y,x,z,-1' of the magnetic space group."""
    dict_msg = take_crystal(cryspy.load_file(F_NAME)).get_dictionary()
    obj = cryspy.load_file(F_NAME)
    crystal = take_crystal(obj)
    crystal.convert_to_superspace()
    for item in crystal.space_group_symop_magn_ssg_operation.items:
        if item.algebraic.endswith(",-1"):
            item.algebraic = item.algebraic.replace(",x4,-1", ",x4+1/2,+1")
    dict_mem = form_mem_crystal_dictionary(crystal.get_dictionary())
    check_same_magnetic_space_group(dict_msg, dict_mem)


def test_mem_superspace_nonzero_k_kept():
    """k ≠ 0: the moments stay modulated (MEM with the Fourier amplitudes)."""
    obj = cryspy.load_file(F_NAME)
    crystal = take_crystal(obj)
    crystal.convert_to_superspace()
    dict_crystal = crystal.get_dictionary()
    dict_crystal["magnetic_k_vector"] = numpy.array([[0.5], [0.], [0.]])
    assert form_mem_crystal_dictionary(dict_crystal) is dict_crystal


def test_mem_magnetization_density_superspace_k0():
    """Same chi_sq as the magnetic space group description. 2.61 with the
    twin domain (k, h, l), fraction 0.124 of main.rcif; MEM ignored the twin
    before and gave 2.81."""
    obj = cryspy.load_file(F_NAME)
    take_crystal(obj).convert_to_superspace()
    res = cryspy.mempy_magnetization_density_reconstruction(obj)
    assert res["chi_sq"] == pytest.approx(2.6129, abs=1e-3)
