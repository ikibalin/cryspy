"""Magnetic structure factors of incommensurate MAGNDATA structures against a
direct sum over the moments of a large crystal.

The .mcif files in ``magndata/`` are unchanged entries 1.1.1, 1.1.2,
1.1.16, 1.1.19, 1.1.48, 1.1.77 and 1.1.84 of the Bilbao MAGNDATA database
(https://www.cryst.ehu.es/magndata, copies from
github.com/hgheiberger/quantum-structure-ml, prop_vector/data). They cover
the operations with epsilon = -1, the phases tau = 1/4, 1/3, 1/2, the terms
h.x in x4 (``-x1-x2+x4+1/2``, ``3x3-x4``), centrings with time reversal,
the rhombohedral centring in hexagonal axes and symmform constraints.

The reference does not use CrysPy: it reads the superspace operations by
itself, gives every atom of N1 x N2 x N3 cells its moment with the msCIF
rule ``M_nu(h.x + epsilon x4 + tau) = theta det(R) R M_mu(x4)`` (checking
that all operations giving the same site give the same moment) and sums
``M exp(2 pi i H.r)`` over the crystal. k is replaced by a nearby rational
value so that the crystal of N cells is periodic.
"""
import os
import warnings
from fractions import Fraction

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_f_m_perp_k_by_dictionary,
    calc_moment_modulation_by_dictionary,
)
from magnetic_k_reference import calc_f_m_perp, generate_moments, take_cells

DIR = os.path.join(os.path.dirname(__file__), "magndata")

# entry, rational k close to the published one
CASES = [
    ("1.1.1", (0, Fraction(9, 20), 0)),
    ("1.1.2", (Fraction(1, 3), Fraction(1, 3), Fraction(1, 2))),
    ("1.1.16", (Fraction(1, 2), Fraction(1, 2), Fraction(3, 10))),
    ("1.1.19", (Fraction(5, 9), 0, 0)),
    ("1.1.48", (0, Fraction(1, 8), 1)),
    ("1.1.77", (Fraction(1, 10), Fraction(1, 10), Fraction(3, 2))),
    ("1.1.84", (0, Fraction(3, 4), Fraction(1, 2))),
]


def set_wave_vector(crystal, k_vector):
    item = crystal.cell_wave_vector.items[0]
    item.x, item.y, item.z = [float(value) for value in k_vector]


def calc_reference(crystal, k_vector, index_hkl):
    """F_M_perp (3, n_hkl) by the direct sum, form factor 1, B = 0."""
    atoms = {}
    for item in crystal.atom_site_moment_fourier.items:
        if item.atom_site_label not in atoms:
            site = crystal.atom_site[item.atom_site_label]
            atoms[item.atom_site_label] = (
                (site.fract_x, site.fract_y, site.fract_z), numpy.zeros((3, 2)))
        atoms[item.atom_site_label][1]["xyz".index(item.axis)] = (item.cos, item.sin)
    n_cell = take_cells(k_vector)
    sites = generate_moments(
        [item.algebraic for item in crystal.space_group_symop_magn_ssg_operation.items],
        [item.algebraic for item in crystal.space_group_symop_magn_ssg_centering.items],
        k_vector, atoms, n_cell)
    cell = crystal.cell
    return calc_f_m_perp(
        sites, n_cell, (cell.length_a, cell.length_b, cell.length_c),
        (cell.angle_alpha, cell.angle_beta, cell.angle_gamma), index_hkl)


@pytest.mark.parametrize("entry, k_vector", CASES)
def test_magndata_structure_factor(entry, k_vector):
    with open(os.path.join(DIR, f"magndata_{entry}.mcif")) as fid:
        text = fid.read()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fourier = cryspy.str_to_globaln(text).items[0].atom_site_moment_fourier
    labels = sorted({item.atom_site_label for item in fourier.items})
    text_scat = "\nloop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n" \
        "_atom_site_scat_lande\n" + "".join(f"{label} Fe3+ 2.0\n" for label in labels)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        crystal = cryspy.str_to_globaln(text + text_scat).items[0]
        set_wave_vector(crystal, k_vector)
        dict_crystal = crystal.get_dictionary()
    # form factor 1 and B = 0, as in the reference
    dict_crystal["mag_atom_j0_parameters"][:] = 0.
    dict_crystal["mag_atom_j0_parameters"][6] = 1.
    dict_crystal["mag_atom_j2_parameters"][:] = 0.
    dict_crystal["atom_b_iso"][:] = 0.

    k_float = numpy.array([float(value) for value in k_vector])
    rng = numpy.random.default_rng(1)
    index_h = rng.integers(-3, 4, size=(3, 10)).astype(float)
    index_hkl = numpy.concatenate(
        [index_h + k_float[:, numpy.newaxis], index_h - k_float[:, numpy.newaxis]], axis=1)

    f_m_perp = calc_f_m_perp_k_by_dictionary(
        dict_crystal, index_hkl, dict_crystal["unit_cell_parameters"])
    f_reference = calc_reference(crystal, k_vector, index_hkl)
    assert numpy.abs(f_reference).max() > 0.05
    # rounded published moments break the symmetry slightly: 1e-3 of 10^-12 cm
    assert numpy.allclose(f_m_perp, f_reference, atol=1e-3)


@pytest.mark.parametrize("entry, k_vector", CASES)
def test_magndata_moments_in_real_space(entry, k_vector):
    """calc_moment_modulation_by_dictionary gives the moments of the
    direct generation in every cell: M(p + L) = Re(A exp(-2 pi i k.L))."""
    with open(os.path.join(DIR, f"magndata_{entry}.mcif")) as fid:
        text = fid.read()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fourier = cryspy.str_to_globaln(text).items[0].atom_site_moment_fourier
    labels = sorted({item.atom_site_label for item in fourier.items})
    text_scat = "\nloop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n" \
        "_atom_site_scat_lande\n" + "".join(f"{label} Fe3+ 2.0\n" for label in labels)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        crystal = cryspy.str_to_globaln(text + text_scat).items[0]
        set_wave_vector(crystal, k_vector)
        modulation = calc_moment_modulation_by_dictionary(crystal.get_dictionary())

    atoms = {}
    for item in crystal.atom_site_moment_fourier.items:
        if item.atom_site_label not in atoms:
            site = crystal.atom_site[item.atom_site_label]
            atoms[item.atom_site_label] = (
                (site.fract_x, site.fract_y, site.fract_z), numpy.zeros((3, 2)))
        atoms[item.atom_site_label][1]["xyz".index(item.axis)] = (item.cos, item.sin)
    sites = generate_moments(
        [item.algebraic for item in crystal.space_group_symop_magn_ssg_operation.items],
        [item.algebraic for item in crystal.space_group_symop_magn_ssg_centering.items],
        k_vector, atoms, take_cells(k_vector))
    k_float = numpy.array([float(value) for value in k_vector])
    n_site = 0
    for l_site in sites.values():
        for position, moment in l_site:
            position = numpy.array([float(value) for value in position])
            cell = numpy.floor(position + 1e-9)
            flags = numpy.all(numpy.isclose(
                modulation["fract_xyz"], (position - cell)[:, numpy.newaxis], atol=1e-4), axis=0)
            assert flags.sum() == 1
            calculated = (modulation["amplitude"][:, flags][:, 0]
                          * numpy.exp(-2j * numpy.pi * k_float @ cell)).real
            numpy.testing.assert_allclose(calculated, numpy.array(moment, dtype=float), atol=1e-4)
            n_site += 1
    assert n_site == len(modulation["label"]) * int(numpy.prod(take_cells(k_vector)))
