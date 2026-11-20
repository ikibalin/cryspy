"""Internal parameters of the items (ATTR_INT_NAMES), '<name>_as_string'
and the columns of the loops are descriptors of the classes
(cif_attributes.py) instead of names answered by __getattr__."""
import numpy
import pytest

import cryspy
from cryspy.B_parent_classes.cif_attributes import (
    CifAsString, CifInternal, LoopColumn, LoopNumpyColumn)
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, AtomSiteL
from cryspy.C_item_loop_classes.cl_1_cell import Cell

ATOM_SITE = """loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0.125(3) 0.5 0. .
O1 O 0.25 0.25 0.25 0.5
"""


def test_classes_have_the_descriptors():
    assert isinstance(Cell.__dict__["m_g"], CifInternal)
    assert isinstance(Cell.__dict__["length_a_as_string"], CifAsString)
    assert isinstance(AtomSiteL.__dict__["fract_x"], LoopColumn)
    assert isinstance(AtomSiteL.__dict__["numpy_fract_x"], LoopNumpyColumn)
    # the internal parameters of the items are columns of the loop too
    assert isinstance(AtomSiteL.__dict__["numpy_scat_length_neutron"], LoopNumpyColumn)


def test_internal_parameter_is_calculated_and_dropped_on_change():
    cell = Cell(length_a=5., length_b=5., length_c=5., angle_alpha=90.,
                angle_beta=90., angle_gamma=90.)
    assert "m_g" not in cell.__dict__
    m_g = cell.m_g
    assert numpy.allclose(m_g, numpy.diag([25., 25., 25.]))
    assert "m_g" in cell.__dict__  # kept: the next reading is the dict
    cell.length_a = 4.
    assert "m_g" not in cell.__dict__
    assert cell.m_g[0, 0] == pytest.approx(16.)


def test_internal_parameter_of_undefined_item():
    """An atom without coordinates: no scattering length yet."""
    atom = AtomSite(label="Fe1")
    with pytest.raises(AttributeError, match="scat_length_neutron"):
        atom.scat_length_neutron
    assert not atom.is_attribute("scat_length_neutron")
    assert not hasattr(atom, "scat_length_neutron")


def test_value_as_string():
    loop = AtomSiteL.from_cif(ATOM_SITE)
    fe1, o1 = loop.items
    assert fe1.fract_x_as_string == "0.1250(30)"  # two digits of sigma
    assert o1.fract_x_as_string == o1.take_value_as_string("fract_x")
    assert fe1.occupancy_as_string == "."  # not given
    assert fe1.label_as_string == "Fe1"


def test_loop_columns():
    loop = AtomSiteL.from_cif(ATOM_SITE)
    assert loop.label == ["Fe1", "O1"]
    assert loop.fract_x == [0.125, 0.25]
    assert loop.fract_x_refinement == [True, False]
    numpy_x = loop.numpy_fract_x
    assert numpy_x.dtype == float and numpy.allclose(numpy_x, [0.125, 0.25])
    assert loop.numpy_fract_x_refinement.dtype == bool
    assert loop.numpy_label.dtype.kind == "U"
    # the array is kept: a change of it goes to the items by numpy_to_items
    assert loop.numpy_fract_x is numpy_x
    numpy_x[1] = 0.3
    loop.numpy_to_items()
    assert loop.items[1].fract_x == pytest.approx(0.3)
    assert "numpy_fract_x" not in loop.__dict__
    assert loop.numpy_fract_x[1] == pytest.approx(0.3)


def test_columns_of_unknown_loop():
    obj = cryspy.str_to_globaln("""global_
data_x
loop_
_publ_author_name
_publ_author_id
'A B' 1
'C D' 2
""")
    loop = obj.data_x.items[0]
    assert loop.name == ["A B", "C D"]
    assert loop.numpy_id.dtype.kind == "i"
