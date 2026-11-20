"""ItemN and LoopN are plain Python objects: every attribute is a
descriptor of the class (cif_attributes.py), there is no __getattr__ that
could answer a misspelt name or loop on an object made by copy or pickle."""
import copy
import pickle

import pytest

import cryspy
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, AtomSiteL
from cryspy.C_item_loop_classes.cl_1_cell import Cell

RCIF = """global_
data_fe
_cell_length_a 2.87
_cell_length_b 2.87
_cell_length_c 2.87
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'I m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0 0 0
"""


def test_no_getattr_in_items_and_loops():
    assert "__getattr__" not in ItemN.__dict__
    assert "__getattr__" not in LoopN.__dict__
    assert "__setattr__" not in LoopN.__dict__


def test_misspelt_names_raise_attribute_error():
    cell = Cell(length_a=3.)
    with pytest.raises(AttributeError):
        cell.lenght_a
    assert not hasattr(AtomSite(label="Fe1"), "fract_x")  # not given
    loop = AtomSiteL(items=[AtomSite(label="Fe1", type_symbol="Fe", fract_x=0.,
                                     fract_y=0., fract_z=0.)])
    with pytest.raises(AttributeError):
        loop.fract_xx
    with pytest.raises(AttributeError):
        loop.numpy_fract_xx


def test_objects_made_without_init():
    """copy and pickle make the object first and ask it for __setstate__."""
    for cls in (Cell, AtomSiteL):
        obj = cls.__new__(cls)
        assert not hasattr(obj, "__setstate__") or callable(obj.__setstate__)
        with pytest.raises(AttributeError):
            obj.items if cls is AtomSiteL else obj.length_a
    rcif = cryspy.str_to_globaln(RCIF)
    for obj in (copy.deepcopy(rcif), pickle.loads(pickle.dumps(rcif))):
        assert obj.to_cif() == rcif.to_cif()


def test_items_of_loop_are_checked():
    atom = AtomSite(label="Fe1", type_symbol="Fe", fract_x=0., fract_y=0., fract_z=0.)
    loop = AtomSiteL()
    loop.items = [atom, Cell(length_a=3.)]
    assert loop.items == [atom]
    loop.loop_name = " second "
    assert loop.loop_name == "second"


def test_blocks_by_name():
    rcif = cryspy.str_to_globaln(RCIF)
    crystal = rcif.crystal_fe
    assert rcif.take_item_by_name("CRYSTAL_FE") is crystal
    assert crystal.cell is crystal.take_item_by_name("cell")
    assert crystal.atom_site.label == ["Fe1"]
    with pytest.raises(AttributeError):
        crystal.atom_site_aniso
    assert crystal.take_item_by_name("atom_site_aniso") is None
    # a data block that is not there
    with pytest.raises(AttributeError, match="crystal_cu"):
        rcif.crystal_cu
    assert not hasattr(rcif, "crystal_cu")
    assert rcif.take_item_by_name("crystal_cu") is None
    assert cryspy.GlobalN().global_name is None  # not given
    assert getattr(rcif, "PREFIX_OF_NOTHING", "default") == "default"
