"""Attributes of the CIF item classes are descriptors made from the ATTR_*
tuples (cryspy/B_parent_classes/cif_attributes.py)."""
import copy
import pickle
import warnings

import numpy
import pytest

import cryspy
from cryspy.B_parent_classes.cif_attributes import (
    CifFlag, CifMark, CifPlain, CifSigma, CifValue)
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.exceptions import CrysPyWarning, RcifFormatError


def take_item_classes() -> list:
    l_class, l_new = [], list(ItemN.__subclasses__())
    while l_new:
        cls = l_new.pop()
        l_class.append(cls)
        l_new.extend(cls.__subclasses__())
    return l_class


ITEM_CLASSES = take_item_classes()


@pytest.mark.parametrize("cls", ITEM_CLASSES, ids=lambda cls: cls.__name__)
def test_every_attribute_is_a_descriptor(cls):
    assert len(cls.ATTR_NAMES) == len(cls.ATTR_TYPES), "a type for every name"
    for name, value_type in zip(cls.ATTR_NAMES, cls.ATTR_TYPES):
        descriptor = getattr(cls, name)
        assert isinstance(descriptor, CifValue) and descriptor.value_type is value_type
    for names, kind in ((cls.ATTR_SIGMA, CifSigma), (cls.ATTR_REF_FLAG, CifFlag),
                        (cls.ATTR_CONSTR_FLAG, CifFlag), (cls.ATTR_CONSTR_MARK, CifMark)):
        for name in names:
            assert isinstance(getattr(cls, name), kind)
    for name in cls.ATTR_INT_PROTECTED_NAMES:
        assert isinstance(getattr(cls, name), CifPlain)


def test_type_limits_and_none():
    cell = cryspy.Cell(length_a=3)
    assert cell.length_a == 3.0 and isinstance(cell.length_a, float)
    cell.length_a = "4.5"
    assert cell.length_a == 4.5
    cell.angle_alpha = 200  # D_MAX
    assert cell.angle_alpha == 180.
    cell.length_b = numpy.array([5.])
    assert cell.length_b == 5. and not isinstance(cell.length_b, numpy.ndarray)
    cell.length_c = "?"
    assert cell.length_c is None and cell.is_attribute("length_c")
    cell.length_a_sigma = "0.1"
    assert cell.length_a_sigma == 0.1
    cell.length_a_refinement = 1
    assert cell.length_a_refinement is True
    with pytest.raises(RcifFormatError, match="'abc' of 'length_a' in Cell is not float"):
        cell.length_a = "abc"


def test_not_given_attribute():
    cell = cryspy.Cell()
    del cell.length_b
    assert not cell.is_attribute("length_b")
    with pytest.raises(AttributeError, match="length_b"):
        cell.length_b
    assert "length_b" not in cell.to_cif()


def test_allowed_values():
    setup = cryspy.Setup(radiation="neutrons")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(CrysPyWarning, match="allowed values"):
            setup.radiation = "light"
    assert setup.radiation == "neutrons"


def test_internal_parameters_are_dropped_by_a_change():
    cell = cryspy.Cell(length_a=3, length_b=4, length_c=5)
    volume = cell.volume
    cell.length_a = 6
    assert "volume" not in cell.__dict__
    assert cell.volume == pytest.approx(2 * volume)


def test_copy_and_pickle():
    cell = cryspy.Cell(length_a=3, length_a_sigma=0.1, length_a_refinement=True)
    for copied in (copy.deepcopy(cell), pickle.loads(pickle.dumps(cell))):
        assert copied.to_cif() == cell.to_cif()
        copied.length_a = 4
        assert cell.length_a == 3.
