"""Hermann-Mauguin symbols written without spaces and unknown symbols."""
import pytest

import cryspy
from cryspy.exceptions import RcifFormatError


@pytest.mark.parametrize("name, name_hm_alt, it_number", [
    ("Ia-3d", "I a -3 d", 230),
    ("Fd-3m", "F d -3 m", 227),
    ("Pnma", "P n m a", 62),
    ("P21/c", "P 21/c", 14),
])
def test_name_without_spaces(name, name_hm_alt, it_number):
    space_group = cryspy.SpaceGroup(name_hm_alt=name)
    assert space_group.name_hm_alt == name_hm_alt
    assert space_group.it_number == it_number


def test_name_without_spaces_from_cif():
    space_group = cryspy.SpaceGroup.from_cif(
        "_space_group_name_H-M_alt Ia-3d\n"
        "_space_group_IT_coordinate_system_code 1\n")
    assert space_group.it_number == 230
    assert space_group.centrosymmetry


def test_unknown_name_gives_clear_error():
    with pytest.warns(UserWarning), pytest.raises(RcifFormatError):
        cryspy.SpaceGroup(name_hm_alt="Xyz")
