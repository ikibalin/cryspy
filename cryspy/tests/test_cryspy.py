import pytest
import math
import numpy
from cryspy import (
    Fitable, 
    Crystal, Cell, SpaceGroup, AtomSiteL)

    


def test_init():
    flag = True
    try:
        vv = Fitable()
        vv = Crystal()
        vv = Cell()
        vv = SpaceGroup()
        vv = AtomSiteL()
    except:
        flag = False
    assert flag