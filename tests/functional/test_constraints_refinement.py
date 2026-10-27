"""Refinement with a linear constraint between two parameters (x+y=b)."""
import os

import numpy

import cryspy

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")


def test_refinement_keeps_linear_constraint():
    obj = cryspy.load_file(F_NAME)
    # only the two coordinates: the extinction radius 0 of the file makes
    # the Hessian of all parameters singular
    obj.fix_variables()
    cryspy.add_constraint(obj, "y1 + y2 = 1.5", {
        "y1": "tbte3.atom_site[Te1].fract_y", "y2": "tbte3.atom_site[Te2].fract_y"})
    atom_site = obj.crystal_tbte3.atom_site
    # the start values do not satisfy the constraint
    assert not numpy.isclose(atom_site["Te1"].fract_y + atom_site["Te2"].fract_y, 1.5)

    cryspy.rhochi_rietveld_refinement(obj)

    y_1, y_2 = atom_site["Te1"].fract_y, atom_site["Te2"].fract_y
    assert numpy.isclose(y_1 + y_2, 1.5, atol=1e-10)
    assert abs(y_1 - 0.92893) < 0.01
    # dependent parameter: sigma from the covariance of the free one
    assert atom_site["Te2"].fract_y_sigma > 0.
    assert numpy.isclose(atom_site["Te1"].fract_y_sigma, atom_site["Te2"].fract_y_sigma)


def test_refinement_with_soft_constraint():
    def refine(sigma):
        obj = cryspy.load_file(F_NAME)
        obj.fix_variables()
        cryspy.add_constraint(obj, "y1 + y2 = 1.5", {
            "y1": "tbte3.atom_site[Te1].fract_y", "y2": "tbte3.atom_site[Te2].fract_y"},
            sigma=sigma)
        cryspy.rhochi_rietveld_refinement(obj)
        atom_site = obj.crystal_tbte3.atom_site
        return atom_site["Te1"].fract_y + atom_site["Te2"].fract_y

    # the data alone move the sum to 1.569 (the other parameters are fixed,
    # chi_sq ~ 1e7): a tight restraint keeps the start sum 1.5, a loose one
    # lets the data decide
    assert abs(refine(1e-5) - 1.5) < 1e-3
    assert refine(1e-2) > 1.55
