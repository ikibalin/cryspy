"""Kinds of the experiments of a CrysPy object and checks of their models.

The kinds are words that recipes refer to in their 'For:' line, so that
the recipes for the experiments of an object are found by code:

    single_crystal, powder, powder_1d, powder_2d, tof  - the experiment
    polarized, unpolarized, flip_ratio                  - the measured data
    nuclear, magnetic                                   - the crystals of it

A crystal is 'magnetic' when it is described by a magnetic space group
(space_group_symop_magn_operation, MCIF), otherwise 'nuclear'.
"""
from typing import Dict, List, Optional

import numpy

from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import (
    AtomSiteSusceptibility, AtomSiteSusceptibilityL)
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import DiffrnRadiation
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import AtomSiteScat, AtomSiteScatL
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

KINDS = ("single_crystal", "powder", "powder_1d", "powder_2d", "tof",
         "polarized", "unpolarized", "flip_ratio", "nuclear", "magnetic")
FLIP_RATIO_TOLERANCE = 1e-6  # |FR_calc - 1| below it: no magnetic scattering
CHI_START = 0.1  # mu_B/T, a start value of a new susceptibility: zero gives no derivative step
# start values of chi_11 (mu_B/T) tried for new magnetic atoms: the flip ratios
# have several minima of chi_sq in chi, the refinement goes to the nearest one
L_CHI_START = (0.01, 0.03, 0.1, 0.3, 1., 3., 10., 30.)


def take_crystals(cryspy_object) -> List[Crystal]:
    return [item for item in getattr(cryspy_object, "items", []) if isinstance(item, Crystal)]


def take_experiments(cryspy_object) -> List[DataN]:
    return [item for item in getattr(cryspy_object, "items", [])
            if isinstance(item, (Diffrn, Pd, Pd2d, TOF))]


def take_phase_labels(experiment: DataN) -> List[str]:
    if not experiment.is_attribute("phase"):
        return []
    phase = experiment.phase
    if isinstance(phase, LoopN):
        return [str(item.label) for item in phase.items]
    return [str(phase.label)]


def take_experiment_crystals(cryspy_object, experiment: DataN) -> List[Crystal]:
    l_label = [label.lower() for label in take_phase_labels(experiment)]
    return [crystal for crystal in take_crystals(cryspy_object)
            if str(crystal.data_name).lower() in l_label]


def flag_magnetic_crystal(crystal: Crystal) -> bool:
    return crystal.is_attribute("space_group_symop_magn_operation")


def flag_polarized_data(block) -> bool:
    """Measured data with intensities for the two polarizations of the beam."""
    l_item = block.items if isinstance(block, LoopN) else [block]
    if len(l_item) == 0:
        return False
    item = l_item[0]
    return any(item.is_attribute(name) for name in item.ATTR_NAMES if "intensity_plus" in name)


def take_experiment_kinds(cryspy_object, experiment: DataN) -> List[str]:
    """Words of KINDS describing the experiment and its crystals."""
    l_kind = []
    if isinstance(experiment, Diffrn):
        l_kind.append("single_crystal")
        refln = experiment.diffrn_refln
        if len(refln.items) > 0 and refln.items[0].is_attribute("fr"):
            l_kind.extend(["polarized", "flip_ratio"])
        else:
            l_kind.append("unpolarized")
    else:
        if isinstance(experiment, Pd):
            l_kind.extend(["powder", "powder_1d"])
            s_meas = "pd_meas"
        elif isinstance(experiment, Pd2d):
            l_kind.extend(["powder", "powder_2d"])
            s_meas = "pd2d_meas"
        else:
            l_kind.append("tof")
            s_meas = "tof_meas"
        flag_polarized = experiment.is_attribute(s_meas) and flag_polarized_data(
            getattr(experiment, s_meas))
        l_kind.append("polarized" if flag_polarized else "unpolarized")
    for crystal in take_experiment_crystals(cryspy_object, experiment):
        kind = "magnetic" if flag_magnetic_crystal(crystal) else "nuclear"
        if kind not in l_kind:
            l_kind.append(kind)
    return l_kind


def take_object_kinds(cryspy_object) -> Dict[str, List[str]]:
    """{name of the experiment: its kinds}."""
    return {str(experiment.data_name): take_experiment_kinds(cryspy_object, experiment)
            for experiment in take_experiments(cryspy_object)}


def take_flip_ratio_experiments(cryspy_object, name: str = "") -> List[Diffrn]:
    l_experiment = [experiment for experiment in take_experiments(cryspy_object)
                    if "flip_ratio" in take_experiment_kinds(cryspy_object, experiment)]
    if name.strip() != "":
        l_experiment = [experiment for experiment in l_experiment
                        if str(experiment.data_name).lower() == name.strip().lower()]
    return l_experiment


def calc_flip_ratio_deviation(experiment: Diffrn) -> Optional[float]:
    """max |FR_calc - 1| of the reflections which are not excluded; None
    before a calculation."""
    l_deviation = []
    for item in experiment.diffrn_refln.items:
        if not item.is_attribute("fr_calc"):
            continue
        if item.is_attribute("excluded") and bool(item.excluded):
            continue
        l_deviation.append(abs(float(item.fr_calc) - 1.))
    if len(l_deviation) == 0:
        return None
    return float(numpy.max(l_deviation))


def take_field(experiment: Diffrn) -> float:
    setup = experiment.setup
    return float(setup.field) if setup.is_attribute("field") else 0.


def take_polarization(experiment: Diffrn) -> float:
    if not experiment.is_attribute("diffrn_radiation"):
        return 0.
    return float(experiment.diffrn_radiation.polarization)


def take_susceptible_atoms(crystal: Crystal) -> List[str]:
    """Labels of the atoms with a susceptibility tensor and a magnetic form
    factor (atom_site_scat): both are needed for the magnetic scattering."""
    if not (crystal.is_attribute("atom_site_susceptibility") and crystal.is_attribute("atom_site_scat")):
        return []
    l_scat = [str(item.label) for item in crystal.atom_site_scat.items]
    return [str(item.label) for item in crystal.atom_site_susceptibility.items
            if str(item.label) in l_scat]


def find_flip_ratio_problems(cryspy_object, experiment: Diffrn) -> List[str]:
    """Reasons why the model gives no magnetic scattering (FR_calc = 1)."""
    l_problem = []
    name = experiment.data_name
    if take_field(experiment) == 0.:
        l_problem.append(f"the magnetic field of {name:} (setup.field) is zero")
    if take_polarization(experiment) == 0.:
        l_problem.append(f"the polarization of the beam of {name:} (diffrn_radiation) is zero")
    for crystal in take_experiment_crystals(cryspy_object, experiment):
        if len(take_susceptible_atoms(crystal)) == 0:
            l_problem.append(f"no atom of the crystal {crystal.data_name:} has a susceptibility "
                             "tensor with a magnetic form factor (atom_site_susceptibility and "
                             "atom_site_scat)")
    return l_problem


def set_field(experiment: Diffrn, field: float):
    experiment.setup.field = field


def set_polarization(experiment: Diffrn, polarization: float):
    if experiment.is_attribute("diffrn_radiation"):
        experiment.diffrn_radiation.polarization = polarization
    else:
        experiment.add_items([DiffrnRadiation(polarization=polarization, efficiency=1.)])


def add_susceptibility(crystal: Crystal, label: str, lande: float) -> str:
    """Make the atom paramagnetic: isotropic susceptibility (Ciso) and the
    magnetic form factor of its type with the Lande factor."""
    if label not in [str(item.label) for item in crystal.atom_site.items]:
        raise ValueError(f"atom '{label:}' is not in the crystal {crystal.data_name:}.")
    ls_out = []
    if not crystal.is_attribute("atom_site_susceptibility"):
        crystal.add_items([AtomSiteSusceptibilityL()])
    l_chi = crystal.atom_site_susceptibility.items
    if label not in [str(item.label) for item in l_chi]:
        l_chi.append(AtomSiteSusceptibility(
            label=label, chi_type="Ciso", chi_11=CHI_START, chi_22=CHI_START, chi_33=CHI_START,
            chi_12=0., chi_13=0., chi_23=0.))
        ls_out.append(f"susceptibility Ciso chi_11 = {CHI_START:} mu_B/T")
    if not crystal.is_attribute("atom_site_scat"):
        crystal.add_items([AtomSiteScatL()])
    l_scat = crystal.atom_site_scat.items
    if label not in [str(item.label) for item in l_scat]:
        l_scat.append(AtomSiteScat(label=label, lande=lande, kappa=1.))
        ls_out.append(f"magnetic form factor with Lande factor {lande:} and kappa 1")
    crystal.apply_constraints()
    if len(ls_out) == 0:
        return f"{crystal.data_name:}: atom {label:} is magnetic already."
    return f"{crystal.data_name:}: atom {label:} gets " + " and ".join(ls_out) + "."


def set_chi_type(cryspy_object, chi_type: str, atoms: str = "") -> List[str]:
    """Set the type of the susceptibility tensor ('Ciso' or 'Cani') of the
    atoms (all atoms with a susceptibility when empty, else labels separated
    by ',' or spaces); returns 'crystal: atom' for each changed atom."""
    l_label = [label.lower() for label in atoms.replace(",", " ").split()]
    l_changed = []
    for crystal in take_crystals(cryspy_object):
        if not crystal.is_attribute("atom_site_susceptibility"):
            continue
        for item in crystal.atom_site_susceptibility.items:
            if len(l_label) > 0 and str(item.label).lower() not in l_label:
                continue
            if not (item.is_attribute("chi_type") and item.chi_type == chi_type):
                item.chi_type = chi_type
                l_changed.append(f"{crystal.data_name:}: {item.label:}")
        crystal.apply_constraints()
    return l_changed


def scan_chi_start(cryspy_object, l_item: List[AtomSiteSusceptibility]) -> float:
    """Set chi_11 of the isotropic tensors l_item (all to the same value)
    to the start value of L_CHI_START with the smallest chi_sq; returns it."""
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
    l_chi_sq = []
    for chi in L_CHI_START:
        for item in l_item:
            item.chi_11 = chi
        for crystal in take_crystals(cryspy_object):
            crystal.apply_constraints()
        chi_sq = float(rhochi_no_refinement(cryspy_object)["chi_sq"])
        l_chi_sq.append(chi_sq if numpy.isfinite(chi_sq) else numpy.inf)
    chi_best = L_CHI_START[int(numpy.argmin(l_chi_sq))]
    for item in l_item:
        item.chi_11 = chi_best
    for crystal in take_crystals(cryspy_object):
        crystal.apply_constraints()
    return chi_best
