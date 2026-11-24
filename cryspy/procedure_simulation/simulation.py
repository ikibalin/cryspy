import os
import numpy





from cryspy.C_item_loop_classes.cl_1_phase import Phase, PhaseL
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeas, PdMeasL
from cryspy.C_item_loop_classes.cl_1_pd_background import PdBackground, PdBackgroundL
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnRefln, DiffrnReflnL
from cryspy.A_functions_base.unit_cell import calc_m_b_by_unit_cell_parameters
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal

from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
from cryspy.procedure_rhochi.rhochi_by_dictionary import rhochi_calc_chi_sq_by_dictionary

from cryspy.A_functions_base.structure_factor import \
    calc_index_hkl_multiplicity_in_range
from cryspy.A_functions_base.structure_factor_magnetic_k import \
    calc_unique_satellite_index_hkl_by_dictionary
        
import cryspy
from cryspy.exceptions import CrysPyError
from cryspy.A_functions_base.magnetic_k_symmetry import snap_k_vector

na = numpy.newaxis

DIR = os.path.dirname(__file__)
# more profile points make the simulation slow and the rcif file huge
N_POINT_MAX = 100000


def simulation_paramagnetic_crystal(cryspy_object: cryspy.GlobalN):
    """Add example of paramagnetic crystal to the cryspy object
    """
    f_crystal_para = os.path.join(DIR, "crystal_Ho2Ti2O7.rcif")
    rcif_crystal = cryspy.load_file(f_crystal_para).crystal_ho2Ti2o7
    cryspy_object.add_items([rcif_crystal, ])


def simulation_magnetic_crystal(cryspy_object: cryspy.GlobalN):
    """Add example of magnetic crystal to the cryspy object
    """
    f_crystal_mag = os.path.join(DIR, "crystal_URu0.96Rh0.04Si2.rcif")
    rcif_crystal = cryspy.load_file(f_crystal_mag).crystal_ururhsi
    cryspy_object.add_items([rcif_crystal, ])


def check_simulation_range(wavelength: float, ttheta_min: float, ttheta_max: float):
    """Raise CrysPyError for a wavelength or a 2theta range that cannot be
    measured."""
    if not wavelength > 0.:
        raise CrysPyError(f"The wavelength must be positive, not {wavelength:}.")
    if not 0. <= ttheta_min < ttheta_max <= 180.:
        raise CrysPyError(
            f"The 2theta range must be within 0..180 degrees and its minimum below "
            f"its maximum, not {ttheta_min:}..{ttheta_max:}.")


def take_crystal_names(cryspy_object: cryspy.GlobalN) -> list:
    """Names of the crystals (Crystal and MagCrystal) of the object."""
    l_name = [item.data_name for item in cryspy_object.items if isinstance(item, cryspy.Crystal)]
    if len(l_name) == 0:
        raise CrysPyError("The object has no crystal to simulate the experiment for.")
    return l_name


def simulate_powder(
        cryspy_object: cryspy.GlobalN, file_name: str, wavelength: float, ttheta_min: float,
        ttheta_max: float, ttheta_step: float, flag_polarized: bool):
    """Add a powder experiment of the template file with the measured
    profile replaced by the calculated one at the given 2theta points
    (degrees) and wavelength (Angstrom); the flat background is kept."""
    check_simulation_range(wavelength, ttheta_min, ttheta_max)
    if not ttheta_step > 0.:
        raise CrysPyError(f"The 2theta step must be positive, not {ttheta_step:}.")
    n_point = int(numpy.floor((ttheta_max - ttheta_min) / ttheta_step + 1e-6)) + 1
    if n_point > N_POINT_MAX:
        raise CrysPyError(
            f"{n_point:} points of the profile: too many, increase the 2theta step.")
    l_name = take_crystal_names(cryspy_object)

    rcif_pd = cryspy.load_file(os.path.join(DIR, file_name)).items[0]
    rcif_pd.setup.wavelength = float(wavelength)
    rcif_pd.range.ttheta_min = float(ttheta_min)
    rcif_pd.range.ttheta_max = float(ttheta_max)
    background = float(numpy.mean(rcif_pd.pd_background.intensity))
    loop_background = PdBackgroundL()
    loop_background.items = [PdBackground(ttheta=float(ttheta), intensity=background)
                             for ttheta in (ttheta_min, ttheta_max)]
    rcif_pd.pd_background = loop_background
    loop_phase = PhaseL()
    loop_phase.items = [Phase(label=name, scale=1.) for name in l_name]
    rcif_pd.phase = loop_phase
    np_tth = numpy.round(ttheta_min + ttheta_step * numpy.arange(n_point), 6)
    # placeholder values, replaced by the calculated profile below
    if flag_polarized:
        l_item_meas = [PdMeas(ttheta=float(tth), intensity_plus=1., intensity_plus_sigma=1.,
                              intensity_minus=1., intensity_minus_sigma=1.) for tth in np_tth]
    else:
        l_item_meas = [PdMeas(ttheta=float(tth), intensity=1., intensity_sigma=1.)
                       for tth in np_tth]
    loop_meas = PdMeasL()
    loop_meas.items = l_item_meas
    rcif_pd.pd_meas = loop_meas
    cryspy_object.add_items([rcif_pd, ])
    rhochi_no_refinement(cryspy_object)

    np_tth = rcif_pd.pd_proc.numpy_ttheta
    np_plus_net = rcif_pd.pd_proc.numpy_intensity_plus_net
    np_minus_net = rcif_pd.pd_proc.numpy_intensity_minus_net
    np_bkg_calc = rcif_pd.pd_proc.numpy_intensity_bkg_calc
    np_sigma = numpy.ones_like(np_plus_net) * 0.1

    l_item_meas = []
    for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
        if flag_polarized:
            item_meas = PdMeas(
                ttheta=tth, intensity_plus=plus + 0.5*bkg_calc, intensity_plus_sigma=sigma,
                intensity_minus=minus + 0.5*bkg_calc, intensity_minus_sigma=sigma,)
        else:
            item_meas = PdMeas(ttheta=tth, intensity=plus + minus + bkg_calc, intensity_sigma=sigma,)
        l_item_meas.append(item_meas)
    loop_meas = PdMeasL()
    loop_meas.items = l_item_meas
    rcif_pd.pd_meas = loop_meas
    rhochi_no_refinement(cryspy_object)


def simulation_polarized_neutron_powder_diffraction(
        cryspy_object: cryspy.GlobalN, wavelength: float = 1.1, ttheta_min: float = 4.,
        ttheta_max: float = 80., ttheta_step: float = 0.2):
    """Simulate a polarized neutron powder diffraction experiment.

    A data block ``pd_pnpd`` with all crystals of the object as phases is
    added; its measured profiles (up and down polarization) are the
    calculated ones on a flat background. Magnetic satellites of crystals
    with a propagation vector are included.

    Parameters
    ----------
    cryspy_object : GlobalN
        Object with the crystals.
    wavelength : float
        Neutron wavelength in Angstrom.
    ttheta_min : float
        Minimal 2theta of the profile in degrees.
    ttheta_max : float
        Maximal 2theta of the profile in degrees.
    ttheta_step : float
        Step of the profile in 2theta in degrees.
    """
    simulate_powder(cryspy_object, "pd_pnpd.rcif", wavelength, ttheta_min, ttheta_max,
                    ttheta_step, True)


def simulation_unpolarized_neutron_powder_diffraction(
        cryspy_object: cryspy.GlobalN, wavelength: float = 1.56, ttheta_min: float = 2.,
        ttheta_max: float = 100., ttheta_step: float = 0.05):
    """Simulate an unpolarized neutron powder diffraction experiment.

    A data block ``pd_unpd`` with all crystals of the object as phases is
    added; its measured profile is the calculated one on a flat
    background. Magnetic satellites of crystals with a propagation vector
    are included.

    Parameters
    ----------
    cryspy_object : GlobalN
        Object with the crystals.
    wavelength : float
        Neutron wavelength in Angstrom.
    ttheta_min : float
        Minimal 2theta of the profile in degrees.
    ttheta_max : float
        Maximal 2theta of the profile in degrees.
    ttheta_step : float
        Step of the profile in 2theta in degrees.
    """
    simulate_powder(cryspy_object, "pd_unpd.rcif", wavelength, ttheta_min, ttheta_max,
                    ttheta_step, False)


def has_k_vector(d_crystal: dict, k_vector) -> bool:
    """True when the crystal dictionary has the propagation vector k (or
    -k, or k shifted by a reciprocal lattice vector)."""
    if "magnetic_k_vector" not in d_crystal:
        return False
    for vector in numpy.asarray(d_crystal["magnetic_k_vector"], dtype=float).transpose():
        for sign in (1., -1.):
            difference = sign * vector - numpy.asarray(k_vector, dtype=float)
            if numpy.allclose(difference, numpy.round(difference), atol=1e-4):
                return True
    return False


def calc_satellites_of_k(crystal, k_vector, sthovl_min: float, sthovl_max: float,
                         unit_cell_parameters) -> numpy.ndarray:
    """Satellites h ± k (3, n_hkl) in the range, one per orbit of the group
    of k of the crystal (and Friedel's law), for a k that is not in the
    crystal."""
    k_group = crystal.take_group_of_k(k_vector)
    d_k = {"magnetic_k_op_rotation": numpy.asarray(k_group["rotation"]),
           "magnetic_k_centring": numpy.asarray(k_group["centring"], dtype=float),
           "magnetic_k_vector": numpy.asarray(k_vector, dtype=float)[:, numpy.newaxis]}
    return calc_unique_satellite_index_hkl_by_dictionary(
        d_k, sthovl_min, sthovl_max, unit_cell_parameters)


def simulation_single_crystal(
        cryspy_object: cryspy.GlobalN, wavelength: float = 1.4, ttheta_min: float = 1.,
        ttheta_max: float = 116., satellites_only: bool = False, k_1: float = 0.,
        k_2: float = 0., k_3: float = 0.):
    """Simulate an unpolarized single crystal experiment.

    For every crystal a data block ``diffrn_sc`` with unpolarized
    intensities of the reflections in the 2theta range is added; a crystal
    with a propagation vector also gets its magnetic satellites h ± k
    (fractional indices, one per symmetry orbit). A propagation vector
    k_1, k_2, k_3 (not zero) adds its satellites too, also to a
    paramagnetic crystal or to a crystal without this k: their intensities
    are then zero, but the list of the satellites in the range (indices,
    angles) is ready for planning an experiment or for the measured data.

    Parameters
    ----------
    cryspy_object : GlobalN
        Object with the crystals.
    wavelength : float
        Neutron wavelength in Angstrom.
    ttheta_min : float
        Minimal 2theta of the reflections in degrees.
    ttheta_max : float
        Maximal 2theta of the reflections in degrees.
    satellites_only : bool
        The magnetic satellites only, without the nuclear reflections; for
        crystals with a propagation vector or with k given.
    k_1 : float
        First component of a propagation vector whose satellites are added
        (reciprocal lattice units); 0, 0, 0 adds none.
    k_2 : float
        Second component of the propagation vector (reciprocal lattice units).
    k_3 : float
        Third component of the propagation vector (reciprocal lattice units).

    Returns
    -------
    str
        Numbers of the simulated reflections and satellites per crystal.
    """
    check_simulation_range(wavelength, ttheta_min, ttheta_max)
    k_extra = None
    if any(abs(float(value)) > 1e-6 for value in (k_1, k_2, k_3)):
        k_extra = snap_k_vector([float(k_1), float(k_2), float(k_3)])[0]
    take_crystal_names(cryspy_object)
    f_diffrn_sc = os.path.join(DIR, "diffrn_sc.rcif")

    l_diffrn, ls_out = [], []
    i_crystal = 0
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            i_crystal += 1
            d_crystal = item.get_dictionary()

            rcif_diffrn = cryspy.load_file(f_diffrn_sc).diffrn_sc
            if i_crystal == 1:
                rcif_diffrn.data_name = "sc"
            else:
                rcif_diffrn.data_name = f"sc{i_crystal}"
            rcif_diffrn.phase.label = item.data_name
            rcif_diffrn.setup.wavelength = float(wavelength)

            unit_cell_parameters = d_crystal["unit_cell_parameters"]
            
            m_b = calc_m_b_by_unit_cell_parameters(unit_cell_parameters)[0]
            reduced_symm_elems = d_crystal.get("reduced_symm_elems", numpy.array([0,0,0,1,1,0,0,0,1,0,0,0,1],dtype=int))
            translation_elems = d_crystal.get("translation_elems", numpy.array([0,0,0,1],dtype=int))
            centrosymmetry = d_crystal.get("centrosymmetry", True)
            sthovl_min = numpy.sin(numpy.radians(0.5 * ttheta_min)) / wavelength
            sthovl_max = numpy.sin(numpy.radians(0.5 * ttheta_max)) / wavelength
            flag_k = "magnetic_k_vector" in d_crystal
            flag_extra = k_extra is not None and not has_k_vector(d_crystal, k_extra)
            if satellites_only and not (flag_k or flag_extra):
                raise CrysPyError(
                    f"The crystal '{item.data_name:}' has no propagation vector, "
                    "so it has no magnetic satellites.")
            np_hkl = numpy.zeros((3, 0), dtype=float)
            if not satellites_only:
                np_hkl = calc_index_hkl_multiplicity_in_range(
                    sthovl_min=sthovl_min, 
                    sthovl_max=sthovl_max, unit_cell_parameters=unit_cell_parameters, 
                    reduced_symm_elems=reduced_symm_elems, translation_elems=translation_elems, 
                    centrosymmetry=centrosymmetry)[0].astype(float)
            if flag_k:
                # magnetic satellites h ± k of a structure with a propagation vector
                np_hkl = numpy.concatenate([
                    np_hkl,
                    calc_unique_satellite_index_hkl_by_dictionary(
                        d_crystal, sthovl_min, sthovl_max, unit_cell_parameters)],
                    axis=1)
            if flag_extra:
                # satellites of a k the crystal does not have: zero intensities
                np_hkl = numpy.concatenate([
                    np_hkl, calc_satellites_of_k(item, k_extra, sthovl_min, sthovl_max,
                                                 unit_cell_parameters)], axis=1)

            rcif_diffrn.diffrn_orient_matrix.ub_11 = float(numpy.round(m_b[0],6))
            rcif_diffrn.diffrn_orient_matrix.ub_12 = float(numpy.round(m_b[1],6))
            rcif_diffrn.diffrn_orient_matrix.ub_13 = float(numpy.round(m_b[2],6))
            rcif_diffrn.diffrn_orient_matrix.ub_21 = float(numpy.round(m_b[3],6))
            rcif_diffrn.diffrn_orient_matrix.ub_22 = float(numpy.round(m_b[4],6))
            rcif_diffrn.diffrn_orient_matrix.ub_23 = float(numpy.round(m_b[5],6))
            rcif_diffrn.diffrn_orient_matrix.ub_31 = float(numpy.round(m_b[6],6))
            rcif_diffrn.diffrn_orient_matrix.ub_32 = float(numpy.round(m_b[7],6))
            rcif_diffrn.diffrn_orient_matrix.ub_33 = float(numpy.round(m_b[8],6))

            n_fractional = int(numpy.any(numpy.abs(np_hkl - numpy.round(np_hkl)) > 1e-4,
                                         axis=0).sum())
            text = (f"{rcif_diffrn.data_name}: {np_hkl.shape[1] - n_fractional} reflections with "
                    f"integer indices and {n_fractional} satellites of crystal {item.data_name} "
                    f"(wavelength {float(wavelength):g} Å, 2theta {ttheta_min:g}–{ttheta_max:g}°).")
            if flag_extra:
                text += (f" k = ({', '.join(f'{float(value):g}' for value in k_extra)}) is not a "
                         "propagation vector of the crystal: its satellites have zero "
                         "intensity (a list for planning the experiment); add_propagation_vector "
                         "gives them a magnetic structure.")
            ls_out.append(text)
            l_item_diffrn_refln = []
            for hkl in np_hkl.T:
                item_diffrn_refln = DiffrnRefln(
                    index_h=hkl[0], index_k=hkl[1], index_l=hkl[2], 
                    intensity=1., intensity_sigma=0.01, 
                    )
                l_item_diffrn_refln.append(item_diffrn_refln)
            loop_diffrn_refln = DiffrnReflnL()
            loop_diffrn_refln.items = l_item_diffrn_refln
            rcif_diffrn.diffrn_refln = loop_diffrn_refln
            l_diffrn.append(rcif_diffrn)

    cryspy_object.add_items(l_diffrn)
    rhochi_no_refinement(cryspy_object)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            item.intensity = float(numpy.round(item.intensity_calc, 5))
    rhochi_no_refinement(cryspy_object)
    return "\n".join(ls_out)


def simulation_fliping_ratio(cryspy_object: cryspy.GlobalN):
    """Simulate FR experiment with default parameters
    """
    f_diffrn_fr = os.path.join(DIR, "diffrn_fr.rcif")

    l_diffrn = []
    i_crystal = 0
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            i_crystal += 1
            d_crystal = item.get_dictionary()

            rcif_diffrn = cryspy.load_file(f_diffrn_fr).diffrn_fr
            if i_crystal == 1:
                rcif_diffrn.data_name = "fr"
            else:
                rcif_diffrn.data_name = f"fr{i_crystal}"
            rcif_diffrn.phase.label = item.data_name
            wavelength = rcif_diffrn.setup.wavelength

            unit_cell_parameters = d_crystal["unit_cell_parameters"]
            
            m_b = calc_m_b_by_unit_cell_parameters(unit_cell_parameters)[0]

            reduced_symm_elems = d_crystal.get("reduced_symm_elems", numpy.array([0,0,0,1,1,0,0,0,1,0,0,0,1],dtype=int))
            translation_elems = d_crystal.get("translation_elems", numpy.array([0,0,0,1],dtype=int))
            centrosymmetry = d_crystal.get("centrosymmetry", True)
            sthovl_min = 0.01/wavelength
            sthovl_max = 0.85/wavelength
            np_hkl = calc_index_hkl_multiplicity_in_range(
                sthovl_min=sthovl_min, 
                sthovl_max=sthovl_max, unit_cell_parameters=unit_cell_parameters, 
                reduced_symm_elems=reduced_symm_elems, translation_elems=translation_elems, 
                centrosymmetry=centrosymmetry)[0]
            
            rcif_diffrn.diffrn_orient_matrix.ub_11 = float(numpy.round(m_b[0],6))
            rcif_diffrn.diffrn_orient_matrix.ub_12 = float(numpy.round(m_b[1],6))
            rcif_diffrn.diffrn_orient_matrix.ub_13 = float(numpy.round(m_b[2],6))
            rcif_diffrn.diffrn_orient_matrix.ub_21 = float(numpy.round(m_b[3],6))
            rcif_diffrn.diffrn_orient_matrix.ub_22 = float(numpy.round(m_b[4],6))
            rcif_diffrn.diffrn_orient_matrix.ub_23 = float(numpy.round(m_b[5],6))
            rcif_diffrn.diffrn_orient_matrix.ub_31 = float(numpy.round(m_b[6],6))
            rcif_diffrn.diffrn_orient_matrix.ub_32 = float(numpy.round(m_b[7],6))
            rcif_diffrn.diffrn_orient_matrix.ub_33 = float(numpy.round(m_b[8],6))

            l_item_diffrn_refln = []
            for hkl in np_hkl.T:
                item_diffrn_refln = DiffrnRefln(
                    index_h=hkl[0], index_k=hkl[1], index_l=hkl[2], 
                    fr=1., fr_sigma=0.01, 
                    )
                l_item_diffrn_refln.append(item_diffrn_refln)
            loop_diffrn_refln = DiffrnReflnL()
            loop_diffrn_refln.items = l_item_diffrn_refln
            rcif_diffrn.diffrn_refln = loop_diffrn_refln
            l_diffrn.append(rcif_diffrn)

    cryspy_object.add_items(l_diffrn)
    remove_extinct_reflections(cryspy_object, l_diffrn)
    rhochi_no_refinement(cryspy_object)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            item.fr = float(numpy.round(item.fr_calc, 5))
    rhochi_no_refinement(cryspy_object)


# Relative intensity below which a reflection is extinct: the intensities of
# systematic absences are round-off, about 1e-30 of the strongest reflection.
EXTINCT_INTENSITY = 1e-10


def remove_extinct_reflections(cryspy_object: cryspy.GlobalN, l_diffrn: list):
    """Remove the reflections without nuclear and magnetic scattering from the
    simulated flipping ratio experiments.

    Their flipping ratio is the ratio of two round-off values: it changes with
    the last bits of the moments from one calculation to the next and cannot
    be measured.
    """
    dict_in_out = {}
    rhochi_calc_chi_sq_by_dictionary(
        cryspy_object.get_dictionary(), dict_in_out=dict_in_out,
        flag_use_precalculated_data=False, flag_calc_analytical_derivatives=False)
    for rcif_diffrn in l_diffrn:
        dict_in_out_diffrn = dict_in_out[rcif_diffrn.get_dictionary()["type_name"]]
        items = rcif_diffrn.diffrn_refln.items
        n_hkl = len(items)
        # twin domains are blocks of n_hkl reflections
        iint = (dict_in_out_diffrn["iint_plus"] + dict_in_out_diffrn["iint_minus"]).reshape(-1, n_hkl).sum(axis=0)
        flags_measurable = iint > EXTINCT_INTENSITY * iint.max()
        loop_diffrn_refln = DiffrnReflnL()
        loop_diffrn_refln.items = [item for item, flag in zip(items, flags_measurable) if flag]
        rcif_diffrn.diffrn_refln = loop_diffrn_refln


def simulation_current_experiments(cryspy_object: cryspy.GlobalN):
    """Simulate existing experiments
    """
    rhochi_no_refinement(cryspy_object)
    l_diffrn = []
    l_pd = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Diffrn):
            l_diffrn.append(item)
        elif isinstance(item, cryspy.Pd):
            l_pd.append(item)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            if item.is_attribute('fr_calc'):
                item.fr = float(numpy.round(item.fr_calc, 5))
            elif item.is_attribute('intensity_calc'):
                item.intensity = float(numpy.round(item.intensity_calc, 5))


    for rcif_pd in l_pd:
        np_tth = rcif_pd.pd_proc.numpy_ttheta
        np_plus_net = rcif_pd.pd_proc.numpy_intensity_plus_net
        np_minus_net = rcif_pd.pd_proc.numpy_intensity_minus_net
        np_bkg_calc = rcif_pd.pd_proc.numpy_intensity_bkg_calc
        np_sigma = numpy.ones_like(np_plus_net) * 0.1

        flag_polarized = rcif_pd.pd_meas.items[0].is_attribute("intensity_plus") and rcif_pd.pd_meas.items[0].is_attribute("intensity_minus")
        flag_unpolarized = rcif_pd.pd_meas.items[0].is_attribute("intensity")

        if flag_polarized:  
            l_item_meas = []
            for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
                item_meas = PdMeas(ttheta=tth, intensity_plus=plus + 0.5*bkg_calc, intensity_plus_sigma=sigma,intensity_minus=minus + 0.5*bkg_calc, intensity_minus_sigma=sigma,)
                l_item_meas.append(item_meas)
        elif flag_unpolarized:              
            l_item_meas = []
            for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
                item_meas = PdMeas(ttheta=tth, intensity=plus + minus + bkg_calc, intensity_sigma=sigma,)
                l_item_meas.append(item_meas)
        loop_meas = PdMeasL()
        loop_meas.items = l_item_meas
        rcif_pd.pd_meas = loop_meas
    rhochi_no_refinement(cryspy_object)
    pass
